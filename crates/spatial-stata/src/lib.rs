//! Stata SPI v3 plugin entry points backed by `spatial-core`.

mod spur_stata;
mod stata_spi;

use spatial_core::{run_model_scpc, CoreError, ErrorKind, Metric, ScpcOutput, ScpcPreparation};
use stata_spi::{
    active_obs, error, plugin_var_count, read_matrix, read_var_matrix, save_scalar, write_matrix,
    STATA_ERR, STATA_OK,
};
use std::ffi::CStr;
use std::os::raw::{c_char, c_int};
use std::panic::{catch_unwind, AssertUnwindSafe};

fn error_kind_label(kind: ErrorKind) -> &'static str {
    match kind {
        ErrorKind::InvalidInput => "invalid input",
        ErrorKind::UnsupportedCase => "unsupported case",
        ErrorKind::TiedEigenspace => "tied eigenspace",
        ErrorKind::InvalidSampledCovariance => "invalid sampled covariance",
        ErrorKind::CalibrationFailed => "calibration failed",
        ErrorKind::DegenerateVariance => "degenerate variance",
        ErrorKind::InternalFailure => "internal failure",
    }
}

pub(crate) fn format_core_error(err: CoreError) -> String {
    format!("{}: {}", error_kind_label(err.kind()), err.message())
}

pub(crate) fn parse_bool_flag(raw: &str, name: &str) -> Result<bool, String> {
    match raw {
        "0" => Ok(false),
        "1" => Ok(true),
        _ => Err(format!(
            "invalid boolean flag {name}='{raw}' (expected 0 or 1)"
        )),
    }
}

pub(crate) fn parse_usize(raw: &str, name: &str) -> Result<usize, String> {
    raw.parse::<usize>()
        .map_err(|_| format!("invalid integer argument {name}='{raw}'"))
}

fn parse_seed(raw: &str) -> Result<u32, String> {
    let val = raw
        .parse::<u64>()
        .map_err(|_| format!("invalid seed argument '{raw}' (expected integer in [0, 2^32))"))?;
    u32::try_from(val).map_err(|_| format!("seed argument '{raw}' must be in [0, 2^32)"))
}

pub(crate) fn parse_nthreads(raw: &str) -> Result<Option<usize>, String> {
    let val = raw
        .parse::<isize>()
        .map_err(|_| format!("invalid nthreads argument '{raw}'"))?;
    if val == -1 {
        Ok(None)
    } else if val >= 1 {
        let nthreads = Some(val as usize);
        spatial_core::threads::validate_nthreads(nthreads).map_err(format_core_error)?;
        Ok(nthreads)
    } else {
        Err(format!("nthreads must be -1 or >= 1, got {val}"))
    }
}

/// Expand the Stata design (endogenous columns first, constant last) into the
/// fitted structural and instrument columns required by the shared operation.
struct StataDesign {
    x: Vec<f64>,
    z: Option<Vec<f64>>,
    z_cols: usize,
}

impl StataDesign {
    fn from_columns(
        rows: usize,
        p: usize,
        has_cons: bool,
        n_endog: usize,
        z_cols: usize,
        raw_x: &[f64],
        excluded: &[f64],
    ) -> Result<Self, String> {
        let nx = p
            .checked_sub(usize::from(has_cons))
            .ok_or("invalid Stata regressor count")?;
        if n_endog > nx
            || raw_x.len() != rows.checked_mul(nx).ok_or("Stata design too large")?
            || excluded.len()
                != rows
                    .checked_mul(z_cols)
                    .ok_or("Stata instrument design too large")?
        {
            return Err("invalid Stata model columns".into());
        }
        let instruments = if n_endog == 0 {
            0
        } else {
            p - n_endog + z_cols
        };
        let mut x = Vec::with_capacity(rows.checked_mul(p).ok_or("Stata design too large")?);
        let mut z = if n_endog == 0 {
            None
        } else {
            Some(Vec::with_capacity(
                rows.checked_mul(instruments)
                    .ok_or("Stata instrument design too large")?,
            ))
        };
        for i in 0..rows {
            let xrow = &raw_x[i * nx..(i + 1) * nx];
            x.extend_from_slice(xrow);
            if has_cons {
                x.push(1.0);
            }
            if let Some(ref mut values) = z {
                values.extend_from_slice(&xrow[n_endog..]);
                values.extend_from_slice(&excluded[i * z_cols..(i + 1) * z_cols]);
                if has_cons {
                    values.push(1.0);
                }
            }
        }
        Ok(Self {
            x,
            z,
            z_cols: instruments,
        })
    }
}

fn cluster_coordinates(
    cl_raw: &[f64],
    coords_obs: &[f64],
    d: usize,
) -> Result<(usize, Vec<usize>, Vec<f64>), String> {
    let m_obs = cl_raw.len();
    if m_obs < 2 {
        return Err(format!(
            "too few active observations ({m_obs}) for clustered SCPC"
        ));
    }
    let mut cluster_id = Vec::with_capacity(m_obs);
    let mut max_id = 0usize;
    for &val in cl_raw {
        if !val.is_finite() || val < 1.0 || val > m_obs as f64 || val != val.round() {
            return Err(format!("invalid cluster index '{val}' in plugin varlist"));
        }
        let g = (val as usize) - 1;
        max_id = max_id.max(g);
        cluster_id.push(g);
    }
    let n_clust = max_id + 1;
    if n_clust < 2 {
        return Err(format!("too few active clusters ({n_clust}) for SCPC"));
    }

    let mut seen = vec![false; n_clust];
    let n_coords = n_clust
        .checked_mul(d)
        .ok_or("cluster coordinate size overflow")?;
    let mut coords_cl = vec![0.0; n_coords];
    for (i, &g) in cluster_id.iter().enumerate() {
        let obs_c = &coords_obs[i * d..(i + 1) * d];
        let cl_c = &mut coords_cl[g * d..(g + 1) * d];
        if !seen[g] {
            seen[g] = true;
            cl_c.copy_from_slice(obs_c);
        } else if obs_c != cl_c {
            return Err("coordinates are not constant within clusters".to_string());
        }
    }
    if seen.iter().any(|&s| !s) {
        return Err("cluster index vector contains empty clusters".to_string());
    }
    Ok((n_clust, cluster_id, coords_cl))
}

struct ModelScpcArgs<'a> {
    d: usize,
    p: usize,
    k_use: usize,
    has_cons: bool,
    avc: f64,
    latlong: bool,
    uncond: bool,
    cvs: bool,
    nthreads: Option<usize>,
    b_mat_name: &'a str,
    bread_mat_name: &'a str,
    stats_mat_name: &'a str,
    cvs_mat_name: &'a str,
    c0_scal_name: &'a str,
    cv_scal_name: &'a str,
    q_scal_name: &'a str,
    tied_scal_name: &'a str,
    has_cluster: bool,
    n_endog: usize,
    n_z_excl: usize,
}

impl<'a> ModelScpcArgs<'a> {
    fn parse(args: &'a [&'a str]) -> Result<Self, String> {
        if args.len() != 20 {
            return Err(format!(
                "model_exact_scpc expects 20 arguments after subcommand, got {}",
                args.len()
            ));
        }

        let d = parse_usize(args[0], "d")?;
        let p = parse_usize(args[1], "p")?;
        let k_use = parse_usize(args[2], "k_use")?;
        let has_cons = parse_bool_flag(args[3], "has_cons")?;
        let avc = args[4]
            .parse::<f64>()
            .map_err(|_| format!("invalid avc '{}'", args[4]))?;
        let latlong = parse_bool_flag(args[5], "latlong")?;
        let uncond = parse_bool_flag(args[6], "uncond")?;
        let cvs = parse_bool_flag(args[7], "cvs")?;
        let nthreads = parse_nthreads(args[8])?;

        let b_mat_name = args[9];
        let bread_mat_name = args[10];
        let stats_mat_name = args[11];
        let cvs_mat_name = args[12];
        let c0_scal_name = args[13];
        let cv_scal_name = args[14];
        let q_scal_name = args[15];
        let tied_scal_name = args[16];
        let has_cluster = parse_bool_flag(args[17], "has_cluster")?;
        let n_endog = parse_usize(args[18], "n_endog")?;
        let n_z_excl = parse_usize(args[19], "n_z_excl")?;

        if d == 0 || p == 0 || k_use == 0 || k_use > p {
            return Err(format!("invalid dimensions: d={d}, p={p}, k_use={k_use}"));
        }
        let nx = p - usize::from(has_cons);
        if n_endog > nx {
            return Err(format!("invalid n_endog={n_endog} > nx={nx}"));
        }
        if latlong && d != 2 {
            return Err(format!("latlong requires d=2 coordinates, got d={d}"));
        }
        if !avc.is_finite() || avc <= 0.001 || avc >= 0.99 {
            return Err(format!("avc must lie strictly in (0.001, 0.99), got {avc}"));
        }
        if b_mat_name.is_empty()
            || bread_mat_name.is_empty()
            || stats_mat_name.is_empty()
            || cvs_mat_name.is_empty()
            || c0_scal_name.is_empty()
            || cv_scal_name.is_empty()
            || q_scal_name.is_empty()
            || tied_scal_name.is_empty()
        {
            return Err("SCPC model matrix and scalar names must be non-empty".to_string());
        }

        Ok(Self {
            d,
            p,
            k_use,
            has_cons,
            avc,
            latlong,
            uncond,
            cvs,
            nthreads,
            b_mat_name,
            bread_mat_name,
            stats_mat_name,
            cvs_mat_name,
            c0_scal_name,
            cv_scal_name,
            q_scal_name,
            tied_scal_name,
            has_cluster,
            n_endog,
            n_z_excl,
        })
    }
}

struct LoadedStataScpc {
    n_sites: usize,
    metric: Metric,
    coords: Vec<f64>,
    b: Vec<f64>,
    cluster_id: Option<Vec<usize>>,
    n_obs: usize,
    scores: Vec<f64>,
    bread: Vec<f64>,
    design: Option<StataDesign>,
}

impl LoadedStataScpc {
    fn prepare_request(&self, k: usize) -> ScpcPreparation<'_> {
        let p = self.b.len();
        ScpcPreparation {
            n_sites: self.n_sites,
            n_obs: self.n_obs,
            n_scores: p,
            coefficients: &self.b[..k],
            bread: &self.bread[..k * p],
            scores: &self.scores,
            design: self.design.as_ref().map(|v| v.x.as_slice()),
            design_cols: self.design.as_ref().map_or(0, |_| p),
            instruments: self.design.as_ref().and_then(|v| v.z.as_deref()),
            instrument_cols: self.design.as_ref().map_or(0, |v| v.z_cols),
            direction_design: None,
            cluster_ids: self.cluster_id.as_deref(),
            fixed_effects: None,
        }
    }
}

fn load_stata_scpc_inputs(parsed: &ModelScpcArgs<'_>) -> Result<LoadedStataScpc, String> {
    let nx = parsed.p - usize::from(parsed.has_cons);
    let cl_vars = usize::from(parsed.has_cluster);
    let iv_vars = if !parsed.uncond && parsed.n_endog > 0 {
        parsed.n_z_excl
    } else {
        0
    };
    let expected_vars =
        parsed.d + parsed.p + cl_vars + if parsed.uncond { 0 } else { nx + iv_vars };
    let actual_vars = plugin_var_count();
    if actual_vars != expected_vars {
        return Err(format!(
            "scpc plugin expected {expected_vars} variables in plugin varlist, got {actual_vars}"
        ));
    }

    let obs = active_obs();
    let m_obs = obs.len();
    if m_obs < 2 {
        return Err(format!("too few active observations ({m_obs}) for SCPC"));
    }

    // Read all Stata dataset variables and matrices on the calling thread before Rayon.
    let coords_obs = read_var_matrix(&obs, 1, parsed.d, "coordinates (s_*)")?;
    let scores = read_var_matrix(&obs, 1 + parsed.d, parsed.p, "scores (scpc_score*)")?;
    let cl_raw = if parsed.has_cluster {
        read_var_matrix(&obs, 1 + parsed.d + parsed.p, 1, "cluster index (cl_idx)")?
    } else {
        Vec::new()
    };
    let x_raw = if !parsed.uncond && nx > 0 {
        read_var_matrix(
            &obs,
            1 + parsed.d + parsed.p + cl_vars,
            nx,
            "conditional regressors",
        )?
    } else {
        Vec::new()
    };
    let z_excl = if !parsed.uncond && parsed.n_endog > 0 && parsed.n_z_excl > 0 {
        read_var_matrix(
            &obs,
            1 + parsed.d + parsed.p + cl_vars + nx,
            parsed.n_z_excl,
            "excluded IV instruments",
        )?
    } else {
        Vec::new()
    };

    let (n_sites, cluster_id_vec, coords) = if parsed.has_cluster {
        let (sites, ids, coords) = cluster_coordinates(&cl_raw, &coords_obs, parsed.d)?;
        (sites, Some(ids), coords)
    } else {
        (m_obs, None, coords_obs)
    };

    let b = read_matrix(parsed.b_mat_name, 1, parsed.p)?;
    let bread = read_matrix(parsed.bread_mat_name, parsed.p, parsed.p)?;

    let design = if parsed.uncond {
        None
    } else {
        Some(StataDesign::from_columns(
            m_obs,
            parsed.p,
            parsed.has_cons,
            parsed.n_endog,
            parsed.n_z_excl,
            &x_raw,
            &z_excl,
        )?)
    };
    let metric = if parsed.latlong {
        Metric::GreatCircleFractionLatLon
    } else {
        Metric::Euclidean
    };

    Ok(LoadedStataScpc {
        n_sites,
        metric,
        coords,
        b,
        cluster_id: cluster_id_vec,
        n_obs: m_obs,
        scores,
        bread,
        design,
    })
}

fn write_stata_scpc_outputs(parsed: &ModelScpcArgs<'_>, out: &ScpcOutput) -> Result<(), String> {
    write_matrix(parsed.stats_mat_name, parsed.k_use, 6, &out.scpcstats)?;
    if parsed.cvs {
        let cvs_mat = out
            .scpccvs
            .as_deref()
            .ok_or_else(|| "missing scpccvs output when cvs=1".to_string())?;
        write_matrix(parsed.cvs_mat_name, parsed.k_use, 4, cvs_mat)?;
    }
    save_scalar(parsed.c0_scal_name, out.c0)?;
    save_scalar(parsed.cv_scal_name, out.cvfin)?;
    save_scalar(parsed.q_scal_name, out.q as f64)?;
    save_scalar(parsed.tied_scal_name, if out.tied { 1.0 } else { 0.0 })
}

fn handle_model_scpc(args: &[&str], approximate: bool) -> Result<(), String> {
    let (parsed, seed) = if approximate {
        if args.len() != 21 {
            return Err(format!(
                "model_approx_scpc expects 21 arguments after subcommand, got {}",
                args.len()
            ));
        }
        let seed = parse_seed(args[20])?;
        (ModelScpcArgs::parse(&args[..20])?, Some(seed))
    } else {
        (ModelScpcArgs::parse(args)?, None)
    };
    let loaded = load_stata_scpc_inputs(&parsed)?;
    match run_model_scpc(
        loaded.prepare_request(parsed.k_use),
        &loaded.coords,
        parsed.d,
        loaded.metric,
        parsed.avc,
        parsed.cvs,
        parsed.nthreads,
        seed,
    ) {
        Ok(out) => write_stata_scpc_outputs(&parsed, &out),
        Err(err) if approximate && err.kind() == ErrorKind::TiedEigenspace => {
            save_scalar(parsed.tied_scal_name, 1.0)
        }
        Err(err) => Err(format_core_error(err)),
    }
}

fn handle_build_info(args: &[&str]) -> Result<(), String> {
    if !args.is_empty() {
        return Err("build_info takes no arguments".to_string());
    }
    stata_spi::display(&format!(
        "spatial-core {} (spec {}, {}, {}-{}, SCPC API {}, SPUR API {})",
        spatial_core::VERSION,
        spatial_core::SPEC_VERSION,
        spatial_core::linalg_provider(),
        std::env::consts::ARCH,
        std::env::consts::OS,
        spatial_core::SCPC_API_VERSION,
        spatial_core::SPUR_API_VERSION,
    ));
    Ok(())
}

fn dispatch(args: &[&str]) -> Result<(), String> {
    let subcmd = *args
        .first()
        .ok_or_else(|| "missing plugin subcommand".to_string())?;
    match subcmd {
        "ping" | "ping_v7" | "ping_spur_v4" => Ok(()),
        "build_info" => handle_build_info(&args[1..]),
        "model_exact_scpc" => handle_model_scpc(&args[1..], false),
        "model_approx_scpc" => handle_model_scpc(&args[1..], true),
        "spur_i0" => spur_stata::handle_spur_test(&args[1..], spur_stata::SpurTestFlavor::I0),
        "spur_i1" => spur_stata::handle_spur_test(&args[1..], spur_stata::SpurTestFlavor::I1),
        "spur_i0_resid" => {
            spur_stata::handle_spur_test(&args[1..], spur_stata::SpurTestFlavor::I0Resid)
        }
        "spur_i1_resid" => {
            spur_stata::handle_spur_test(&args[1..], spur_stata::SpurTestFlavor::I1Resid)
        }
        "spur_halflife" => spur_stata::handle_spur_halflife(&args[1..]),
        "spur_transform" => spur_stata::handle_spur_transform(&args[1..]),
        "spur_pipeline" => spur_stata::handle_spur_pipeline(&args[1..]),
        other => Err(format!("unknown spatial-stata subcommand '{other}'")),
    }
}

/// Entry point invoked from `c/stata_bridge.c`.
///
/// # Safety
/// `argv` must point to at least `argc` valid NUL-terminated C string pointers when `argc > 0`.
#[no_mangle]
pub unsafe extern "C" fn rs_stata_call(argc: c_int, argv: *const *const c_char) -> c_int {
    let res = catch_unwind(AssertUnwindSafe(|| {
        if argc <= 0 || argv.is_null() {
            return Err("plugin called without subcommand arguments".to_string());
        }
        let mut args = Vec::with_capacity(argc as usize);
        for idx in 0..(argc as usize) {
            let ptr = unsafe { *argv.add(idx) };
            if ptr.is_null() {
                return Err(format!("null argv[{idx}] pointer"));
            }
            let arg = unsafe { CStr::from_ptr(ptr) }
                .to_str()
                .map_err(|_| format!("argv[{idx}] is not valid UTF-8"))?;
            args.push(arg);
        }
        dispatch(&args)
    }));

    match res {
        Ok(Ok(())) => STATA_OK,
        Ok(Err(msg)) => {
            error(&format!("spatial Rust plugin error: {msg}"));
            STATA_ERR
        }
        Err(_) => {
            error("spatial Rust plugin panicked unexpectedly");
            STATA_ERR
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spatial_core::preparation::{prepare_scpc_inputs, PreparedScpcInputs};

    #[allow(clippy::too_many_arguments)]
    fn fixture_model(
        n_sites: usize,
        n_obs: usize,
        has_cons: bool,
        n_endog: usize,
        z_cols: usize,
        b: &[f64],
        bread: &[f64],
        scores: &[f64],
        raw_x: &[f64],
        excluded: &[f64],
        ids: Option<&[usize]>,
    ) -> Result<LoadedStataScpc, String> {
        let design =
            StataDesign::from_columns(n_obs, b.len(), has_cons, n_endog, z_cols, raw_x, excluded)?;
        Ok(LoadedStataScpc {
            n_sites,
            n_obs,
            metric: Metric::Euclidean,
            coords: Vec::new(),
            b: b.to_vec(),
            bread: bread.to_vec(),
            scores: scores.to_vec(),
            design: Some(design),
            cluster_id: ids.map(<[usize]>::to_vec),
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn fixture_inputs(
        n_sites: usize,
        n_obs: usize,
        has_cons: bool,
        n_endog: usize,
        z_cols: usize,
        b: &[f64],
        bread: &[f64],
        scores: &[f64],
        raw_x: &[f64],
        excluded: &[f64],
        ids: Option<&[usize]>,
    ) -> Result<PreparedScpcInputs, String> {
        let model = fixture_model(
            n_sites, n_obs, has_cons, n_endog, z_cols, b, bread, scores, raw_x, excluded, ids,
        )?;
        prepare_scpc_inputs(model.prepare_request(b.len())).map_err(format_core_error)
    }

    fn ols_fixture_basis(
        x: &[f64],
        n: usize,
        nx: usize,
        has_cons: bool,
    ) -> Result<(Vec<f64>, usize), String> {
        let p = nx + usize::from(has_cons);
        let res = fixture_inputs(
            n,
            n,
            has_cons,
            0,
            0,
            &vec![0.; p],
            &vec![0.; p * p],
            &vec![0.; n * p],
            x,
            &[],
            None,
        )?;
        Ok((res.proj_u.unwrap(), res.proj_rank))
    }

    fn load_bench_200_csv() -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let path = format!(
            "{}/../spatial-core/tests/fixtures/bench_200.csv",
            env!("CARGO_MANIFEST_DIR")
        );
        let text = std::fs::read_to_string(path).expect("bench_200.csv");
        let mut y = Vec::new();
        let mut x1 = Vec::new();
        let mut x2 = Vec::new();
        let mut coords_lon_lat = Vec::new();
        for line in text.lines().skip(1) {
            if line.trim().is_empty() {
                continue;
            }
            let cols: Vec<f64> = line
                .split(',')
                .map(|s| s.trim().parse::<f64>().unwrap())
                .collect();
            y.push(cols[1]);
            x1.push(cols[2]);
            x2.push(cols[3]);
            coords_lon_lat.extend_from_slice(&[cols[6], cols[7]]);
        }
        (y, x1, x2, coords_lon_lat)
    }

    fn invert_symmetric(mat: &[f64], p: usize) -> Vec<f64> {
        let mut a = mat.to_vec();
        let mut inv = vec![0.0; p * p];
        for i in 0..p {
            inv[i * p + i] = 1.0;
        }
        for col in 0..p {
            let mut pivot = col;
            for row in (col + 1)..p {
                if a[row * p + col].abs() > a[pivot * p + col].abs() {
                    pivot = row;
                }
            }
            if pivot != col {
                for k in 0..p {
                    a.swap(col * p + k, pivot * p + k);
                    inv.swap(col * p + k, pivot * p + k);
                }
            }
            let diag = a[col * p + col];
            for k in 0..p {
                a[col * p + k] /= diag;
                inv[col * p + k] /= diag;
            }
            for row in 0..p {
                if row == col {
                    continue;
                }
                let factor = a[row * p + col];
                for k in 0..p {
                    a[row * p + k] -= factor * a[col * p + k];
                    inv[row * p + k] -= factor * inv[col * p + k];
                }
            }
        }
        inv
    }

    #[test]
    fn orthonormal_basis_handles_cons_and_nocons() {
        let n = 4;
        let x_raw = vec![1.0, 2.0, 3.0, 4.0];
        let (q_cons, rank_cons) = ols_fixture_basis(&x_raw, n, 1, true).unwrap();
        assert_eq!(rank_cons, 2);
        // Check Q' Q == I_2
        for r1 in 0..2 {
            for r2 in 0..2 {
                let dot: f64 = (0..n)
                    .map(|i| q_cons[i * 2 + r1] * q_cons[i * 2 + r2])
                    .sum();
                let exp = if r1 == r2 { 1.0 } else { 0.0 };
                assert!((dot - exp).abs() < 1e-14);
            }
        }
        let (q_nocons, rank_nocons) = ols_fixture_basis(&x_raw, n, 1, false).unwrap();
        assert_eq!(rank_nocons, 1);
        let norm2: f64 = q_nocons.iter().map(|v| v * v).sum();
        assert!((norm2 - 1.0).abs() < 1e-14);

        // Intercept-only (nx = 0, has_cons = true)
        let (q_int, rank_int) = ols_fixture_basis(&[], n, 0, true).unwrap();
        assert_eq!(rank_int, 1);
        for val in q_int {
            assert!((val.abs() - 0.5).abs() < 1e-14);
        }
    }

    #[test]
    fn ols_preparation_matches_bench_200_euclid_and_latlong() {
        let (y, x1, x2, coords_lon_lat) = load_bench_200_csv();
        let n = y.len();
        let p = 3;
        let nx = 2;
        let mut x_raw = Vec::with_capacity(n * nx);
        let mut x_design = Vec::with_capacity(n * p);
        for i in 0..n {
            x_raw.extend_from_slice(&[x1[i], x2[i]]);
            x_design.extend_from_slice(&[x1[i], x2[i], 1.0]);
        }

        let mut xtx = vec![0.0; p * p];
        let mut xty = vec![0.0; p];
        for i in 0..n {
            for r in 0..p {
                xty[r] += x_design[i * p + r] * y[i];
                for c in 0..p {
                    xtx[r * p + c] += x_design[i * p + r] * x_design[i * p + c];
                }
            }
        }
        let bread = invert_symmetric(&xtx, p);
        let mut b = vec![0.0; p];
        for r in 0..p {
            b[r] = (0..p).map(|c| bread[r * p + c] * xty[c]).sum();
        }
        let mut scores = vec![0.0; n * p];
        for i in 0..n {
            let yhat: f64 = (0..p).map(|c| x_design[i * p + c] * b[c]).sum();
            let resid = y[i] - yhat;
            for c in 0..p {
                scores[i * p + c] = x_design[i * p + c] * resid;
            }
        }

        let model =
            fixture_model(n, n, true, 0, 0, &b, &bread, &scores, &x_raw, &[], None).unwrap();
        let out_euclid = run_model_scpc(
            model.prepare_request(p),
            &coords_lon_lat,
            2,
            Metric::Euclidean,
            0.03,
            true,
            Some(2),
            None,
        )
        .unwrap();
        assert_eq!(out_euclid.q, 19);
        let dist_euclid =
            spatial_core::distance::distance_matrix(&coords_lon_lat, n, 2, Metric::Euclidean)
                .unwrap();
        let mut dvec_euclid = Vec::with_capacity(n * (n - 1) / 2);
        for i in 1..n {
            dvec_euclid.extend_from_slice(&dist_euclid[i * n..i * n + i]);
        }
        assert!(
            (spatial_core::calibration::avg_correlation(out_euclid.c0, &dvec_euclid) - 0.03).abs()
                < 1e-12
        );
        assert!((out_euclid.c0 - 3.500409719601247).abs() < 1e-3);
        assert!((out_euclid.cvfin - 2.1008365544083087).abs() < 1e-3);
        assert!(!out_euclid.tied);
        // Row 0 is x1, row 1 is x2, row 2 is _cons
        assert!((out_euclid.scpcstats[0] - 1.5473280763431876).abs() < 1e-10);
        assert!((out_euclid.scpcstats[12] - 1.4834167407175964).abs() < 1e-10);

        // Stata latlong: s_1=lat, s_2=lon
        let mut coords_lat_lon = Vec::with_capacity(n * 2);
        for pair in coords_lon_lat.chunks_exact(2) {
            coords_lat_lon.extend_from_slice(&[pair[1], pair[0]]);
        }
        let out_geo = run_model_scpc(
            model.prepare_request(p),
            &coords_lat_lon,
            2,
            Metric::GreatCircleFractionLatLon,
            0.03,
            true,
            Some(2),
            None,
        )
        .unwrap();
        assert_eq!(out_geo.q, 18);
        let dist_geo = spatial_core::distance::distance_matrix(
            &coords_lat_lon,
            n,
            2,
            Metric::GreatCircleFractionLatLon,
        )
        .unwrap();
        let mut dvec_geo = Vec::with_capacity(n * (n - 1) / 2);
        for i in 1..n {
            dvec_geo.extend_from_slice(&dist_geo[i * n..i * n + i]);
        }
        assert!(
            (spatial_core::calibration::avg_correlation(out_geo.c0, &dvec_geo) - 0.03).abs()
                < 1e-12
        );
        assert!(
            (out_geo.c0 - 1396.4193247376622).abs() < 1.0,
            "geo c0={}",
            out_geo.c0
        );
        assert!(
            (out_geo.cvfin - 2.1009220402410382).abs() < 5e-3,
            "geo cvfin={}",
            out_geo.cvfin
        );
        assert!(!out_geo.tied);

        // First 20 rows of bench_200 exhibit tied eigenspace for automatic fallback.
        let n_tie = 20;
        let mut xtx_20 = vec![0.0; p * p];
        let mut xty_20 = vec![0.0; p];
        for i in 0..n_tie {
            for r in 0..p {
                xty_20[r] += x_design[i * p + r] * y[i];
                for c in 0..p {
                    xtx_20[r * p + c] += x_design[i * p + r] * x_design[i * p + c];
                }
            }
        }
        let bread_20 = invert_symmetric(&xtx_20, p);
        let mut b_20 = vec![0.0; p];
        for r in 0..p {
            b_20[r] = (0..p).map(|c| bread_20[r * p + c] * xty_20[c]).sum();
        }
        let mut scores_20 = vec![0.0; n_tie * p];
        for i in 0..n_tie {
            let yhat: f64 = (0..p).map(|c| x_design[i * p + c] * b_20[c]).sum();
            let resid = y[i] - yhat;
            for c in 0..p {
                scores_20[i * p + c] = x_design[i * p + c] * resid;
            }
        }
        let model_20 = fixture_model(
            n_tie,
            n_tie,
            true,
            0,
            0,
            &b_20,
            &bread_20,
            &scores_20,
            &x_raw[..n_tie * nx],
            &[],
            None,
        )
        .unwrap();
        let out_20 = run_model_scpc(
            model_20.prepare_request(p),
            &coords_lon_lat[..n_tie * 2],
            2,
            Metric::Euclidean,
            0.03,
            false,
            Some(1),
            None,
        )
        .unwrap();
        assert_eq!(out_20.q, 7);
        assert!(out_20.tied);
    }

    #[test]
    fn ols_preparation_is_scale_invariant_and_rejects_collinear_selected_columns() {
        let n = 5;
        let base = [1.0, -2.0, 0.5, 3.0, -1.5];
        // Small-unit regressor 1e-12 * x and extreme tiny 1e-200 * x stay full rank with _cons.
        for scale in [1e-12, 1e-200, 1e12] {
            let x_raw: Vec<f64> = base.iter().map(|&v| v * scale).collect();
            let (q, rank) = ols_fixture_basis(&x_raw, n, 1, true).unwrap();
            assert_eq!(rank, 2, "scale={scale}");
            for r1 in 0..2 {
                for r2 in 0..2 {
                    let dot: f64 = (0..n).map(|i| q[i * 2 + r1] * q[i * 2 + r2]).sum();
                    let exp = if r1 == r2 { 1.0 } else { 0.0 };
                    assert!((dot - exp).abs() < 1e-14, "scale={scale}, dot={dot}");
                }
            }
        }

        let mut x_raw = Vec::new();
        for &v in &base {
            x_raw.extend_from_slice(&[1e-12 * v, 2e-12 * v, 0.0]);
        }
        assert!(ols_fixture_basis(&x_raw, n, 3, true)
            .unwrap_err()
            .contains("rank deficient"));
    }

    #[test]
    fn dispatch_supports_only_the_model_scpc_abi() {
        for cmd in ["ping", "ping_v7"] {
            assert!(dispatch(&[cmd]).is_ok());
        }
        for cmd in [
            "ping_v2",
            "ping_v3",
            "ping_v4",
            "ping_v5",
            "ping_v6",
            "exact_scpc",
            "approx_scpc",
            "unknown",
        ] {
            assert!(dispatch(&[cmd])
                .unwrap_err()
                .contains("unknown spatial-stata subcommand"));
        }
        assert!(dispatch(&[]).is_err());

        let args = [
            "2", "3", "3", "1", "0.03", "0", "0", "1", "-1", "b", "bread", "stats", "cvs", "c0",
            "cv", "q", "tied", "1", "0", "0",
        ];
        let parsed = ModelScpcArgs::parse(&args).unwrap();
        assert_eq!(parsed.tied_scal_name, "tied");
        assert!(parsed.has_cluster);

        // Stubs expose no variables. A valid request reaches the SPI varlist check.
        let err = handle_model_scpc(&args, false).unwrap_err();
        assert!(err.contains("variables in plugin varlist"), "err={err}");
        assert!(ModelScpcArgs::parse(&args[..19]).is_err());
        for (position, value) in [(16, ""), (17, "extra"), (4, "0.99")] {
            let mut bad = args;
            bad[position] = value;
            assert!(ModelScpcArgs::parse(&bad).is_err());
        }
        let mut bad_latlong = args;
        bad_latlong[0] = "1";
        bad_latlong[5] = "1";
        assert!(ModelScpcArgs::parse(&bad_latlong).is_err());
    }

    #[test]
    fn error_mapping_and_c_entrypoint_behave_safely() {
        assert_eq!(
            format_core_error(CoreError::calibration("bad target")),
            "calibration failed: bad target"
        );
        assert_eq!(
            format_core_error(CoreError::degenerate("zero var")),
            "degenerate variance: zero var"
        );
        assert_eq!(
            format_core_error(CoreError::tied("tie")),
            "tied eigenspace: tie"
        );
        assert_eq!(
            format_core_error(CoreError::unsupported("unsupported")),
            "unsupported case: unsupported"
        );
        assert_eq!(
            format_core_error(CoreError::invalid_input("invalid")),
            "invalid input: invalid"
        );
        assert_eq!(
            format_core_error(CoreError::invalid_sampled_covariance("not psd")),
            "invalid sampled covariance: not psd"
        );
        assert_eq!(
            format_core_error(CoreError::internal("boom")),
            "internal failure: boom"
        );

        let ping_cstr = std::ffi::CString::new("ping_v7").unwrap();
        let argv_ok = [ping_cstr.as_ptr()];
        assert_eq!(unsafe { rs_stata_call(1, argv_ok.as_ptr()) }, STATA_OK);
        assert_eq!(unsafe { rs_stata_call(0, std::ptr::null()) }, STATA_ERR);
        let argv_null = [std::ptr::null()];
        assert_eq!(unsafe { rs_stata_call(1, argv_null.as_ptr()) }, STATA_ERR);
    }

    #[test]
    fn arg_parsers_validate_inputs() {
        assert_eq!(parse_nthreads("-1").unwrap(), None);
        assert_eq!(parse_nthreads("1").unwrap(), Some(1));
        assert_eq!(parse_nthreads("4").unwrap(), Some(4));
        assert!(parse_nthreads("0").is_err());
        assert!(parse_nthreads("-2").is_err());
        assert!(parse_nthreads("1025").is_err());
        assert!(!parse_bool_flag("0", "flag").unwrap());
        assert!(parse_bool_flag("1", "flag").unwrap());
        assert!(parse_bool_flag("2", "flag").is_err());
    }

    #[test]
    fn iv_preparation_matches_direct_2sls_projection() {
        let m_obs = 8;
        let p = 3; // [x_endog, w_exog, _cons]
        let nx = 2;
        let mut x_raw = Vec::with_capacity(m_obs * nx);
        let mut z_excl = Vec::with_capacity(m_obs);
        let mut y = Vec::with_capacity(m_obs);
        for i in 0..m_obs {
            let fi = i as f64;
            let wi = (fi * 0.7).sin();
            let zi = (fi * 1.1).cos() + 0.2 * fi;
            let xi = 0.8 * zi + 0.3 * wi + 0.25 * (fi * 1.9).sin();
            let yi = 1.0 + 1.4 * xi + 0.5 * wi + 0.3 * (fi * 1.9).sin();
            x_raw.extend_from_slice(&[xi, wi]);
            z_excl.push(zi);
            y.push(yi);
        }
        let b = [1.4, 0.5, 1.0];
        let mut bread = vec![0.0; p * p];
        for i in 0..p {
            bread[i * p + i] = 1.0;
        }
        let mut scores = vec![0.0; m_obs * p];
        for i in 0..m_obs {
            scores[i * p] = z_excl[i] * 0.1;
            scores[i * p + 1] = x_raw[i * nx + 1] * 0.1;
            scores[i * p + 2] = 0.1;
        }
        let inputs = fixture_inputs(
            m_obs, m_obs, true, 1, 1, &b, &bread, &scores, &x_raw, &z_excl, None,
        )
        .unwrap();
        let u = inputs.proj_u.as_ref().unwrap();
        let v = inputs.proj_v.as_ref().unwrap();
        let xj = inputs.xj.as_ref().unwrap();
        let xjs = inputs.xjs.as_ref().unwrap();
        assert_eq!(inputs.proj_rank, p);
        // Verify the normalized factor orientation M' = I - U V' with columns of V scaled:
        // V' * U == I_p and each xjs[i] == sign(xj[i]).
        for r in 0..p {
            for c in 0..p {
                let dot: f64 = (0..m_obs).map(|i| v[i * p + r] * u[i * p + c]).sum();
                let exp = if r == c { 1.0 } else { 0.0 };
                assert!((dot - exp).abs() < 1e-11, "V'U[{r},{c}]={dot}");
            }
        }
        // Check M' projects onto P_Z X's nullspace via column orthogonality: for any A, X' (A - U V' A) == 0.
        for i in 0..m_obs {
            for j in 0..p {
                let expected_sign = if xj[i * p + j] < 0.0 { -1.0 } else { 1.0 };
                assert_eq!(xjs[i * p + j], expected_sign);
            }
        }
        let expected_scales = [0, 1].map(|col| {
            (0..m_obs)
                .map(|i| x_raw[i * nx + col].abs())
                .fold(0.0_f64, f64::max)
        });
        for i in 0..m_obs {
            assert_eq!(v[i * p], x_raw[i * nx] / expected_scales[0]);
            assert_eq!(v[i * p + 1], x_raw[i * nx + 1] / expected_scales[1]);
            assert_eq!(v[i * p + 2], 1.0);
        }
    }

    #[test]
    fn iv_svd_avoids_normal_equations_and_rejects_irrelevant_instruments() {
        // Case 1: kappa(Q_Z' X) ~ 2e7, so kappa^2 ~ 4e14 would destroy normal equations.
        let m_obs = 6;
        let p = 2; // [x1, x2], has_cons = false
        let eps = 1e-7_f64;
        // Build orthonormal rows/vectors on R^6 via simple basis vectors e_0, e_1:
        // x1 = e_0, x2 = e_0 + eps * e_1 + 0.5 * e_2; z1 = e_0, z2 = e_1 -> P_Z X_tilde = [e_0, e_0 + eps*e_1]
        let mut x_raw = vec![0.0; m_obs * p];
        let mut z_excl = vec![0.0; m_obs * p];
        z_excl[0] = 1.0; // obs 0, z1
        z_excl[p + 1] = 1.0; // obs 1, z2
        x_raw[0] = 1.0;
        x_raw[1] = 1.0;
        x_raw[p + 1] = eps;
        x_raw[2 * p + 1] = 0.5;
        let b = [1.0, 2.0];
        let bread = [1.0, 0.0, 0.0, 1.0];
        let scores = vec![0.1; m_obs * p];

        let inputs = fixture_inputs(
            m_obs, m_obs, false, 2, p, &b, &bread, &scores, &x_raw, &z_excl, None,
        )
        .expect("ill-conditioned kappa ~ 2e7 should succeed without normal equations");
        let u = inputs.proj_u.as_ref().unwrap();
        let v = inputs.proj_v.as_ref().unwrap();
        for r in 0..p {
            for c in 0..p {
                let dot: f64 = (0..m_obs).map(|i| v[i * p + r] * u[i * p + c]).sum();
                let exp = if r == c { 1.0 } else { 0.0 };
                assert!((dot - exp).abs() < 1e-7, "V'U[{r},{c}]={dot}");
            }
        }

        // Case 2: Orthogonal / irrelevant instrument (e_0 vs e_1, e_2 -> zero projected column norm)
        let mut z_irrel = vec![0.0; m_obs * p];
        z_irrel[2 * p] = 1.0;
        z_irrel[3 * p + 1] = 1.0;
        let x_irrel = {
            let mut m = vec![0.0; m_obs * p];
            m[0] = 1.0;
            m[p + 1] = 1.0;
            m
        };
        let err = fixture_inputs(
            m_obs, m_obs, false, 2, p, &b, &bread, &scores, &x_irrel, &z_irrel, None,
        )
        .err()
        .expect("irrelevant instrument must fail rank check");
        assert!(err.contains("rank deficient"), "unexpected err: {err}");

        // Case 3: Pure numerical roundoff overlap (~1e-13) must also fail 1e-10 SVD rank check
        let mut z_roundoff = z_irrel.clone();
        z_roundoff[0] = 1e-13;
        let err_roundoff = fixture_inputs(
            m_obs,
            m_obs,
            false,
            2,
            p,
            &b,
            &bread,
            &scores,
            &x_irrel,
            &z_roundoff,
            None,
        )
        .err()
        .expect("1e-13 roundoff overlap must fail 1e-10 relative SVD check");
        assert!(
            err_roundoff.contains("rank deficient"),
            "unexpected err: {err_roundoff}"
        );
    }

    #[test]
    fn iv_svd_rejects_a_uniformly_negligible_projection() {
        let err = fixture_inputs(
            3,
            3,
            false,
            1,
            1,
            &[0.0],
            &[1.0],
            &[0.1; 3],
            &[1.0, 0.0, 0.0],
            &[1e-13, 1.0, 0.0],
            None,
        )
        .err()
        .expect("one nearly orthogonal instrument must fail rank validation");
        assert!(err.contains("rank deficient"), "unexpected err: {err}");
    }

    #[test]
    fn approx_scpc_args_normalization_permutation_and_cluster_remapping() {
        assert!(parse_seed("1").unwrap() == 1);
        assert!(parse_seed("4294967295").unwrap() == u32::MAX);
        assert!(parse_seed("-1").is_err());
        assert!(parse_seed("4294967296").is_err());

        let base20 = [
            "2", "3", "3", "1", "0.05", "0", "0", "1", "-1", "b", "bread", "stats", "cvs", "c0",
            "cv", "q", "tied", "1", "0", "0",
        ];
        assert!(handle_model_scpc(&base20[..19], true).is_err());
        let err = handle_model_scpc(&base20, true).unwrap_err();
        assert!(err.contains("expects 21 arguments"), "err={err}");

        let mut args21 = base20.to_vec();
        args21.push("7");
        let err21 = handle_model_scpc(&args21, true).unwrap_err();
        assert!(
            err21.contains("variables in plugin varlist"),
            "err21={err21}"
        );

        // Run the complete raw-model operation for a permuted, clustered design.
        let n_sites = 12usize;
        let per_cl = 3usize;
        let m_obs = n_sites * per_cl;
        let p = 2usize;
        let mut coords_cl = Vec::with_capacity(n_sites * 2);
        for g in 0..n_sites {
            // Reversed order so normalize_coords produces a non-identity permutation.
            let gf = (n_sites - 1 - g) as f64;
            coords_cl.push(gf * 0.25);
            coords_cl.push((gf * 1.7).sin() * 0.3);
        }
        let mut cid = Vec::with_capacity(m_obs);
        let mut x_raw = Vec::with_capacity(m_obs);
        let mut scores_obs = Vec::with_capacity(m_obs * p);
        for g in 0..n_sites {
            for r in 0..per_cl {
                let obs_i = g * per_cl + r;
                let xv = (obs_i as f64 * 0.41).cos() + 0.2 * (g as f64);
                let u = ((obs_i + 1) as f64 * 0.73).sin() * 0.5;
                cid.push(g);
                x_raw.push(xv);
                scores_obs.extend([u * xv, u]);
            }
        }
        let b = [1.2, 0.4];
        let bread = [1.0 / m_obs as f64, 0.0, 0.0, 1.0 / m_obs as f64];
        let design = StataDesign::from_columns(m_obs, p, true, 0, 0, &x_raw, &[]).unwrap();
        let loaded = LoadedStataScpc {
            n_sites,
            n_obs: m_obs,
            coords: coords_cl,
            metric: Metric::GreatCircleFractionLatLon,
            b: b.to_vec(),
            bread: bread.to_vec(),
            cluster_id: Some(cid),
            scores: scores_obs,
            design: Some(design),
        };
        let run = |threads| {
            run_model_scpc(
                loaded.prepare_request(p),
                &loaded.coords,
                2,
                loaded.metric,
                0.1,
                true,
                Some(threads),
                Some(5),
            )
            .unwrap()
        };
        let out_1 = run(1);
        let out_2 = run(2);
        assert_eq!(out_1.q, out_2.q);
        for (&a, &b) in out_1.scpcstats.iter().zip(out_2.scpcstats.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }
}
