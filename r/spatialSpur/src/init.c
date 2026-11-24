#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP wrap__rust_build_info(void);

extern SEXP wrap__rust_run_spur_i0(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP emat,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_i1(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP emat,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_i0_resid(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP x_in,
    SEXP emat,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_i1_resid(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP x_in,
    SEXP emat,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_halflife(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP emat,
    SEXP level,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_transform(
    SEXP kind,
    SEXP n,
    SEXP y,
    SEXP coords,
    SEXP latlong,
    SEXP radius,
    SEXP cluster_id,
    SEXP nthreads
);

extern SEXP wrap__rust_run_spur_pipeline(
    SEXP coords,
    SEXP latlong,
    SEXP y,
    SEXP x_in,
    SEXP emat,
    SEXP transform_vars,
    SEXP nthreads
);

extern SEXP wrap__rust_run_model_scpc(
    SEXP coords, SEXP latlong, SEXP avc, SEXP cvs, SEXP coef, SEXP scores,
    SEXP bread, SEXP design, SEXP instruments, SEXP cluster_id, SEXP fe_id, SEXP nthreads,
    SEXP inference, SEXP method, SEXP seed, SEXP landmarks, SEXP samples,
    SEXP cv_tolerance, SEXP residualizer
);

static const R_CallMethodDef CallEntries[] = {
    {"wrap__rust_run_model_scpc", (DL_FUNC) &wrap__rust_run_model_scpc, 19},
    {"wrap__rust_build_info", (DL_FUNC) &wrap__rust_build_info, 0},
    {"wrap__rust_run_spur_i0", (DL_FUNC) &wrap__rust_run_spur_i0, 5},
    {"wrap__rust_run_spur_i1", (DL_FUNC) &wrap__rust_run_spur_i1, 5},
    {"wrap__rust_run_spur_i0_resid", (DL_FUNC) &wrap__rust_run_spur_i0_resid, 6},
    {"wrap__rust_run_spur_i1_resid", (DL_FUNC) &wrap__rust_run_spur_i1_resid, 6},
    {"wrap__rust_run_spur_halflife", (DL_FUNC) &wrap__rust_run_spur_halflife, 6},
    {"wrap__rust_run_spur_transform", (DL_FUNC) &wrap__rust_run_spur_transform, 8},
    {"wrap__rust_run_spur_pipeline", (DL_FUNC) &wrap__rust_run_spur_pipeline, 7},
    {NULL, NULL, 0}
};

SEXP sp_r_nil(void) { return R_NilValue; }
SEXP sp_r_names_sym(void) { return R_NamesSymbol; }
SEXP sp_r_class_sym(void) { return R_ClassSymbol; }
SEXP sp_r_global_env(void) { return R_GlobalEnv; }

/* Some macOS R builds leave attribute_visible empty even under Clang. */
#if defined(__GNUC__) && !defined(_WIN32) && !defined(HAVE_VISIBILITY_ATTRIBUTE)
__attribute__((visibility("default")))
#endif
void attribute_visible R_init_spatialSpur(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
