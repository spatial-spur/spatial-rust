\name{spatialSpur}
\alias{spatialSpur}
\alias{scpc}
\alias{spur}
\alias{spurtest}
\alias{spurtest_i0}
\alias{spurtest_i1}
\alias{spurtest_i0resid}
\alias{spurtest_i1resid}
\alias{spurhalflife}
\alias{spurtransform}
\alias{load_chetty_data}
\alias{spatial_build_info}
\title{Spatial Correlation-Robust Inference and Spatial Unit Root Diagnostics}
\description{
Calculate spatial correlation-robust intervals for fitted unweighted regressions,
run spatial unit root diagnostics, invert spatial half-life tests, and transform
spatial observations. The functions use a shared native numerical engine.
}
\usage{
scpc(model, data = NULL, lon = NULL, lat = NULL, coords_euclidean = NULL,
  cluster = NULL, ncoef = NULL, avc = 0.03, method = "auto",
  large_n_seed = 1, uncond = FALSE, cvs = FALSE, backend = NULL,
  nthreads = NULL, coords = NULL, geographic = FALSE, row_id = NULL,
  units = NULL, inference = NULL, landmarks = 1000L,
  landmark_samples = 20L, cv_tolerance = 0.001)

spur(formula, data, lon = NULL, lat = NULL, coords_euclidean = NULL,
  q = 15, nrep = NULL, seed = 42L, avc = 0.03, uncond = FALSE,
  cvs = FALSE, backend = NULL, nthreads = NULL, draws = NULL,
  tuning_draws = NULL, method = "auto", large_n_seed = 1)

spurtest(formula, data, test, lon = NULL, lat = NULL, coords_euclidean = NULL,
  q = 15, nrep = NULL, seed = 42L, backend = NULL, nthreads = NULL,
  draws = NULL, tuning_draws = NULL)
spurtest_i0(formula, data, ...)
spurtest_i1(formula, data, ...)
spurtest_i0resid(formula, data, ...)
spurtest_i1resid(formula, data, ...)

spurhalflife(formula, data, lon = NULL, lat = NULL, coords_euclidean = NULL,
  q = 15, nrep = NULL, level = 95, normdist = FALSE, seed = 42L,
  backend = NULL, nthreads = NULL, draws = NULL)

spurtransform(formula, data, lon = NULL, lat = NULL, coords_euclidean = NULL,
  prefix = "h_", transformation = "lbmgls", radius = NULL,
  clustvar = NULL, backend = NULL, nthreads = NULL, replace = FALSE)

load_chetty_data()
spatial_build_info()
}
\arguments{
\item{model}{An unweighted \code{lm} or a \code{fixest::feols} OLS or 2SLS fit.
For fixest, retain fitted matrices with \code{demeaned = TRUE}. Save source rows
with \code{data.save = TRUE} to look up coordinates by name.}
\item{data}{A data frame containing coordinates and, for formula interfaces,
analysis variables.}
\item{lon, lat}{The names of geographic coordinate columns in degrees.}
\item{coords_euclidean}{Names of Euclidean coordinate columns.}
\item{coords}{A numeric coordinate matrix in fitted observation order.}
\item{geographic}{Whether \code{coords} contains longitude and latitude.}
\item{row_id}{A shared unique fitted row identifier column.}
\item{cluster}{A site label column or a fitted-order site label vector.}
\item{units}{A descriptive label for the supplied coordinate units.}
\item{ncoef}{The number of coefficients to report.}
\item{avc}{The upper bound on average pairwise correlation.}
\item{method}{Spatial geometry choice \code{"auto"}, \code{"exact"}, or \code{"nystrom"}.}
\item{large_n_seed, landmarks, landmark_samples}{Settings for the Nyström ensemble.}
\item{uncond, inference}{Select conditional OLS, conditional ordinary 2SLS, or
unconditional inference.}
\item{cvs}{Whether to return additional critical values.}
\item{cv_tolerance}{The critical-value tolerance in t units.}
\item{backend}{Native backend selector \code{"rust"} or \code{"auto"}.}
\item{nthreads}{Optional positive native thread budget.}
\item{formula}{A regression formula, one-response formula, or supplied variable
names, depending on the procedure.}
\item{test}{One of \code{"i0"}, \code{"i1"}, \code{"i0resid"}, or \code{"i1resid"}.}
\item{q}{The number of low-frequency spatial directions.}
\item{nrep}{Requested simulation draws. Defaults to 100,000 or the width of
the supplied draw matrix.}
\item{seed}{Seed for local draws. The caller's R RNG state is restored. Supply
draw matrices explicitly if the caller uses R's Box-Muller normal generator.}
\item{draws, tuning_draws}{Optional matrices with \code{q} rows of independent
standard Gaussian drivers for inference and tuning. Half-life uses one draw matrix.}
\item{level}{Confidence level in percent.}
\item{normdist}{Return half-life in sample diameters when true.}
\item{prefix}{Prefix for transformed columns.}
\item{transformation}{One of \code{"lbmgls"}, \code{"nn"}, \code{"iso"}, or \code{"cluster"}.}
\item{radius}{A strictly positive ISO radius in input units or meters for geography.}
\item{clustvar}{Site label column for cluster demeaning.}
\item{replace}{Allow replacement of existing output columns in the returned data.}
\item{...}{Additional arguments to \code{spurtest}.}
}
\value{
\code{scpc} returns coefficient estimates, uncertainty, effective critical values
and a spatial basis in the reported site order. Diagnostics return likelihood
ratios, simulation cutoffs, and p-values. Half-life returns a matrix of confidence
components with units and normalization. The pipeline returns diagnostics and
both regressions. The transform and data loader return data frames.
\code{spatial_build_info} returns the compiled methods version and provider.
}
\examples{
d <- load_chetty_data()
d <- d[complete.cases(d[c("am", "fracblack", "lat", "lon")]), ][1:60, ]
scpc(lm(am ~ fracblack, d), d, lon = "lon", lat = "lat", nthreads = 1L)
}
