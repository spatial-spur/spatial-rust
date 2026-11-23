spatial_build_info <- function() .Call(wrap__rust_build_info)

.native_call <- function(symbol, ...) {
  result <- .Call(symbol, ...)
  if (inherits(result, "extendr_error")) {
    detail <- result$value
    if (inherits(detail, "condition")) stop(detail)
    if (is.list(detail) && !is.null(detail$message)) {
      stop(structure(list(message = detail$message, call = NULL),
                     class = c(detail$class, "error", "condition")))
    }
    stop(as.character(detail), call. = FALSE)
  }
  result
}

.execution <- function(backend, nthreads) {
  if (!is.null(backend) && !identical(backend, "rust") && !identical(backend, "auto")) {
    stop("backend must be 'rust' or 'auto'", call. = FALSE)
  }
  threads <- if (is.null(nthreads)) NULL else as.integer(.integer_option(nthreads, "nthreads"))
  c(spatial_build_info(), list(nthreads_requested = threads, backend = "rust"))
}

scpc <- function(model, data = NULL, lon = NULL, lat = NULL, coords_euclidean = NULL,
                 cluster = NULL, ncoef = NULL, avc = 0.03, method = "auto",
                 large_n_seed = 1, uncond = FALSE, cvs = FALSE, backend = NULL,
                 nthreads = NULL, coords = NULL, geographic = FALSE, row_id = NULL,
                 units = NULL, inference = NULL, landmarks = 1000L,
                 landmark_samples = 20L, cv_tolerance = 0.001) {
  execution <- .execution(backend, nthreads)
  fit <- .fitted_model(model)
  estimator <- if (is.null(fit$z)) "ols" else "2sls"
  conditional <- paste0("conditional_", estimator)
  if (is.null(inference)) inference <- if (uncond) "unconditional" else conditional
  if (!inference %in% c(conditional, "unconditional") ||
      (uncond && inference != "unconditional")) {
    stop("The requested inference is incompatible with this fitted ", estimator, " model",
         call. = FALSE)
  }
  rows <- fit$rows
  if (is.null(coords)) {
    columns <- .coordinate_columns(lon, lat, coords_euclidean)
    cluster_column <- is.character(cluster) && length(cluster) == 1L
    aligned <- .match_frame(fit, data, row_id, c(columns, if (cluster_column) cluster))
    rows <- aligned$rows
    geographic <- !is.null(lon)
    coords <- .read_coordinates(aligned$data, columns, geographic)
    if (cluster_column) {
      if (!cluster %in% names(aligned$data)) {
        stop("Spatial cluster column '", cluster, "' is missing", call. = FALSE)
      }
      cluster <- aligned$data[[cluster]]
    }
  } else {
    if (!all(vapply(list(data, lon, lat, coords_euclidean, row_id), is.null, logical(1)))) {
      stop("With explicit coords, supply coordinates and labels in fitted row order", call. = FALSE)
    }
    coords <- .canonical_coordinates(coords, geographic)
  }
  sites <- .make_sites(coords, rows, cluster, geographic)
  if (identical(method, "approx")) method <- "nystrom"
  method <- match.arg(method, c("auto", "exact", "nystrom"))
  if (method == "auto") method <- if (nrow(sites$coords) >= 4500L) "nystrom" else "exact"
  approximate <- method == "nystrom"
  seed <- if (approximate) as.double(.integer_option(large_n_seed, "large_n_seed", 0, 2^32 - 1))
  if (approximate) {
    landmarks <- min(.integer_option(landmarks, "landmarks", 2), nrow(sites$coords))
    landmark_samples <- .integer_option(landmark_samples, "landmark_samples")
  }
  use <- seq_len(if (is.null(ncoef)) length(fit$estimates) else {
    .integer_option(ncoef, "ncoef", 1, length(fit$estimates))
  })
  names <- names(fit$estimates)[use]
  conditional <- inference != "unconditional"
  out <- .native_call(wrap__rust_run_model_scpc,
    sites$coords, geographic, as.double(avc), cvs, as.double(fit$estimates[use]),
    .real_matrix(fit$scores, "Fitted scores"), fit$bread[use, , drop = FALSE],
    if (conditional) fit$x, if (conditional) fit$z,
    if (!is.null(sites$group)) sites$group - 1L,
    if (conditional && !is.null(fit$fe)) fit$fe - 1L,
    execution$nthreads_requested, inference, method,
    seed, as.integer(landmarks), as.integer(landmark_samples), as.double(cv_tolerance),
    if (conditional) fit$residualize)
  dimnames(out$scpcstats) <- list(names, c("coef", "se", "t", "p", "ci_lower", "ci_upper"))
  if (!is.null(out$scpccvs)) dimnames(out$scpccvs) <- list(names, c("32%", "10%", "5%", "1%"))
  names(out$critical_95) <- names
  rownames(out$wfin) <- as.character(sites$ids)
  out$W <- out$wfin
  out$wfin <- NULL
  out$cv <- out$cvfin
  out$cvfin <- NULL
  out$inference <- inference
  out$estimator <- estimator
  out$method <- method
  out$avc <- avc
  out$observations <- rows
  out$sites <- sites$ids
  out$coordinates <- sites$coords
  out$site_membership <- if (is.null(sites$group)) seq_along(rows) else sites$group
  out$n_obs <- length(rows)
  out$n_sites <- length(sites$ids)
  out$coordinate_units <- if (is.null(units)) {
    if (geographic) "fraction of Earth circumference" else "coordinate units"
  } else units
  out$approximation <- if (approximate) list(seed = seed, landmarks = landmarks, samples = landmark_samples)
  out$cv_tolerance <- cv_tolerance
  out$diagnostics <- list(tied_eigenvalues = out$tied)
  out$execution <- execution
  out$backend <- "rust"
  structure(out, class = "scpc")
}

print.scpc <- function(x, ...) {
  cat("SCPC", x$inference, x$method, "q =", x$q, "\n")
  print(cbind(x$scpcstats, critical_95 = x$critical_95), ...)
  invisible(x)
}
