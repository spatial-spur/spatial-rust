.spatial_runtime_env <- new.env(parent = emptyenv())
.spatial_runtime_env$init_pid <- NULL

.onLoad <- function(libname, pkgname) {
  .spatial_runtime_env$init_pid <- Sys.getpid()
}

.is_forked_child <- function() {
  init <- .spatial_runtime_env$init_pid
  if (is.null(init)) {
    .spatial_runtime_env$init_pid <- Sys.getpid()
    return(FALSE)
  }
  !identical(Sys.getpid(), init)
}

.is_r_check_limit_cores <- function() {
  val <- trimws(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""))
  nzchar(val) && !tolower(val) %in% c("false", "0", "no")
}

spatial_build_info <- function() .Call(wrap__rust_build_info)

.handle_native_result <- function(result) {
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
  threads <- if (!is.null(nthreads)) {
    as.integer(.integer_option(nthreads, "nthreads"))
  } else if (.is_forked_child()) {
    1L
  } else {
    NULL
  }
  if (.is_r_check_limit_cores() && (is.null(threads) || threads > 2L)) {
    threads <- 2L
  }
  c(spatial_build_info(), list(nthreads_requested = threads, backend = "rust"))
}

scpc <- function(model, data = NULL, lon = NULL, lat = NULL, coords_euclidean = NULL,
                 cluster = NULL, ncoef = NULL, avc = 0.03, method = "auto",
                 large_n_seed = 1, uncond = FALSE, cvs = FALSE, backend = NULL,
                 nthreads = NULL, coords = NULL, geographic = FALSE, row_id = NULL,
                 units = NULL, inference = NULL, landmarks = 1000L,
                 landmark_samples = 20L, cv_tolerance = 0.001) {
  execution <- .execution(backend, nthreads)
  uncond <- .scalar_logical(uncond, "uncond")
  cvs <- .scalar_logical(cvs, "cvs")
  geographic <- .scalar_logical(geographic, "geographic")
  avc <- .scalar_double(avc, "avc", gt = 0.001, lt = 0.99)
  cv_tolerance <- .scalar_double(cv_tolerance, "cv_tolerance", gt = 0)
  landmarks_val <- .integer_option(landmarks, "landmarks", 2)
  samples_val <- .integer_option(landmark_samples, "landmark_samples")
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
    if (!is.null(coords_euclidean) && geographic) {
      stop("Do not set geographic = TRUE with Euclidean coordinates", call. = FALSE)
    }
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
  sites <- .make_sites(coords, rows, cluster)
  if (identical(method, "approx")) method <- "nystrom"
  method <- match.arg(method, c("auto", "exact", "nystrom"))
  if (method == "auto") method <- if (nrow(sites$coords) >= 4500L) "nystrom" else "exact"
  approximate <- method == "nystrom"
  seed <- if (approximate) as.double(.integer_option(large_n_seed, "large_n_seed", 0, 2^32 - 1))
  if (approximate) {
    landmarks_val <- min(landmarks_val, nrow(sites$coords))
  }
  use <- seq_len(if (is.null(ncoef)) length(fit$estimates) else {
    .integer_option(ncoef, "ncoef", 1, length(fit$estimates))
  })
  coef_names <- names(fit$estimates)[use]
  conditional <- inference != "unconditional"
  out <- .handle_native_result(.Call(wrap__rust_run_model_scpc,
    sites$coords, geographic, avc, cvs, as.double(fit$estimates[use]),
    .real_matrix(fit$scores, "Fitted scores"), fit$bread[use, , drop = FALSE],
    if (conditional) fit$x, if (conditional) fit$z,
    if (!is.null(sites$group)) sites$group - 1L,
    if (conditional && !is.null(fit$fe)) fit$fe - 1L,
    execution$nthreads_requested, inference, method,
    seed, as.integer(landmarks_val), as.integer(samples_val), cv_tolerance,
    if (conditional) fit$residualize))
  dimnames(out$scpcstats) <- list(coef_names, c("coef", "se", "t", "p", "ci_lower", "ci_upper"))
  if (!is.null(out$scpccvs)) dimnames(out$scpccvs) <- list(coef_names, c("32%", "10%", "5%", "1%"))
  names(out$critical_95) <- coef_names
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
  out$approximation <- if (approximate) {
    list(seed = seed, landmarks = landmarks_val, samples = samples_val)
  }
  out$cv_tolerance <- cv_tolerance
  out$diagnostics <- list(tied_eigenvalues = out$tied)
  out$execution <- execution
  out$backend <- "rust"
  structure(out, class = "scpc")
}

print.scpc <- function(x, ...) {
  cat("SCPC ", x$inference, ", ", x$method, ", q = ", x$q, ", avc = ", format(x$avc),
      ", backend = rust (", x$execution$linalg_provider, ")\n", sep = "")
  if (identical(x$method, "nystrom") && !is.null(x$approximation)) {
    cat("Nystrom approximation: seed = ", x$approximation$seed,
        ", landmarks = ", x$approximation$landmarks,
        ", samples = ", x$approximation$samples, "\n", sep = "")
  }
  print(cbind(x$scpcstats, critical_95 = x$critical_95), ...)
  if (!is.null(x$scpccvs)) {
    cat("\nTwo-sided critical values:\n")
    print(x$scpccvs, ...)
  }
  if (isTRUE(x$diagnostics$tied_eigenvalues)) {
    cat("\nNote: Leading spatial eigenvalues exhibited a numerical tie.\n")
  }
  invisible(x)
}

coef.scpc <- function(object, ...) {
  object$scpcstats[, "coef"]
}

confint.scpc <- function(object, parm, level = 0.95, ...) {
  level <- .scalar_double(level, "level", gt = 0, lt = 1)
  cf <- object$scpcstats[, "coef"]
  pnames <- names(cf)
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }
  idx <- match(parm, pnames)
  if (anyNA(idx)) {
    stop("Unknown coefficient in parm", call. = FALSE)
  }
  if (abs(level - 0.95) < 1e-12) {
    ci <- object$scpcstats[idx, c("ci_lower", "ci_upper"), drop = FALSE]
  } else {
    cv_map <- c("0.68" = 1L, "0.9" = 2L, "0.95" = 3L, "0.99" = 4L)
    key <- as.character(round(level, 6))
    if (is.null(object$scpccvs) || !key %in% names(cv_map)) {
      stop("Confidence level ", level, " is not available (pass cvs = TRUE for 68%/90%/99%)",
           call. = FALSE)
    }
    cv <- object$scpccvs[idx, cv_map[[key]]]
    se <- object$scpcstats[idx, "se"]
    est <- cf[idx]
    ci <- cbind(est - cv * se, est + cv * se)
  }
  probs <- c((1 - level) / 2, (1 + level) / 2)
  colnames(ci) <- paste0(formatC(100 * probs, format = "fg", digits = 3), " %")
  rownames(ci) <- pnames[idx]
  ci
}

summary.scpc <- function(object, ...) {
  structure(list(object = object), class = "summary.scpc")
}

print.summary.scpc <- function(x, ...) {
  obj <- x$object
  cat("SCPC Inference Summary\n")
  cat("Inference:", obj$inference, "| Estimator:", obj$estimator, "| Method:", obj$method, "\n")
  cat("Observations:", obj$n_obs, "| Sites:", obj$n_sites, "| q:", obj$q, "| avc:", obj$avc, "\n")
  cat("Backend: rust (", obj$execution$linalg_provider, ")\n\n", sep = "")
  print(cbind(obj$scpcstats, critical_95 = obj$critical_95), ...)
  if (!is.null(obj$scpccvs)) {
    cat("\nTwo-sided critical values:\n")
    print(obj$scpccvs, ...)
  }
  if (isTRUE(obj$diagnostics$tied_eigenvalues)) {
    cat("\nNote: Leading spatial eigenvalues exhibited a numerical tie.\n")
  }
  invisible(x)
}
