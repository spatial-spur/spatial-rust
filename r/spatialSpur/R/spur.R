.local_draws <- function(rows, nrep, seed) {
  if (RNGkind()[2L] == "Box-Muller") {
    stop("Supply draws explicitly when using R's Box-Muller generator", call. = FALSE)
  }
  old <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  })
  if (!is.null(seed)) set.seed(.integer_option(seed, "seed", 0))
  matrix(stats::rnorm(rows * as.double(nrep)), nrow = rows, ncol = nrep)
}

.draw_blocks <- function(q, nrep, seed, draws, tuning_draws, tune) {
  q <- .integer_option(q, "q")
  if (is.null(draws)) {
    if (!is.null(tuning_draws)) stop("Supply an inference draw block", call. = FALSE)
    nrep <- if (is.null(nrep)) 100000L else .integer_option(nrep, "nrep")
    values <- .local_draws(if (tune) 2L * q else q, nrep, seed)
  } else {
    draws <- .real_matrix(draws, "Inference draws")
    if (nrow(draws) != q || (!is.null(nrep) &&
        .integer_option(nrep, "nrep") != ncol(draws))) {
      stop("Inference draws must have dimensions q by nrep", call. = FALSE)
    }
    seed <- NULL
    nrep <- ncol(draws)
    if (tune) {
      tuning_draws <- .real_matrix(tuning_draws, "Tuning draws")
      if (!identical(dim(draws), dim(tuning_draws))) {
        stop("Inference and tuning draws must have the same dimensions", call. = FALSE)
      }
      values <- rbind(draws, tuning_draws)
    } else {
      values <- draws
    }
  }
  list(values = values, simulation = list(q = q, requested_draws = nrep,
    inference_draws = nrep, tuning_draws = if (tune) nrep else 0L, seed = seed))
}

.spur_coords <- function(data, lon, lat, coords_euclidean) {
  .read_coordinates(data, .coordinate_columns(lon, lat, coords_euclidean), !is.null(lon))
}

.response <- function(data, name) {
  values <- data[[name]]
  if (!is.numeric(values) || !is.null(dim(values)) || any(!is.finite(values))) {
    stop("Response '", name, "' must name a finite numeric column", call. = FALSE)
  }
  as.double(values)
}

.levels_name <- function(formula) {
  if (is.character(formula) && length(formula) == 1L && !grepl("~", formula, fixed = TRUE)) {
    return(formula)
  }
  formula <- stats::as.formula(formula)
  if (length(formula) == 3L && !identical(as.character(formula[[3L]]), "1")) {
    stop("Levels diagnostics take one response column", call. = FALSE)
  }
  lhs <- formula[[2L]]
  if (!is.symbol(lhs)) stop("Specify a response column name", call. = FALSE)
  as.character(lhs)
}

.spur_fit <- function(formula, data) {
  formula <- stats::as.formula(formula)
  if (.has_offset(formula)) stop("SPUR does not support a regression offset", call. = FALSE)
  model <- stats::lm(formula, data = data, x = TRUE, y = TRUE, na.action = stats::na.omit)
  x <- .real_matrix(model$x, "Fitted design")
  if (!is.numeric(model$y) || !is.null(dim(model$y)) || any(!is.finite(model$y))) {
    stop("SPUR requires one finite numeric response", call. = FALSE)
  }
  .coefficient_bread(x, colnames(x))
  list(model = model, x = x, y = as.double(model$y),
       sample = data[match(rownames(x), rownames(data)), , drop = FALSE])
}

.diagnostic <- function(out, type, rows, simulation, longitude, execution) {
  if (is.null(out$cvalue)) out$cvalue <- out$cv_vec
  out$test <- type
  out$observations <- rows
  out$n_obs <- length(rows)
  out$n_sites <- length(rows)
  out$coordinate_units <- if (longitude) "meters" else "coordinate units"
  out$normalization <- "sample diameter"
  out$simulation <- simulation
  out$execution <- execution
  out$backend <- "rust"
  out
}

spurtest <- function(formula, data, test, lon = NULL, lat = NULL, coords_euclidean = NULL,
                     q = 15, nrep = NULL, seed = 42L, backend = NULL, nthreads = NULL,
                     draws = NULL, tuning_draws = NULL) {
  test <- match.arg(test, c("i0", "i1", "i0resid", "i1resid"))
  execution <- .execution(backend, nthreads)
  residual <- endsWith(test, "resid")
  if (residual) {
    fit <- .spur_fit(formula, data)
    sample <- fit$sample
    y <- fit$y
  } else {
    sample <- data
    y <- .response(sample, .levels_name(formula))
  }
  coords <- .spur_coords(sample, lon, lat, coords_euclidean)
  q <- .integer_option(q, "q", 2, nrow(sample) - if (residual) ncol(fit$x) else 1L)
  sim <- .draw_blocks(q, nrep, seed, draws, tuning_draws, TRUE)
  args <- list(coords, !is.null(lon), matrix(y, ncol = 1L))
  if (residual) args <- c(args, list(fit$x))
  symbol <- switch(test,
    i0 = wrap__rust_run_spur_i0, i1 = wrap__rust_run_spur_i1,
    i0resid = wrap__rust_run_spur_i0_resid, i1resid = wrap__rust_run_spur_i1_resid)
  out <- do.call(.native_call, c(list(symbol), args, list(sim$values, execution$nthreads_requested)))
  .diagnostic(out, test, rownames(sample), sim$simulation, !is.null(lon), execution)
}

spurtest_i0 <- function(formula, data, ...) spurtest(formula, data, "i0", ...)
spurtest_i1 <- function(formula, data, ...) spurtest(formula, data, "i1", ...)
spurtest_i0resid <- function(formula, data, ...) spurtest(formula, data, "i0resid", ...)
spurtest_i1resid <- function(formula, data, ...) spurtest(formula, data, "i1resid", ...)

spurhalflife <- function(formula, data, lon = NULL, lat = NULL, coords_euclidean = NULL,
                         q = 15, nrep = NULL, level = 95, normdist = FALSE, seed = 42L,
                         backend = NULL, nthreads = NULL, draws = NULL) {
  execution <- .execution(backend, nthreads)
  q <- .integer_option(q, "q", 2, min(100L, nrow(data) - 1L))
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 100) {
    stop("level must lie strictly between 0 and 100", call. = FALSE)
  }
  geographic <- !is.null(lon)
  coords <- .spur_coords(data, lon, lat, coords_euclidean)
  y <- .response(data, .levels_name(formula))
  sim <- .draw_blocks(q, nrep, seed, draws, NULL, FALSE)
  out <- .native_call(wrap__rust_run_spur_halflife,
    coords, geographic, y, sim$values, as.double(level / 100), execution$nthreads_requested)
  out$max_dist <- out$max_dist * if (geographic) execution$earth_circumference_m else 1
  divisor <- if (normdist) 1 else out$max_dist
  out$components <- out$components * divisor
  colnames(out$components) <- c("lower", "upper")
  out$hl_ci_lower <- out$hl_ci_lower * divisor
  out$hl_ci_upper <- out$hl_ci_upper * divisor
  out$hl_grid_ho <- out$hl_grid_ho * divisor
  out$level <- level
  out$normdist <- normdist
  out$distance_units <- if (geographic) "meters" else "coordinate units"
  out$normalization <- "sample diameter"
  out$simulation <- sim$simulation
  out$observations <- rownames(data)
  out$n_obs <- nrow(data)
  out$n_sites <- nrow(data)
  out$backend <- "rust"
  out$execution <- execution
  structure(out, class = "spatial_spur_halflife")
}

print.spatial_spur_halflife <- function(x, ...) {
  units <- if (x$normdist) "sample diameters" else x$distance_units
  cat(x$level, "% spatial half-life confidence set (", units, ")\n", sep = "")
  print(x$components, ...)
  invisible(x)
}

spurtransform <- function(formula, data, lon = NULL, lat = NULL, coords_euclidean = NULL,
                          prefix = "h_", transformation = "lbmgls", radius = NULL,
                          clustvar = NULL, backend = NULL, nthreads = NULL, replace = FALSE) {
  execution <- .execution(backend, nthreads)
  if (inherits(formula, "formula") && .has_offset(formula)) {
    stop("SPUR transformation does not support offsets", call. = FALSE)
  }
  columns <- if (inherits(formula, "formula")) all.vars(formula) else formula
  outputs <- paste0(prefix, columns)
  if (!is.character(columns) || length(columns) == 0L || anyDuplicated(columns) ||
      !is.character(prefix) || length(prefix) != 1L || !nzchar(prefix)) {
    stop("Supply distinct response columns and a nonempty output prefix", call. = FALSE)
  }
  if (!replace && any(outputs %in% names(data))) {
    stop("The requested output columns already exist", call. = FALSE)
  }
  transformation <- match.arg(transformation, c("lbmgls", "nn", "iso", "cluster"))
  if ((transformation == "iso") != !is.null(radius)) {
    stop("Supply radius only for the ISO transformation", call. = FALSE)
  }
  if ((transformation == "cluster") != !is.null(clustvar)) {
    stop("Supply clustvar only for the cluster transformation", call. = FALSE)
  }
  y <- do.call(cbind, lapply(columns, function(column) .response(data, column)))
  coords <- groups <- NULL
  geographic <- !is.null(lon) && is.null(clustvar)
  if (is.null(clustvar)) {
    coords <- .spur_coords(data, lon, lat, coords_euclidean)
  } else {
    labels <- data[[clustvar]]
    if (is.null(labels) || anyNA(labels)) {
      stop("clustvar must name a column of nonmissing site labels", call. = FALSE)
    }
    groups <- match(labels, unique(labels)) - 1L
  }
  if (!is.null(radius)) radius <- as.double(radius) / if (geographic) execution$earth_circumference_m else 1
  transformed <- .native_call(wrap__rust_run_spur_transform, transformation, nrow(data),
    y, coords, geographic, radius, groups, execution$nthreads_requested)$transformed
  data[outputs] <- as.data.frame(transformed)
  data
}

spur <- function(formula, data, lon = NULL, lat = NULL, coords_euclidean = NULL,
                 q = 15, nrep = NULL, seed = 42L, avc = 0.03, uncond = FALSE,
                 cvs = FALSE, backend = NULL, nthreads = NULL, draws = NULL,
                 tuning_draws = NULL, method = "auto", large_n_seed = 1) {
  execution <- .execution(backend, nthreads)
  fit <- .spur_fit(formula, data)
  intercept <- which(attr(fit$x, "assign") == 0L)
  if (length(intercept) != 1L) stop("The SPUR pipeline requires a fitted intercept", call. = FALSE)
  q <- .integer_option(q, "q", 2, nrow(fit$x) - ncol(fit$x))
  numeric_inputs <- cbind(fit$y, fit$x[, -intercept, drop = FALSE])
  coords <- .spur_coords(fit$sample, lon, lat, coords_euclidean)
  geographic <- !is.null(lon)
  sim <- .draw_blocks(q, nrep, seed, draws, tuning_draws, TRUE)
  out <- .native_call(wrap__rust_run_spur_pipeline, coords, geographic, fit$y, fit$x,
    sim$values, numeric_inputs, execution$nthreads_requested)
  diagnostics <- setNames(lapply(c("i0", "i1", "i0resid", "i1resid"), function(type) {
    .diagnostic(out[[type]], type, rownames(fit$sample), sim$simulation, geographic, execution)
  }), c("i0", "i1", "i0resid", "i1resid"))
  transformed <- fit$x
  transformed[, -intercept] <- out$transformed_vars[, -1L, drop = FALSE]
  transformed_fit <- stats::lm.fit(transformed, out$transformed_vars[, 1L])
  transformed_fit$x <- transformed
  transformed_fit$terms <- fit$model$terms
  class(transformed_fit) <- "lm"
  levels <- scpc(fit$model, coords = coords, geographic = geographic, avc = avc,
    uncond = uncond, cvs = cvs, backend = backend, nthreads = nthreads,
    method = method, large_n_seed = large_n_seed)
  adjusted <- scpc(transformed_fit, coords = coords, geographic = geographic, avc = avc,
    uncond = uncond, cvs = cvs, backend = backend, nthreads = nthreads,
    method = method, large_n_seed = large_n_seed)
  list(tests = diagnostics, levels = levels, transformed = adjusted,
       transformed_vars = out$transformed_vars, observations = rownames(fit$sample))
}

load_chetty_data <- function() {
  utils::read.csv(system.file("extdata", "chetty.csv", package = "spatialSpur"))
}
