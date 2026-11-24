# Fitted regression state and the distinct, explicitly ordered spatial sites.

.real_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) == 0L || ncol(x) == 0L || any(!is.finite(x))) {
    stop(name, " must be a nonempty finite numeric matrix", call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE", call. = FALSE)
  }
  x
}

.scalar_double <- function(x, name, gt = NULL, lt = NULL) {
  if (!is.numeric(x) || is.logical(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(name, " must be a single finite number", call. = FALSE)
  }
  val <- as.double(x)
  if (!is.null(gt) && !(val > gt)) stop(name, " must be greater than ", gt, call. = FALSE)
  if (!is.null(lt) && !(val < lt)) stop(name, " must be less than ", lt, call. = FALSE)
  val
}

.require_data_frame <- function(data) {
  if (!is.data.frame(data) || anyDuplicated(names(data))) {
    stop("data must be a data frame with unique column names", call. = FALSE)
  }
  data
}

.integer_option <- function(x, name, min = 1, max = .Machine$integer.max) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != trunc(x) || x < min || x > max) {
    stop(name, " must be an integer in the supported range", call. = FALSE)
  }
  x
}

.scaled_svd <- function(x, right = FALSE) {
  scale <- apply(abs(x), 2L, max)
  scale[scale == 0] <- 1
  decomp <- svd(sweep(x, 2L, scale, "/"), nv = if (right) ncol(x) else 0L)
  decomp$rank <- sum(decomp$d > 1e-10 * decomp$d[1L])
  decomp$scale <- scale
  decomp
}

.coefficient_bread <- function(x, coefficients) {
  decomp <- .scaled_svd(x, right = TRUE)
  if (decomp$rank < ncol(x)) {
    omitted <- seq.int(decomp$rank + 1L, ncol(x))
    aliased <- apply(abs(decomp$v[, omitted, drop = FALSE]), 1L, max) > 1e-7
    stop("Fitted coefficients are not separately identified: ",
         paste(coefficients[aliased], collapse = ", "), call. = FALSE)
  }
  factor <- sweep(sweep(decomp$v, 2L, decomp$d, "/"), 1L, decomp$scale, "/")
  tcrossprod(factor)
}

.has_offset <- function(formula, data = NULL) {
  terms_obj <- if (is.null(data)) {
    tryCatch(stats::terms(formula), error = function(e) NULL)
  } else {
    stats::terms(formula, data = data)
  }
  !is.null(terms_obj) && length(attr(terms_obj, "offset")) > 0L
}

.extract_lm <- function(model) {
  if (inherits(model, "mlm")) {
    stop("SCPC supports single-response unweighted lm fits", call. = FALSE)
  }
  if (!is.null(model$weights)) stop("SCPC requires an unweighted fit", call. = FALSE)
  if (!is.null(model$offset) ||
      (!is.null(model$terms) && .has_offset(stats::formula(model), model$model))) {
    stop("SCPC does not support regression offsets", call. = FALSE)
  }
  if (is.null(model$model) && is.null(model[["x"]])) {
    stop("Retain the fitted model frame or design with lm(model = TRUE) or lm(x = TRUE)",
         call. = FALSE)
  }
  coefficients <- stats::coef(model)
  x <- .real_matrix(if (is.null(model[["x"]])) stats::model.matrix(model) else model[["x"]],
                    "Fitted lm design")
  bread <- .coefficient_bread(x, names(coefficients))
  if (anyNA(coefficients)) {
    aliased <- names(coefficients)[is.na(coefficients)]
    stop("Fitted coefficients are not separately identified: ",
         paste(aliased, collapse = ", "), call. = FALSE)
  }
  # Use stored residuals. residuals() can reinstate omitted rows with na.exclude.
  list(x = x, z = NULL, estimates = coefficients, scores = x * model$residuals,
       bread = bread, rows = rownames(x), source = model$model, host = "lm")
}

.align_design <- function(x, columns, name) {
  if (is.null(x) || !all(columns %in% colnames(x))) {
    stop("Cannot align ", name, " to the fitted coefficients", call. = FALSE)
  }
  .real_matrix(x[, columns, drop = FALSE], name)
}

.extract_fixest <- function(model) {
  if (!identical(model$method, "feols") || !is.null(model$weights)) {
    stop("SCPC supports unweighted fixest::feols() and ordinary 2SLS", call. = FALSE)
  }
  if (!is.null(model$offset) || .has_offset(model$fml_all$linear, model$data)) {
    stop("SCPC does not support regression offsets", call. = FALSE)
  }
  if (!is.null(model$slope_flag) && any(model$slope_flag != 0L)) {
    stop("SCPC does not support varying-slope absorbed effects", call. = FALSE)
  }
  if (is.null(model$X_demeaned)) {
    stop("Retain the fitted design with fixest::feols(demeaned = TRUE)", call. = FALSE)
  }
  estimates <- stats::coef(model)
  if (anyNA(estimates)) {
    aliased <- names(estimates)[is.na(estimates)]
    stop("Fitted coefficients are not separately identified: ",
         paste(aliased, collapse = ", "), call. = FALSE)
  }
  direction <- .align_design(model$X_demeaned, names(estimates), "Stored fixest design")
  x <- direction
  z <- NULL
  if (isTRUE(model[["is_iv"]]) || isTRUE(model[["iv"]])) {
    if (is.null(model$iv_first_stage) || !identical(as.integer(model$iv_stage), 2L)) {
      stop("Retain the fitted fixest 2SLS model and its first stages", call. = FALSE)
    }
    endog <- model$iv_endo_names
    fitted_names <- model$iv_endo_names_fit
    for (i in seq_along(endog)) {
      first <- model$iv_first_stage[[endog[i]]]
      if (is.null(first$X_demeaned) || is.null(first$y_demeaned) ||
          length(first$y_demeaned) != nrow(x) || !fitted_names[i] %in% colnames(x)) {
        stop("First-stage matrices do not identify the fitted 2SLS sample", call. = FALSE)
      }
      if (is.null(z)) z <- .real_matrix(first$X_demeaned, "Stored fixest instruments")
      if (!isTRUE(all.equal(z, first$X_demeaned))) {
        stop("The fitted first stages use inconsistent instrument designs", call. = FALSE)
      }
      x[, fitted_names[i]] <- first$y_demeaned
    }
    z_svd <- .scaled_svd(z)
    qz <- z_svd$u[, seq_len(z_svd$rank), drop = FALSE]
    fitted <- qz %*% crossprod(qz, x)
    if (!isTRUE(all.equal(unname(fitted), unname(direction), tolerance = 1e-7))) {
      stop("Stored first stages do not reproduce the fitted 2SLS projections", call. = FALSE)
    }
  }
  groups <- model$fixef_id
  if (length(groups) > 1L) {
    residualize <- function(x) {
      .real_matrix(fixest::demean(x, f = groups, tol = 1e-10,
                                 iter = max(10000L, model$fixef.iter),
                                 nthreads = 1L), "Within spatial directions")
    }
  } else {
    residualize <- NULL
  }
  source <- model$data
  if (!is.null(source) && !is.null(model$obs_selection)) {
    sel <- model$obs_selection
    nms <- names(sel)
    for (i in seq_along(sel)) {
      if (!is.null(nms) && identical(nms[i], "subset")) next
      source <- source[sel[[i]], , drop = FALSE]
    }
  }
  if (!is.null(source) && nrow(source) != nrow(x)) {
    stop("Stored fixest estimation rows do not match the fitted design", call. = FALSE)
  }
  list(x = x, z = z, estimates = estimates, scores = direction * model$residuals,
       bread = .coefficient_bread(direction, names(estimates)), source = source,
       rows = if (is.null(source)) seq_len(nrow(x)) else rownames(source), host = "fixest",
       fe = if (length(groups) == 1L) match(groups[[1L]], unique(groups[[1L]])) else NULL,
       residualize = residualize)
}

.fitted_model <- function(model) {
  if (inherits(model, "lm") && !inherits(model, "glm")) return(.extract_lm(model))
  if (inherits(model, "fixest")) return(.extract_fixest(model))
  stop("SCPC supports unweighted lm, fixest OLS, and fixest ordinary 2SLS", call. = FALSE)
}

.coordinate_columns <- function(lon, lat, coords_euclidean) {
  if (!is.null(coords_euclidean)) {
    if (!is.null(lon) || !is.null(lat) || !is.character(coords_euclidean) ||
        length(coords_euclidean) == 0L || anyNA(coords_euclidean) ||
        any(!nzchar(coords_euclidean)) || anyDuplicated(coords_euclidean)) {
      stop("Supply distinct valid coordinate column names for coords_euclidean", call. = FALSE)
    }
    return(coords_euclidean)
  }
  if (!is.character(lon) || !is.character(lat) || length(lon) != 1L ||
      length(lat) != 1L || is.na(lon) || is.na(lat) || !nzchar(lon) || !nzchar(lat) ||
      identical(lon, lat)) {
    stop("Supply both lon and lat, or coords_euclidean", call. = FALSE)
  }
  c(lon, lat)
}

.canonical_coordinates <- function(coords, geographic) {
  coords <- .real_matrix(coords, "Coordinates")
  if (geographic) {
    if (ncol(coords) != 2L || any(abs(coords[, 2L]) > 90)) {
      stop("Geographic coordinates require longitude and latitude within [-90, 90]",
           call. = FALSE)
    }
    longitude <- coords[, 1L]
    wrap <- longitude < -180 | longitude >= 180
    longitude[wrap] <- (longitude[wrap] + 180) %% 360 - 180
    longitude[abs(coords[, 2L]) == 90] <- 0
    coords[, 1L] <- longitude
  }
  coords
}

.read_coordinates <- function(data, columns, geographic) {
  values <- lapply(columns, function(name) {
    column <- data[[name]]
    if (!is.numeric(column) || !is.null(dim(column)) || any(!is.finite(column))) {
      stop("Coordinate column '", name, "' must be numeric and finite", call. = FALSE)
    }
    as.double(column)
  })
  .canonical_coordinates(do.call(cbind, values), geographic)
}

.match_frame <- function(fit, data, row_id, spatial) {
  source <- fit$source
  if (is.null(source) && fit$host == "fixest") {
    stop("Save fitted row identity with fixest::feols(data.save = TRUE), or supply coords in fitted order",
         call. = FALSE)
  }
  if (is.null(data)) {
    if (is.null(source)) stop("Supply data or coords in fitted order", call. = FALSE)
    data <- source
  }
  .require_data_frame(data)
  if (!is.null(row_id)) {
    if (!is.character(row_id) || length(row_id) != 1L || is.na(row_id) || !nzchar(row_id)) {
      stop("row_id must be a single column name", call. = FALSE)
    }
    if (!row_id %in% names(data)) {
      stop("Supplied data does not contain row_id '", row_id, "'", call. = FALSE)
    }
    if (!is.null(source) && row_id %in% names(source)) {
      rows <- source[[row_id]]
      labels <- data[[row_id]]
    } else if (fit$host == "lm" && !anyNA(match(fit$rows, rownames(data))) &&
               !anyNA(data[[row_id]]) && !anyDuplicated(data[[row_id]])) {
      rows <- fit$rows
      labels <- rownames(data)
    } else {
      stop("Fitted model data does not contain row_id '", row_id,
           "'; retain it in the fitted data or match by row names", call. = FALSE)
    }
  } else {
    rows <- fit$rows
    labels <- rownames(data)
  }
  if (length(rows) != nrow(fit$x) || anyNA(rows) || anyNA(labels) ||
      anyDuplicated(rows) || anyDuplicated(labels)) {
    stop("The row identifiers must uniquely identify observations in both frames", call. = FALSE)
  }
  position <- match(rows, labels)
  if (anyNA(position)) stop("Supplied data does not identify every fitted row", call. = FALSE)
  out <- data[position, , drop = FALSE]
  common <- setdiff(intersect(names(source), names(out)), spatial)
  if (length(common) && !all(vapply(common, function(name) {
    isTRUE(all.equal(source[[name]], out[[name]], check.attributes = FALSE, tolerance = 0))
  }, logical(1)))) {
    stop("Supplied data values disagree with the fitted observations", call. = FALSE)
  }
  list(data = out, rows = if (!is.null(row_id)) out[[row_id]] else rows)
}

.make_sites <- function(coords, rows, group) {
  if (nrow(coords) != length(rows)) {
    stop("Coordinates must have one row per fitted observation", call. = FALSE)
  }
  if (is.null(group)) return(list(coords = coords, ids = rows, group = NULL))
  if (length(group) != length(rows) || !is.null(dim(group)) || anyNA(group)) {
    stop("Supply one nonmissing cluster label per fitted observation", call. = FALSE)
  }
  ids <- unique(group)
  group_id <- match(group, ids)
  sites <- coords[match(ids, group), , drop = FALSE]
  if (any(coords != sites[group_id, , drop = FALSE])) {
    stop("Coordinates must be constant within each spatial cluster", call. = FALSE)
  }
  list(coords = sites, ids = ids, group = group_id)
}
