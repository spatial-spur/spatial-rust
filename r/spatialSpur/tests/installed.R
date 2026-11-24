library(spatialSpur)

assert_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

set.seed(937)
n <- 84L
d <- data.frame(y = rnorm(n), x = rnorm(n), w = rnorm(n), z = rnorm(n),
                lon = -104 + 8 * runif(n), lat = 32 + 10 * runif(n),
                g = rep(1:7, 12), h = rep(1:6, each = 14),
                id = 200 + seq_len(n), row.names = paste0("id", seq_len(n)))
d$u <- 1.7 * d$z + rnorm(n)
d$y <- d$y + d$x + d$u
d$y[seq(3, n, by = 5)] <- NA
fit <- lm(y ~ x, d, na.action = na.exclude)
one <- scpc(fit, d, lon = "lon", lat = "lat", nthreads = 1L, backend = "rust", cvs = TRUE)
reorder <- sample(n)
two <- scpc(fit, d[reorder, ], lon = "lon", lat = "lat", nthreads = 1L, cvs = TRUE)
stopifnot(identical(one$scpcstats, two$scpcstats), one$n_obs == sum(!is.na(d$y)),
          identical(one$observations, rownames(d)[!is.na(d$y)]),
          max(abs(one$critical_95 - one$scpccvs[, "5%"])) < 1e-12,
          identical(one$execution$spec_version, spatial_build_info()$spec_version))

with_alias <- transform(d, alias = 2 * x)
for (uncond in c(FALSE, TRUE)) {
  assert_error(scpc(lm(y ~ x + alias, with_alias), with_alias, lon = "lon", lat = "lat",
                    uncond = uncond), "not separately identified.*x, alias")
}
assert_error(scpc(lm(y ~ x + offset(w), d), d, lon = "lon", lat = "lat"),
             "does not support regression offsets")
assert_error(spurtest_i1resid(y ~ x + offset(w), d, lon = "lon", lat = "lat"),
             "does not support a regression offset")

if (requireNamespace("fixest", quietly = TRUE)) {
  d$clipw <- pmax(quantile(d$w, .1), pmin(d$w, quantile(d$w, .9)))
  a <- fixest::feols(y ~ x + pmax(quantile(w, .1), pmin(w, quantile(w, .9))) | g + h,
                     d, demeaned = TRUE, data.save = TRUE, notes = FALSE)
  rows <- fixest::obs(a)
  # Compare with a literal precomputed column on the sample over which fixest evaluated it.
  b <- fixest::feols(y ~ x + clipw | g + h, d, demeaned = TRUE, data.save = TRUE, notes = FALSE)
  stopifnot(isTRUE(all.equal(unname(a$X_demeaned), unname(b$X_demeaned), tolerance = 1e-12)))
  opts <- list(lon = "lon", lat = "lat", nthreads = 2L, backend = "rust")
  result <- do.call(scpc, c(list(a, d[reorder, ], row_id = "id"), opts))
  literal <- do.call(scpc, c(list(b, d), opts))
  stopifnot(isTRUE(all.equal(unname(result$scpcstats), unname(literal$scpcstats), tolerance = 1e-11)),
            identical(result$observations, d$id[rows]))
  iv <- fixest::feols(y ~ x | g + h | u ~ z, d, demeaned = TRUE, data.save = TRUE, notes = FALSE)
  fit_iv <- do.call(scpc, c(list(iv, d), opts))
  stopifnot(fit_iv$estimator == "2sls", fit_iv$inference == "conditional_2sls")
  fi <- iv$iv_first_stage[["u"]]
  x <- iv$X_demeaned
  x[, "fit_u"] <- fi$y_demeaned
  direction <- qr.fitted(qr(fi$X_demeaned), x)
  bread <- solve(crossprod(direction))
  target <- sqrt(fit_iv$n_sites / fit_iv$q) *
    sqrt(colSums(crossprod(fit_iv$W[, -1], (direction * iv$residuals) %*% bread)^2))
  stopifnot(max(abs(target - fit_iv$scpcstats[, "se"])) < 1e-12)
}

opts <- list(lon = "lon", lat = "lat", q = 4L, nrep = 200L, nthreads = 1L, backend = "rust")
set.seed(384)
control <- replicate(4, rnorm(n))
set.seed(384)
observed <- replicate(4, {
  values <- rnorm(n)
  do.call(spurtest_i0, c(list("y", transform(d, y = values)), opts))
  values
})
stopifnot(identical(control, observed))
rm(".Random.seed", envir = .GlobalEnv)
invisible(do.call(spurtest_i0, c(list("x", d), opts)))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
set.seed(52)
caller <- .Random.seed
assert_error(do.call(spurtest_i0, c(list("x", transform(d, x = 7)), opts)), "uninformative|constant|variation")
stopifnot(identical(caller, .Random.seed))
assert_error(do.call(spurtest_i0, c(list("x", d, seed = -1), opts)), "seed")
stopifnot(identical(caller, .Random.seed))

distances <- do.call(spurhalflife, c(list("x", d), opts))
relative <- do.call(spurhalflife, c(list("x", d, normdist = TRUE), opts))
stopifnot(distances$distance_units == "meters", relative$normdist,
          identical(distances$components, relative$components * distances$max_dist),
          any(grepl("sample diameters", capture.output(print(relative)))))
p <- do.call(spur, c(list(y ~ x, d), opts))
stopifnot(p$levels$n_obs == p$transformed$n_obs, p$levels$backend == "rust",
          identical(rownames(p$transformed$scpcstats), "x"),
          identical(caller, .Random.seed))

regular <- data.frame(y = 0:16, lon = -120 + (0:16)/4, lat = 0)
tr <- spurtransform(~y, regular, lon = "lon", lat = "lat", transformation = "nn")
stopifnot(identical(tr$h_y, c(-1, rep(0, 15), 1)))
regular$y <- 1e12
tr <- spurtransform(~y, regular, lon = "lon", lat = "lat", transformation = "iso", radius = 1e6)
stopifnot(all(tr$h_y == 0), nrow(load_chetty_data()) == 741L)

stopifnot(identical(spatialSpur:::spatial_build_info()$package_version, as.character(packageVersion("spatialSpur"))))

dot_fit <- lm(y ~ ., d[, c("y", "x", "u")])
dot_scpc <- scpc(dot_fit, d, lon = "lon", lat = "lat", nthreads = 1L, cvs = TRUE)
stopifnot(length(coef(dot_scpc)) == 3L,
          all(dim(confint(dot_scpc, level = 0.90)) == c(3L, 2L)),
          any(grepl("avc", capture.output(summary(dot_scpc)))))
assert_error(scpc(lm(cbind(y, x) ~ u, d), d, lon = "lon", lat = "lat"), "single-response")

id_fit <- lm(y ~ x, d)
id_scpc <- scpc(id_fit, d[reorder, ], row_id = "id", lon = "lon", lat = "lat", nthreads = 1L)
stopifnot(identical(id_scpc$observations, d$id[!is.na(d$y)]))

if (requireNamespace("fixest", quietly = TRUE)) {
  sub_fit <- fixest::feols(y ~ x, d, subset = d$g > 0, demeaned = TRUE, data.save = TRUE, notes = FALSE)
  sub_scpc <- scpc(sub_fit, d, lon = "lon", lat = "lat", nthreads = 1L)
  stopifnot(identical(sub_scpc$observations, rownames(d)[fixest::obs(sub_fit)]))
  split_fits <- fixest::feols(y ~ x, d, subset = d$h <= 4L, split = ~g, demeaned = TRUE, data.save = TRUE, notes = FALSE)
  for (sf in split_fits) {
    sf_scpc <- scpc(sf, d, lon = "lon", lat = "lat", nthreads = 1L)
    stopifnot(identical(sf_scpc$observations, rownames(d)[fixest::obs(sf)]))
  }
  d_sing <- d
  d_sing$g_sing <- d_sing$g
  keep_idx <- which(d_sing$g_sing == 1L & !is.na(d_sing$y))[1L]
  d_sing$g_sing[d_sing$g_sing == 1L & seq_len(nrow(d_sing)) != keep_idx] <- 2L
  sing_fits <- fixest::feols(y ~ x | g_sing, d_sing, subset = d_sing$h <= 4L, split = ~h,
                             fixef.rm = "singleton", demeaned = TRUE, data.save = TRUE, notes = FALSE)
  sing_one <- sing_fits[[1L]]
  sing_scpc <- scpc(sing_one, d_sing, lon = "lon", lat = "lat", nthreads = 1L)
  stopifnot(identical(sing_scpc$observations, rownames(d_sing)[fixest::obs(sing_one)]))
  glm_fit <- fixest::feglm(I(y > 0) ~ x, d, family = "binomial", notes = FALSE)
  assert_error(scpc(glm_fit, d, lon = "lon", lat = "lat"), "feols")
  slope_fit <- fixest::feols(y ~ 1 | g[x], d, demeaned = TRUE, data.save = TRUE, notes = FALSE)
  assert_error(scpc(slope_fit, d, lon = "lon", lat = "lat", uncond = TRUE), "varying-slope")
}

set.seed(123)
before_null <- .Random.seed
invisible(do.call(spurtest_i0, c(list("x", d), within(opts, { seed <- NULL }))))
stopifnot(!identical(before_null, .Random.seed))

tr_expr <- spurtransform(log(abs(y) + 1) ~ I(x^2), d, lon = "lon", lat = "lat", nthreads = 1L)
d_manual <- transform(d, logy = log(abs(y) + 1), x2 = x^2)
tr_manual <- spurtransform(logy ~ x2, d_manual, lon = "lon", lat = "lat", nthreads = 1L)
expr_cols <- c("h_log(abs(y) + 1)", "h_I(x^2)")
stopifnot(all(expr_cols %in% colnames(tr_expr)),
          all(is.na(tr_expr[is.na(d$y), expr_cols[1L]])),
          isTRUE(all.equal(unname(as.matrix(tr_expr[expr_cols])),
                           unname(as.matrix(tr_manual[c("h_logy", "h_x2")])), tolerance = 1e-12)))
assert_error(spurtransform(y ~ x + offset(w), d, lon = "lon", lat = "lat"), "offset")

d1d <- data.frame(s = seq(0, 1, length.out = 24L), y = sin(seq(0.2, 2.7, length.out = 24L)))
d1d_rev <- d1d[rev(seq_len(nrow(d1d))), ]
set.seed(81)
z4 <- matrix(rnorm(4L * 25L), nrow = 4L)
t4 <- matrix(rnorm(4L * 80L), nrow = 4L)
d4 <- matrix(rnorm(4L * 80L), nrow = 4L)
zeros4 <- matrix(0, nrow = 4L, ncol = 25L)
for (fn in list(spurtest_i0, spurtest_i1)) {
  assert_error(fn("y", d1d, coords_euclidean = "s", q = 4L, draws = zeros4, tuning_draws = z4, nthreads = 1L),
               "positive column norms")
  assert_error(fn("y", d1d, coords_euclidean = "s", q = 4L, draws = z4, tuning_draws = zeros4, nthreads = 1L),
               "positive column norms")
  a_perm <- fn("y", d1d, coords_euclidean = "s", q = 4L, draws = d4, tuning_draws = t4, nthreads = 1L)
  b_perm <- fn("y", d1d_rev, coords_euclidean = "s", q = 4L, draws = d4, tuning_draws = t4, nthreads = 1L)
  stopifnot(abs(a_perm$lr - b_perm$lr) < 1e-11,
            abs(a_perm$pvalue - b_perm$pvalue) < 1e-12,
            abs(a_perm$ha_parm - b_perm$ha_parm) < 1e-11)
}
stopifnot(!is.null(colnames(p$transformed_vars)), !is.null(rownames(p$transformed_vars)))
w_out <- list()
tuned_25 <- withCallingHandlers(
  spurtest_i0("x", d, lon = "lon", lat = "lat", q = 4L, nrep = 25L, nthreads = 1L),
  warning = function(w) {
    w_out[[length(w_out) + 1L]] <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
stopifnot(any(grepl("unattainable", w_out)),
          is.numeric(tuned_25$best_tuning_power),
          is.finite(tuned_25$best_tuning_power))

cat("R installed integration checks passed\n")
