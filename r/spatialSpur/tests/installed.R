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
          identical(caller, .Random.seed))

regular <- data.frame(y = 0:16, lon = -120 + (0:16)/4, lat = 0)
tr <- spurtransform(~y, regular, lon = "lon", lat = "lat", transformation = "nn")
stopifnot(identical(tr$h_y, c(-1, rep(0, 15), 1)))
regular$y <- 1e12
tr <- spurtransform(~y, regular, lon = "lon", lat = "lat", transformation = "iso", radius = 1e6)
stopifnot(all(tr$h_y == 0), nrow(load_chetty_data()) == 741L)
cat("R installed integration checks passed\n")
