library(spatialSpur)

d <- read.csv("common.csv")
q <- 6L
nrep <- 400L
read_draws <- function(file, rows) {
  matrix(readBin(file, "double", n = rows * nrep, size = 8L, endian = "little"),
         nrow = rows, byrow = TRUE)
}
normals <- read_draws("diagnostics.bin", 2L * q)
geo <- list(lon = "lon", lat = "lat", nthreads = 1L, backend = "rust")
p <- do.call(spur, c(list(y ~ x + x2, d, q = q, draws = normals[1:q, ],
                         tuning_draws = normals[q + (1:q), ], cvs = TRUE), geo))
stopifnot(p$levels$backend == "rust")
fit_rows <- function(fit, phase, case) {
  b <- cbind(fit$scpcstats, fit$critical_95)
  labels <- sub("^\\(Intercept\\)$", "_cons", rownames(b))
  data.frame(case = case, stage = phase, term = labels, n = fit$n_obs, q = fit$q,
             coef = b[, 1], se = b[, 2], t = b[, 3], p = b[, 4],
             lo = b[, 5], hi = b[, 6], cv = b[, 7])
}
cl <- read.csv("clustered.csv")
cl_fit <- do.call(scpc, c(list(lm(y ~ x, cl), cl, cluster = "cl"), geo))
fits <- list(
  fit_rows(p$levels, "levels", "common"),
  fit_rows(p$transformed, "adjusted", "common"),
  fit_rows(cl_fit, "levels", "clustered")
)
cv_rows <- function(fit, phase) {
  m <- fit$scpccvs
  labels <- sub("^\\(Intercept\\)$", "_cons", rownames(m))
  data.frame(stage = phase, term = labels, cv32 = m[, 1], cv10 = m[, 2], cv05 = m[, 3], cv01 = m[, 4])
}
write.csv(rbind(cv_rows(p$levels, "levels"), cv_rows(p$transformed, "adjusted")), "cvs-r.csv", row.names = FALSE)
lh <- do.call(spurhalflife, c(list("y", d, q = q, draws = read_draws("halflife.bin", q),
                                  normdist = TRUE), geo))
lh_raw <- do.call(spurhalflife, c(list("y", d, q = q, draws = read_draws("halflife.bin", q),
                                      normdist = FALSE), geo))
diagnostic <- do.call(rbind, lapply(p$tests, function(test) {
  data.frame(test = test$test, statistic = test$LR, p = test$pvalue, ha = test$ha_parm,
             cv01 = test$cvalue[1], cv05 = test$cvalue[2], cv10 = test$cvalue[3])
}))
write.csv(diagnostic, "diagnostics-r.csv", row.names = FALSE)
write.csv(data.frame(lh$components), "halflife-r.csv", row.names = FALSE)
write.csv(data.frame(lh_raw$components), "halflife_raw-r.csv", row.names = FALSE)
tr_nn <- do.call(spurtransform, c(list(~y, d, transformation = "nn"), geo))
tr_iso <- do.call(spurtransform, c(list(~y, d, transformation = "iso", radius = 500000), geo))
tr_cl <- spurtransform(~y, cl, transformation = "cluster", clustvar = "cl", nthreads = 1L, backend = "rust")
write.csv(data.frame(h_y = tr_nn$h_y), "tr_nn_r.csv", row.names = FALSE)
write.csv(data.frame(h_y = tr_iso$h_y), "tr_iso_r.csv", row.names = FALSE)
write.csv(data.frame(h_y = tr_cl$h_y), "tr_cl_r.csv", row.names = FALSE)

if (file.exists("chetty-cases.csv")) {
  cases <- read.csv("chetty-cases.csv")
  for (key in unique(cases$term)) {
    frame <- cases[cases$term == key, ]
    coords <- as.matrix(frame[c("lon", "lat")])
    a <- scpc(lm(y ~ 0 + x, frame), coords = coords, geographic = TRUE, nthreads = 1L)
    tr <- do.call(spurtransform, c(list(y ~ x, frame), geo))
    b <- scpc(lm(h_y ~ 0 + h_x, tr), coords = coords, geographic = TRUE, nthreads = 1L)
    rownames(b$scpcstats) <- "x"
    fits <- c(fits, list(fit_rows(a, "levels", key), fit_rows(b, "adjusted", key)))
  }
}
write.csv(do.call(rbind, fits), "scpc-r.csv", row.names = FALSE)
cat("R parity completed\n")
