# Check registered routines, backend identity, and the dynamic exports of installed R packages.
packages <- commandArgs(trailingOnly = TRUE)
mode <- Sys.getenv("SPATIAL_RUST")
if (!length(packages) || !mode %in% c("0", "1")) {
  stop("Usage: set SPATIAL_RUST to 0 or 1 and run Rscript check_r_library.R PACKAGE ...")
}

symbols <- function(path) {
  os <- Sys.info()[["sysname"]]
  run <- function(tool, flags) {
    result <- system2(tool, c(flags, shQuote(path)), stdout = TRUE)
    if (!is.null(attr(result, "status"))) stop("Could not inspect ", path)
    result
  }
  if (os == "Darwin") {
    return(sub("^_", "", run("nm", "-gUj")))
  }
  if (os == "Windows") {
    tool <- Sys.which("llvm-readobj")
    if (!nzchar(tool)) tool <- file.path(Sys.getenv("ProgramFiles"), "LLVM/bin/llvm-readobj.exe")
    lines <- grep("^[[:space:]]*Name: ", run(tool, "--coff-exports"), value = TRUE)
    return(sub("^[[:space:]]*Name: ", "", lines))
  }
  sub("[[:space:]].*", "", run("nm", c("-D", "--defined-only", "--format=posix")))
}

for (pkg in packages) {
  loadNamespace(pkg)
  operation <- switch(pkg, scpcR = "wrap__rust_run_model_scpc", spuR = "wrap__rust_run_spur_i0",
                      stop("Unknown package ", pkg))
  native <- mode == "1"
  stopifnot(is.loaded(operation, PACKAGE = pkg) == native)
  if (native) {
    info <- .Call("wrap__rust_build_info", PACKAGE = pkg)
    stopifnot(info$linalg_provider %in% c("accelerate", "faer"))
    message(pkg, " uses ", info$linalg_provider)
  }
  path <- getLoadedDLLs()[[pkg]][["path"]]
  actual <- symbols(path)
  expected <- paste0("R_init_", pkg)
  if (!identical(unname(actual), expected)) {
    stop(path, " exported ", length(actual), " symbols (expected only ", expected,
         "). First symbols: ", paste(head(actual, 10L), collapse = ", "))
  }
  message(pkg, " exports only ", expected)
}
