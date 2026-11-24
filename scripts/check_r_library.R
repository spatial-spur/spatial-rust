# Check registered routines, backend identity, and the dynamic exports of the installed R package.
packages <- commandArgs(trailingOnly = TRUE)
if (!length(packages)) packages <- "spatialSpur"

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
  for (sym in c("wrap__rust_run_model_scpc", "wrap__rust_run_spur_i0")) {
    stopifnot(is.loaded(sym, PACKAGE = pkg))
  }
  info <- .Call("wrap__rust_build_info", PACKAGE = pkg)
  stopifnot(info$linalg_provider %in% c("accelerate", "faer"))
  message(pkg, " uses ", info$linalg_provider, " (spec ", info$spec_version, ")")
  path <- getLoadedDLLs()[[pkg]][["path"]]
  actual <- symbols(path)
  expected <- paste0("R_init_", pkg)
  if (!identical(unname(actual), expected)) {
    stop(path, " exported ", length(actual), " symbols (expected only ", expected,
         "). First symbols: ", paste(head(actual, 10L), collapse = ", "))
  }
  message(pkg, " exports only ", expected)
  os <- Sys.info()[["sysname"]]
  if (os %in% c("Darwin", "Linux")) {
    flags <- if (os == "Darwin") "-u" else c("-D", "-u")
    undef <- system2("nm", c(flags, shQuote(path)), stdout = TRUE)
    if (!is.null(attr(undef, "status"))) stop("Could not inspect undefined symbols of ", path)
    if (os == "Darwin") undef <- sub("^_", "", trimws(undef))
    else undef <- trimws(sub("^[[:space:]]*U[[:space:]]+", "", grep("^[[:space:]]*U[[:space:]]+", undef, value = TRUE)))
    banned <- c("BODY", "CLOENV", "ENCLOS", "FORMALS", "R_BaseNamespace", "R_GetCurrentEnv", "R_MakeUnwindCont")
    found <- intersect(banned, undef)
    if (length(found)) {
      stop(path, " imports non-API R symbols: ", paste(found, collapse = ", "))
    }
    message(pkg, " imports zero forbidden R non-API symbols")
  }
}
