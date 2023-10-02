.onLoad <- function(libname, pkgname) {
  if (is.na(Sys.getenv("COLUMNS", unset = NA))) {
    Sys.setenv(COLUMNS = getOption("width"))
  }

  cores <- detectCores()
  cl <- makeSOCKcluster(round(cores / 2))
  registerDoSNOW(cl)
  assign("cl", cl, envir = parent.env(environment()))

  UKMotor <- ChainLadder::UKMotor
  assign("UKMotor", UKMotor, envir = parent.env(environment()))

  key <- list(
    single = 1L,
    calendar = 2L,
    origin = 3L,
    normal = 1L,
    gamma = 2L,
    poisson = 3L,
    parametric = 1L,
    residuals = 2L,
    pairs = 3L,
    standardised = 1L,
    studentised = 2L,
    lognormal = 3L,
    conditional = 1L,
    unconditional = 2L
  )

  assign("key", key, envir = parent.env(environment()))
}

.onUnload <- function(libpath) {
  stopCluster(claimsBoot:::cl)
}

.onDetach <- function(libpath) {
  stopCluster(claimsBoot:::cl)
}