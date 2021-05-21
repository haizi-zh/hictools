get_juicer_tools <- function() {
  system.file("extdata", "juicer_tools_1.22.01.jar", package = "hictools")
}


# Is the resolution valid? Only support a limited number of resolution choices
is_valid_resol <- function(resol) {
  allowed_resol <- c(2.5e6L,
                       1e6L,
                       500e3L,
                       250e3L,
                       100e3L,
                       50e3L,
                       25e3L,
                       10e3L,
                       5e3L,
                       2.5e3L,
                       1e3L)
  valid <- is.integer(resol) && resol %in% allowed_resol

  # NA is also not valid
  return(!is.na(valid) && valid)
}
