get_juicer_tools <- function() {
  system.file("extdata", "juicer_tools_1.22.01.jar", package = "hictools")
}


# Is the resolution valid? Only support a limited number of resolution choices
is_valid_resol <- function(resol) {
  allowed_resol <- as.integer(c(5, 10, 25, 50, 100, 250, 500, 1000, 2500) * 1e3)
  valid <- is.integer(resol) && resol %in% allowed_resol

  # NA is also not valid
  return(!is.na(valid) && valid)
}
