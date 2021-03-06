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


#' Concatenate a list of `ht_table` objects
#' 
#' Hi-C data in the list should be compatible. In other words, they should have
#' the same `resol`, `norm`, `type`, and optionally, `genome`.
#' @return A concatenated `ht_table`
#' @export
concat_hic <- function(hic_list) {
  assert_that(rlang::is_list(hic_list))
  hic_list %>% walk( ~ assert_that(is(., "ht_table")))
  
  if (length(hic_list) == 0)
    return(hic_list)
  
  # Should be compatible
  resol <- hic_list %>% map_int( ~ attr(., "resol")) %>% unique()
  assert_that(is_scalar_integer(resol))
  norm <- hic_list %>% map_chr( ~ attr(., "norm")) %>% unique()
  assert_that(is_scalar_character(norm))
  type <- hic_list %>% map_chr( ~ attr(., "type")) %>% unique()
  assert_that(is_scalar_character(type))
  # genome can be NULL
  genome <-
    hic_list %>% map( ~ attr(., "genome")) %>% unlist() %>% unique()
  assert_that(is_null(genome) || is_scalar_character(genome))
  
  data.table::rbindlist(l = hic_list) %>% hictools::ht_table(
    resol = resol,
    type = type,
    norm = norm,
    genome = genome
  )
}