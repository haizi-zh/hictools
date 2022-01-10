#' Return the path of juicer_tools which is shipped along with hictools
get_juicer_tools <- function() {
  system.file("extdata", "juicer_tools_1.22.01.jar", package = "hictools")
}


#' Is the resolution valid? Only support a limited number of resolution choices
is_valid_resol <- function(resol) {
  return(is_scalar_integer(resol) && resol %in% allowed_resol())
}


#' Get a list of all allowed resolutions
allowed_resol <- function() {
  return(c(
    2.5e6L,
    1e6L,
    500e3L,
    250e3L,
    100e3L,
    50e3L,
    25e3L,
    10e3L,
    5e3L,
    2.5e3L,
    1e3L
  ))
}


#' Concatenate a list of `ht_table`
#' 
#' Under the hood, an `ht_table` is a `data.table`, so you can
#' concatenate a list of `ht_table` by `data.tables::rbindlist`.
#' However, the result is no long an `ht_table`, but a `data.table`. In order to
#' preserve `ht_table` properties, `concat_hic` is preferable.
#' 
#' Hi-C data in the list should be compatible. In other words, they should have
#' the same `resol`, `norm`, `type`, and optionally, `genome`.
#'
#' @param hic_list A list of `ht_table`
#' 
#' @return A concatenated `ht_table`
#' 
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
  
  data.table::rbindlist(l = hic_list) %>%
    ht_table(
      resol = resol,
      type = type,
      norm = norm,
      genome = genome
    )
}

#' Set Hi-C normalization method
#' 
#' @param hic a `ht_table` object
#' @export
`hic_norm<-` <- function(hic, value = c("NONE", "KR", "VC", "VC_SQRT", "SCALE")) {
  assert_that(is(hic, "ht_table"))
  value <- match.arg(value)
  
  data.table::setattr(hic, "norm", value)
  return(hic)
}

#' Get Hi-C normalization method
#' 
#' @param hic a `ht_table` object
#' @export
hic_norm <- function(hic) {
  assert_that(is(hic, "ht_table"))
  return(attr(hic, "norm"))
}

#' Get Hi-C data type
#' 
#' @param hic a `ht_table` object
#' @export
hic_type <- function(hic) {
  assert_that(is(hic, "ht_table"))
  return(attr(hic, "type"))
}

#' Set Hi-C data type
#' 
#' @param hic a `ht_table` object
#' @export
`hic_type<-` <- function(hic, value = c("observed", "oe", "expected", "pearson", "cofrag")) {
  assert_that(is(hic, "ht_table"))
  value <- match.arg(value)
  
  data.table::setattr(hic, "type", value)
  return(hic)
}


#' Get Hi-C genome
#' 
#' @param hic a `ht_table` object
#' @export
hic_genome <- function(hic) {
  assert_that(is(hic, "ht_table"))
  return(attr(hic, "genome"))
}


#' Set Hi-C genome
#' 
#' @param hic a `ht_table` object
#' @export
`hic_genome<-` <- function(hic, value) {
  assert_that(is(hic, "ht_table"))
  # Make sure the genome name is valid
  assert_that(is_scalar_character(value))
  genome <- bedtorch::get_seqinfo(value)
  assert_that(!is.null(genome))
  
  assert_that(
    unique(c(hic$chrom1, hic$chrom2)) %in% GenomicRanges::seqnames(genome),
    msg = str_interp("Hi-C records are not compatible with genome ${value}")
  )
  
  data.table::setattr(hic, "genome", value)
  return(hic)
}


#' Get Hi-C resolution
#' 
#' @param hic a `ht_table` object
#' @export
hic_resol <- function(hic) {
  assert_that(is(hic, "ht_table"))
  return(as.integer(attr(hic, "resol")))
}


#' Get Hi-C sample name
#' 
#' @param hic a `ht_table` object
#' @export
hic_sample <- function(hic) {
  assert_that(is(hic, "ht_table"))
  return(attr(hic, "sample"))
}


#' Set Hi-C sample name
#' 
#' @param hic a `ht_table` object
#' @export
`hic_sample<-` <- function(hic, value) {
  assert_that(is(hic, "ht_table"), is_scalar_character(value))
  
  data.table::setattr(hic, "sample", value)
  return(hic)
}



#' Apply a "decay curve" to a Hi-C matrix
#' 
#' This may be useful when `hic_matrix` is a cofrag dataset
#'
#' @param decay_curve A data frame with two columns: `dist` and `factor`
#' @export
decay <- function(hic_matrix, decay_curve) {
  assert_that(is(decay_curve, "data.frame"))
  assert_that(is(hic_matrix, "data.frame"))
  
  hic <- data.table::copy(hic_matrix)
  decay_curve <- data.table::as.data.table(decay_curve)
  decay_curve <- decay_curve[, .(dist, factor)]
  hic[, dist := abs(pos2 - pos1)]
  hic <- merge(hic, decay_curve, all.x = TRUE, by = "dist")[!is.na(factor)]
  hic[, `:=`(score = score * factor, dist = NULL, factor = NULL)]
  return(ht_table(hic, copy_from = hic_matrix))
}


#' Test if `x` is a vector of whole numbers
is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (!is.numeric(x) || length(x) == 0)
    return(FALSE)
  
  return(abs(x - round(x)) < tol)
}