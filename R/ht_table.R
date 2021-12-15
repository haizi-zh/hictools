# Constructor of ht_table, the core data type for hictools
#
# @param dt A `data.table` or `data.frame` input.
# @param genome A character value specifying the genome name.
new_ht_table <-
  function(dt,
           resol,
           type = c("observed", "oe", "expected", "cofrag"),
           norm = c("NONE", "KR", "VC", "VC_SQRT"),
           genome = NULL) {
    assert_that(is(dt, "data.frame"))
    assert_that(is_null(genome) || is_character(genome))
    resol <- as.integer(resol)
    type <- match.arg(type)
    norm <- match.arg(norm)
    
    dt <- data.table::copy(dt) %>% data.table::as.data.table()
    dt <- dt[order(chrom1, pos1, chrom2, pos2)]

    data.table::setattr(dt, "class", c("ht_table", class(dt)))
    data.table::setkey(dt, "chrom1", "pos1", "chrom2", "pos2")
    data.table::setattr(dt, "resol", resol)
    data.table::setattr(dt, "type", type)
    data.table::setattr(dt, "norm", norm)
    data.table::setattr(dt, "genome", genome)

    dt
  }


validate_ht_table <- function(ht) {
  assert_that(!is_null(ht))
  assert_that(is(ht, "data.table"))

  ht_names <- names(ht)
  assert_that(are_equal(ht_names[1:5], c("chrom1", "pos1", "chrom2", "pos2", "score")))
  # Test column data types
  c(1, 3) %>%
    walk(function(idx) assert_that(is_character(ht[[idx]])))
  c(2, 4) %>%
    walk(function(idx) assert_that(is_integer(ht[[idx]])))
  assert_that(is_double(ht[[5]]))
  
  resol <- attr(ht, "resol")
  assert_that(is_valid_resol(resol))

  type <- attr(ht, "type")
  assert_that(is_scalar_character(type) && type %in% c("observed", "oe", "expected", "cofrag"))
  norm <- attr(ht, "norm")
  assert_that(is_scalar_character(norm) && norm %in% c("NONE", "KR", "VC", "VC_SQRT"))
  genome <- attr(ht, "genome")
  assert_that(is_null(genome) || is_scalar_character(genome))

  ht
}


#' Construct an `ht_table` object, which represents a Hi-C-like dataset.
#' 
#' @param dt A `data.frame` input.
#' @param resol A positive integer for the resolution (or bin size).
#' @param type `observed` for raw counts, `oe` for "observation/expected" value,
#'   and `cofrag` for cofragmentation contact score.
#' @param norm Indicate whether any normalization method has been applied.
#' @param genome A character scalar specifying the genome name. For example:
#'   `hg37-1kg`. Can be `NULL` if no genome is specified.
#' @export
ht_table <-
  function(dt,
           resol,
           type,
           norm,
           genome) {
    validate_ht_table(new_ht_table(dt, resol, type, norm, genome))
  }


#' @export
print.ht_table <- function(x, ...) {
  NextMethod()

  type <- attr(x, "type")
  norm <- attr(x, "norm")
  resol <- attr(x, "resol")
  genome <- attr(x, "genome") %||% "unspecified"

  cat("-------\n")
  cat(str_interp("Resolution: ${resol}\n"))
  cat(str_interp("Type: ${type}\n"))
  cat(str_interp("Norm: ${norm}\n"))
  cat(str_interp("Reference genome: ${genome}\n"))
}
