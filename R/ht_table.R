# Constructor of ht_table, the core data type for hictools
#
# @param dt A `data.table` or `data.frame` input.
# @param genome A character value specifying the genome name.
new_ht_table <-
  function(dt,
           resol,
           type = c("observed", "oe"),
           norm = c("NONE", "KR", "VC", "VC_SQRT")) {
    if (is(dt, "ht_table"))
      return(dt)

    supported_resol <- c(2.5e6L,
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
    stopifnot(!is.null(resol))
    resol <- as.integer(resol)
    stopifnot(resol %in% supported_resol)

    type <- match.arg(type)
    norm <- match.arg(norm)
    if (!is(dt, "data.table"))
      dt <- data.table::as.data.table(dt)

    dt <- dt[order(chrom1, pos1, chrom2, pos2)]

    data.table::setattr(dt, "class", c("ht_table", class(dt)))
    data.table::setkey(dt, "chrom1", "pos1", "chrom2", "pos2")
    data.table::setattr(dt, "resol", resol)
    data.table::setattr(dt, "type", type)
    data.table::setattr(dt, "norm", norm)
    data.table::setattr(dt, "chrom", unique(as.character(c(
      dt$chrom1, dt$chrom2
    ))))

    dt
  }


validate_ht_table <- function(ht) {
  stopifnot(!is.null(ht))
  stopifnot(is(ht, "data.table"))

  ht_names <- names(ht)
  stopifnot(identical(ht_names[1:5], c("chrom1", "pos1", "chrom2", "pos2", "score")))
  # Test column data types
  stopifnot(colnames(ht)[c(1, 3)] %>% sapply(function(x)
    is.character(ht[[x]])) %>% all())
  stopifnot(colnames(ht)[c(2, 4)] %>% sapply(function(x)
    is.integer(ht[[x]])) %>% all())
  stopifnot(colnames(ht)[5] %>% sapply(function(x)
    is.numeric(ht[[x]])) %>% all())

  chrom <- attr(ht, "chrom")
  stopifnot(!is.null(chrom) && is.character(chrom) && length(chrom) == 1)

  resol <- attr(ht, "resol")
  supported_resol <- c(2.5e6L,
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
  stopifnot(!is.null(resol) && is.integer(resol) && resol %in% supported_resol)

  type <- attr(ht, "type")
  stopifnot(!is.null(type) && type %in% c("observed", "oe"))

  norm <- attr(ht, "norm")
  stopifnot(!is.null(norm) && norm %in% c("NONE", "KR", "VC", "VC_SQRT"))

  ht
}


#' @export
ht_table <-
  function(dt,
           resol,
           type = c("observed", "oe"),
           norm = c("NONE", "KR", "VC", "VC_SQRT")) {
    validate_ht_table(new_ht_table(dt, resol, type, norm))
  }


#' @export
print.ht_table <- function(x, ...) {
  NextMethod()

  chrom <- paste(attr(x, "chrom"), collapse = " ")
  type <- attr(x, "type")
  norm <- attr(x, "norm")
  resol <- attr(x, "resol")

  cat("-------\n")
  cat(str_interp("Resolution: ${resol}\n"))
  cat(str_interp("Chrom: ${chrom}\n"))
  cat(str_interp("Type: ${type}\n"))
  cat(str_interp("Norm: ${norm}\n"))
}
