#' @importFrom readr read_lines read_delim
#' @importFrom stringr str_interp str_split str_detect str_sort str_match str_trim str_remove
#' @importFrom purrr map map_int map_dbl map_lgl map_chr map_dfr walk pwalk pmap_chr pmap_dfr `%||%` discard keep
#' @importFrom tidyr expand_grid
#' @importFrom dplyr `%>%` select filter rename mutate as_tibble tibble arrange desc
#' @importFrom magrittr `%<>%` set_attr
#' @importFrom rhdf5 h5read
#' @import ggplot2
#' @importFrom data.table setnames setkey fread fwrite shift setDT data.table rbindlist as.data.table
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps mcols `mcols<-` pintersect ranges `ranges<-` width
#' @importFrom GenomeInfoDb seqnames `seqnames<-` seqlengths seqinfo `seqinfo<-` seqlevels `seqlevels<-`
#' @importFrom S4Vectors from to queryHits subjectHits
#' @importFrom rlang is_null is_character is_integer is_double is_scalar_character is_scalar_integer is_scalar_logical
#' @importFrom assertthat assert_that are_equal
NULL
