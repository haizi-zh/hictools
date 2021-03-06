% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{concat_hic}
\alias{concat_hic}
\title{Concatenate a list of \code{ht_table} objects}
\usage{
concat_hic(hic_list)
}
\value{
A concatenated \code{ht_table}
}
\description{
Hi-C data in the list should be compatible. In other words, they should have
the same \code{resol}, \code{norm}, \code{type}, and optionally, \code{genome}.
}
