% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{pearson_ht}
\alias{pearson_ht}
\title{Calculate Pearson}
\usage{
pearson_ht(hic_matrix, chrom, method = c("lieberman", "none"))
}
\arguments{
\item{method}{\code{lieberman} indicates to calculate Pearson correlation matrix
based on observed/expected values, using the method described in Lieberman
2009. \code{none} indicates directly calculate Pearson correlation from the Hi-C
matrix.}
}
\description{
Calculate Pearson
}
