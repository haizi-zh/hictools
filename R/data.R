#' Gene density
#' 
#' Gene density data from `gencode.v30.b37.gene_density`. For a given interval,
#' indicate the number of genes overlapping with that interval.
#' 
#' The dataset comes with various interval widths: 50 kbp, 100 kbp, 250 kbp, 500
#' kbp, 1000 kbp and 2500 kbp.
#' 
#' To prepare the dataset:
#' 
#' `bedtools makewindows -g human_g1k_v37.chrom.sizes -w 500000 |`
#' `bedtools map -b gencode.v30.b37.gene_density.bedgraph -a - -c 4 -o count`
#' 
#' @examples 
#' 
#' data(gencode.v30.b37.gene_density.50kbp)
"gencode.v30.b37.gene_density.50kbp"

#' @rdname gencode.v30.b37.gene_density.50kbp
"gencode.v30.b37.gene_density.100kbp"

#' @rdname gencode.v30.b37.gene_density.50kbp
"gencode.v30.b37.gene_density.250kbp"

#' @rdname gencode.v30.b37.gene_density.50kbp
"gencode.v30.b37.gene_density.500kbp"

#' @rdname gencode.v30.b37.gene_density.50kbp
"gencode.v30.b37.gene_density.1000kbp"

#' @rdname gencode.v30.b37.gene_density.50kbp
"gencode.v30.b37.gene_density.250kbp"


#' GC content
#' 
#' Fixed-interval GC content for a number of reference genomes, with different bin sizes.
"gc.GRCh37.100kbp"

#' @rdname gc.GRCh37.100kbp
"gc.GRCh37.250kbp"

#' @rdname gc.GRCh37.100kbp
"gc.GRCh37.500kbp"

#' @rdname gc.GRCh37.100kbp
"gc.GRCh37.1000kbp"

#' @rdname gc.GRCh37.100kbp
"gc.GRCh37.2500kbp"

#' @rdname gc.GRCh37.100kbp
"gc.GRCh38.100kbp"

#' @rdname gc.GRCh37.100kbp
"gc.GRCh38.250kbp"

#' @rdname gc.GRCh37.100kbp
"gc.GRCh38.500kbp"

#' @rdname gc.GRCh37.100kbp
"gc.GRCh38.1000kbp"

#' @rdname gc.GRCh37.100kbp
"gc.GRCh38.2500kbp"