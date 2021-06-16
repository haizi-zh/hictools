# Genearte gene density tracks for various interval widths
# Data source: gencode.v30.b37.gene_density.bedgraph

library(tidyverse)

load("data/chrom_sizes.hs37-1kg.rda")

chrom_sizes_file <- tempfile(fileext = ".tsv")
data.table::fwrite(`chrom_sizes.hs37-1kg`, file = chrom_sizes_file, sep = "\t", col.names = FALSE)

c(50, 100, 250, 500, 1000, 2500, 5000) %>%
  walk(function(kbp) {
    track_bed <- tempfile(fileext = ".bed")
    on.exit(rm(track_bed), add = TRUE)
    cmd <-
      str_interp(
        paste0(
          "bedtools makewindows -g ${chrom_sizes_file} ",
          "-w ${kbp}000 | bedtools map -b ",
          "data-raw/gencode.v30.b37.gene_density.bedgraph ",
          "-a - -c 4 -o count > ${track_bed}"))
    system(cmd)
    
    track_name <- str_interp("gencode.v30.b37.gene_density.${kbp}kbp")
    
    env <- rlang::env()
    env[[track_name]] <- bedtorch::read_bed(track_bed, genome = "hs37-1kg")
    env$gencode.v30.b37.gene_density.500kbp
    
    save(list = track_name, file = str_interp("data/gencode.v30.b37.gene_density.${kbp}kbp.rda"), envir = env)
  })
