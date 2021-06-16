# Genearte white blood cell compartment tracks for various resolutions.
# Data source: wbc.rep1.mapq30.hic
# Method: Juicer tools

library(tidyverse)

devtools::load_all()

expand_grid(kbp = c(500, 1000, 2500, 5000), norm = c("NONE", "VC", "VC_SQRT", "KR")) %>% 
  pwalk(function(kbp, norm) {
    kbp <- as.integer(kbp)
    resol <- as.integer(kbp * 1e3)
    
    env <- rlang::env()
    track_name <-
      str_interp("gencode.v30.b37.gene_density.${kbp}kbp")
    data(list = track_name,
         package = "hictools",
         envir = env)
    gene_density <- env[[track_name]]
    
    comp <- c(1:22, "X", "Y") %>%
    # comp <- c("21", "22") %>%
      map(function(chrom) {
        logging::loginfo(str_interp("Processing ${kbp}kbp, ${norm}, chr${chrom}"))
        comp <-
          hictools::compartment_juicer(
            "data-raw/wbc.rep1.mapq30.hic",
            norm = norm,
            chrom = chrom,
            resol = resol
          )
        comp <-
          flip_compartment(comp,
                           standard = gene_density %>% bedtorch::as.bedtorch_table(),
                           resol = resol)
        comp
      }) %>% data.table::rbindlist()
    
    attr(comp, "genome") <- "hs37-1kg"
    comp %<>% bedtorch::as.GenomicRanges() %>% GenomicRanges::trim()
    track_name <- str_interp("wbc.rep1.compartment.hs37-1kg.${norm}.${kbp}kbp")
    env[[track_name]] <- comp
    
    save(
      list = track_name,
      file = str_interp("data/${track_name}.rda"),
      envir = env
    )
  })
