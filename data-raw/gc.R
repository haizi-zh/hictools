# Generate GC content dataset for common reference genomes

library(tidyverse)
library(GenomicRanges)

# Prepare Seqinfo objects
hs37d5 <- BSgenome::getBSgenome(genome = "BSgenome.Hsapiens.1000genomes.hs37d5", load.only = TRUE)

GRCh38 <- BSgenome::getBSgenome(genome = "BSgenome.Hsapiens.NCBI.GRCh38", load.only = TRUE)

# Make fixed-size windows
make_windows <- function(genome, width, chrom = NULL) {
  windows <- GenomicRanges::tileGenome(
    seqlengths = BSgenome::seqinfo(genome),
    tilewidth = width,
    cut.last.tile.in.chrom = TRUE
  )
  
  if (!is.null(chrom)) {
    windows <-
      windows[as.character(GenomicRanges::seqnames(windows)) %in% chrom]
  }
  
  return(windows)
}

# Calculate GC contents of each fragment
calc_gc <- function(genome, interval, mask_threshold = NULL) {
  # To reduce memory footprint, do the calculation for each chromosome separately
  result <-
    lapply(unique(GenomicRanges::seqnames(interval)), function(chrom) {
      interval <- interval[GenomicRanges::seqnames(interval) == chrom]
      seqs <- BSgenome::getSeq(genome, interval)
      gc <-
        (Biostrings::letterFrequency(seqs, letters = "CG"))[, 1]
      fourbases <-
        (Biostrings::letterFrequency(seqs, letters = "ATCG"))[, 1]
      gc <- gc / fourbases
      gc[is.nan(gc)] <- NA
      
      # Discard a bin if it contains too many Ns
      freq_n <- Biostrings::letterFrequency(seqs, letters = "N", as.prob = TRUE)[, 1]
      if (!is.null(mask_threshold)) {
        gc[freq_n > mask_threshold] <- NA
      }
      
      interval$gc <- gc
      interval$freq_n <- freq_n
      return(interval[!is.na(interval$gc)])
    })
  
  return(do.call(c, result))
}

# Generate the dataset
gc_result <- c(100e3L, 250e3L, 500e3L, 1000e3L, 2500e3L) %>%
  map(function(bin_size) {
    logging::loginfo(bin_size)
    calc_gc(hs37d5,
            make_windows(hs37d5, width = bin_size, chrom = c(1:22, "X")))
  })

gc.GRCh37.100kbp <- gc_result[[1]]
gc.GRCh37.250kbp <- gc_result[[2]]
gc.GRCh37.500kbp <- gc_result[[3]]
gc.GRCh37.1000kbp <- gc_result[[4]]
gc.GRCh37.2500kbp <- gc_result[[5]]
usethis::use_data(
  gc.GRCh37.100kbp,
  gc.GRCh37.250kbp,
  gc.GRCh37.500kbp,
  gc.GRCh37.1000kbp,
  gc.GRCh37.2500kbp,
  overwrite = TRUE
)

gc_result <- c(100e3L, 250e3L, 500e3L, 1000e3L, 2500e3L) %>%
  map(function(bin_size) {
    logging::loginfo(bin_size)
    calc_gc(GRCh38,
            make_windows(GRCh38, width = bin_size, chrom = c(1:22, "X")))
  })

gc.GRCh38.100kbp <- gc_result[[1]]
gc.GRCh38.250kbp <- gc_result[[2]]
gc.GRCh38.500kbp <- gc_result[[3]]
gc.GRCh38.1000kbp <- gc_result[[4]]
gc.GRCh38.2500kbp <- gc_result[[5]]
usethis::use_data(
  gc.GRCh38.100kbp,
  gc.GRCh38.250kbp,
  gc.GRCh38.500kbp,
  gc.GRCh38.1000kbp,
  gc.GRCh38.2500kbp,
  overwrite = TRUE
)
