## code to prepare `DATASET` dataset goes here

chrom_sizes <- read_tsv(
  "/Users/haizi/SynologyDrive/Research/CCHMC/dev/Juicebox/src/juicebox/tools/chrom/sizes/hg19.chrom.sizes",
  col_names = c("chrom", "size"),
  col_types = "ci"
) %>% mutate(ref_genome = "hg19")

usethis::use_data(chrom_sizes, overwrite = TRUE)
#
# chrom_sizes_hg19 <- chrom_sizes %>% select(-ref_genome)
#
# usethis::use_data(chrom_sizes_hg19, overwrite = TRUE, internal = TRUE)
