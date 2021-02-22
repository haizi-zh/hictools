hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.observed.short")

hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.hic", resol = 500e3L, chrom = "22", matrix = "oe", norm = "NONE")

hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.observed.short") %>%
  strat_dist()

hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.observed.short") %>%
  oe_normalize(method = "juicer", juicertools = "~/CCHMC/tools/juicer_tools_1.22.01.jar") -> tmp_juicer_oe

hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.observed.short") %>%
  oe_normalize(method = "obs_exp") -> tmp2

cor(tmp1$score, tmp2$score, method = "spearman")

expand_grid(
  method = c("lieberman", "average", "obs_exp", "nonzero"),
  cor_type = c("spearman", "pearson")
) %>% pmap_dfr(function(method, cor_type) {
  hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.observed.short") %>%
    oe_normalize(method = method, smoothing = TRUE, min_nonzero = 10) %>%
    inner_join(
      x = .,
      y = tmp_juicer_oe,
      by = c("chr1", "pos1", "chr2", "pos2")
    ) %>% with(., cor(score.x, score.y, method = cor_type)) -> v
  tibble(method = method, cor_type = cor_type, value = v)
})


hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.observed.short") %>%
  hictools::calc_compartment(method = "juicer", juicertools = "~/CCHMC/tools/juicer_tools_1.22.01.jar") -> tmp_comp_juicer

hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.observed.short") %>%
  calc_compartment(method = "pca", oe_method = "lieberman", matrix = "observed") %>%
  inner_join(
    x = .,
    y = tmp_juicer,
    by = c("chrom", "start", "end")
  ) %>% na.omit() %>% with(., cor(score.x, score.y, method = "pearson"))

expand_grid(
  oe_method = c("lieberman", "average", "obs_exp", "nonzero"),
  cor_type = c("spearman", "pearson")
) %>% pmap_dfr(function(oe_method, cor_type) {
  hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.observed.short") %>%
    calc_compartment(method = "pca", oe_method = oe_method, matrix = "observed", smoothing = TRUE, min_nonzero = 4) %>%
    inner_join(
      x = .,
      y = tmp_juicer,
      by = c("chrom", "start", "end")
    ) %>% with(., cor(score.x, score.y, method = cor_type, use = "complete.obs")) -> v
  tibble(oe_method = oe_method, cor_type = cor_type, value = v)
})

hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.oe.txt", format = "juicer_dump", chrom = "22") %>%
  # oe_normalize(method = "juicer", juicertools = "~/CCHMC/tools/juicer_tools_1.22.01.jar") %>%
  calc_compartment(method = "pca", oe_method = "lieberman", matrix = "oe") %>%
  inner_join(
    x = .,
    y = tmp_juicer,
    by = c("chrom", "start", "end")
  ) %>% na.omit() %>% with(., cor(score.x, score.y, method = "pearson"))


# Stratify & smoothing
hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.observed.short") %>%
  strat_dist() %>%
  arrange(nonzero) -> tmp_strat


# Visualization
hictools::load_hic(file_path = "../../cofrag/sandbox/wbc.rep1.chr22.observed.short") %>%
  calc_compartment(method = "juicer", juicertools = "~/CCHMC/tools/juicer_tools_1.22.01.jar") %>%
  na.omit() %>%
  plot_compartment()
