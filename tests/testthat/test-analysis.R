example_path <- here("data-raw/example.bed")

test_that("Distance-based summary/stratification works", {
  strata <- load_hic_genbed(example_path) %>% strat_dist()
  expect_equal(nrow(strata), 378)

  # Is sorted
  index_1 <- strata %>% select(chrom, dist)
  index_2 <- index_1 %>% arrange(chrom, dist)
  expect_true(all(c(index_1 == index_2)))

  expect_equal(mean(strata$total), 76386.796, tolerance = 1e-3)

  strata <- load_hic_genbed(example_path) %>% strat_dist(smoothing = FALSE)
  expect_equal(mean(strata$total), 76388.87, tolerance = 1e-3)
})


test_that("Distance-based summary/stratification for selected chromosomes works", {
  chrom <- c("22", "20")
  strata <- load_hic_genbed(example_path, chrom = chrom) %>% strat_dist()
  expect_equal(strata$chrom %>% unique, sort(chrom))

  expect_equal(nrow(strata), 252)

  # Is sorted
  index_1 <- strata %>% select(chrom, dist)
  index_2 <- index_1 %>% arrange(chrom, dist)
  expect_true(all(c(index_1 == index_2)))

  expect_equal(mean(strata$total), 56696, tolerance = 1e-3)

  strata <- load_hic_genbed(example_path, chrom = chrom) %>% strat_dist(smoothing = FALSE)
  expect_equal(mean(strata$total), 56698, tolerance = 1e-3)
})
