example_hic_path <- "data-raw/example.hic"
example_short_path <- "data-raw/example.short"
example_dump_path <- "data-raw/example.txt"
example_bed_path <- "data-raw/example.bed"
example_cool_path <- "data-raw/example.cool"

test_that("Loading .hic file for single chromosome works", {
  resol <- 500e3
  hic <- load_juicer_hic(here(example_hic_path), chrom = "22", resol = resol)
  expect_true(check_valid_hic(hic))
  expect_identical(attr(hic, "resol"), as.integer(resol))
  expect_equal(attr(hic, "type"), "observed")
  expect_equal(attr(hic, "norm"), "NONE")
  expect_equal(nrow(hic), 2498)
  expect_equal(sum(hic$score), 2425325)
  expect_equal(sd(hic$score), 3843.397, tolerance = 0.001)

  type <- "oe"
  norm <- "KR"
  hic <- load_juicer_hic(here(example_hic_path), chrom = "22", resol = resol, matrix = type, norm = norm)
  expect_true(check_valid_hic(hic))
  expect_identical(attr(hic, "resol"), as.integer(resol))
  expect_equal(attr(hic, "type"), type)
  expect_equal(attr(hic, "norm"), norm)
  expect_equal(nrow(hic), 2498)
  expect_equal(sum(hic$score), 5019.03, tolerance = 0.01)
  expect_equal(sd(hic$score), 1.524, tolerance = 0.001)
})


test_that("Loading .hic requires valid resol", {
  expect_error(load_juicer_hic(example_hic_path, chrom = "22", resol = NULL))
  expect_error(load_juicer_hic(example_hic_path, chrom = "22", resol = 1L))
})


test_that("Loading .hic for missing chromosome works", {
  expect_warning(load_juicer_hic(example_hic_path, chrom = "X", resol = 500e3L))
})


test_that("Loading .hic only works for observed or oe matrix types", {
  expect_error(load_juicer_hic(example_hic_path, chrom = "22", resol = 500e3L, matrix = "unknown"))
})

test_that("Loading .hic only works for certain matrix types", {
  expect_error(load_juicer_hic(example_hic_path, chrom = "22", resol = 500e3L, norm = "unknown"))
})

test_that("Loading .hic for multiple chromosomes works", {
  resol <- 500e3
  chrom <- c("22", "21")
  hic <- load_juicer_hic(here(example_hic_path), chrom = chrom, resol = resol)
  expect_true(check_valid_hic(hic))
  expect_equal(attr(hic, "chrom"), c("21", "22"))
  expect_identical(attr(hic, "resol"), as.integer(resol))
  expect_equal(attr(hic, "type"), "observed")
  expect_equal(attr(hic, "norm"), "NONE")
  expect_equal(nrow(hic), 5173)
  expect_equal(sum(hic$score), 4906415)
  expect_equal(sd(hic$score), 3826.457, tolerance = 0.001)
})


test_that("Loading .short/.cool/GENBED file for single chromosome works", {
  chrom = "22"
  hic_list <- list(
    load_juicer_short(here(example_short_path), chrom = chrom, matrix = "observed", norm = "NONE"),
    load_hic_genbed(here(example_bed_path), chrom = chrom, matrix = "observed", norm = "NONE"),
    load_hic_cool(here(example_cool_path), chrom = chrom, matrix = "observed", norm = "NONE")
  )
  for (hic in hic_list) {
    expect_true(check_valid_hic(hic))
    expect_identical(attr(hic, "resol"), 500e3L)
    expect_equal(attr(hic, "type"), "observed")
    expect_equal(attr(hic, "norm"), "NONE")
    expect_equal(nrow(hic), 2498)
    expect_equal(sum(hic$score), 2425325)
    expect_equal(sd(hic$score), 3843.397, tolerance = 0.001)
    expect_equal(attr(hic, "chrom"), chrom)
  }
})

test_that("Loading .short/.cool/GENBED file for multiple chromosome works", {
  chrom <- c("22", "21")
  hic_list <- list(
    load_juicer_short(here(example_short_path), chrom = chrom, matrix = "observed", norm = "NONE"),
    load_hic_genbed(here(example_bed_path), chrom = chrom, matrix = "observed", norm = "NONE"),
    load_hic_cool(here(example_cool_path), chrom = chrom, matrix = "observed", norm = "NONE")
  )
  for (hic in hic_list) {
    expect_true(check_valid_hic(hic))
    expect_identical(attr(hic, "resol"), 500e3L)
    expect_equal(attr(hic, "type"), "observed")
    expect_equal(attr(hic, "norm"), "NONE")
    expect_equal(nrow(hic), 5173)
    expect_equal(sum(hic$score), 4906415)
    expect_equal(sd(hic$score), 3826.457, tolerance = 0.001)
    expect_equal(attr(hic, "chrom"), c("21", "22"))
  }
})


test_that("Loading Juicer-dumped file for single chromosome works", {
  chrom = "22"
  hic <- load_juicer_dump(here(example_dump_path), chrom = "22", matrix = "observed", norm = "NONE")
  expect_true(check_valid_hic(hic))
  expect_identical(attr(hic, "resol"), 500e3L)
  expect_equal(attr(hic, "type"), "observed")
  expect_equal(attr(hic, "norm"), "NONE")
  expect_equal(nrow(hic), 2498)
  expect_equal(sum(hic$score), 2425325)
  expect_equal(sd(hic$score), 3843.397, tolerance = 0.001)
  expect_equal(attr(hic, "chrom"), chrom)
})

test_that("Loading Juicer-dumped file must have exactly one chrom", {
  expect_error(load_juicer_dump(here(example_dump_path), chrom = NULL))
  expect_error(load_juicer_dump(here(example_dump_path), chrom = c("21", "22")))
})
