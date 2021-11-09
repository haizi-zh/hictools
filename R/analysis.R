# MIT License
#
# Copyright (c) 2021 Haizi Zheng
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Author: Haizi Zheng
# Copyright: Copyright 2021, Haizi Zheng
# Email: haizi.zh@gmail.com
# License: MIT
#
# Provides functions for analyzing Hi-C dataset

#' Calculate distance-stratified summary
#'
#' @param hic_matrix A Hi-C object
#' @param smoothing Whether to apply smoothing for large-distance strata.
#' Default is \code{TRUE}.
#' @param min_nonzero If smoothing, in each stratum, the minimal number of
#' non-zero entries. Default is 10.
#' @export
strat_dist <-
  function(hic_matrix,
           smoothing = TRUE,
           min_nonzero = 10) {
    resol <- attr(hic_matrix, "resol")
    chroms <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2))

    hic_matrix <- as_tibble(hic_matrix)

    chroms %>% map_dfr(function(chrom) {
      strat_data <- hic_matrix %>%
        mutate(dist = abs(pos1 - pos2)) %>%
        group_by(dist) %>%
        summarize(
          # avg = mean(score, na.rm = TRUE),
          total = sum(score, na.rm = TRUE),
          nonzero = sum(score > 0)
        ) %>%
        mutate(distbin = dist %/% resol) %>%
        arrange(distbin)

      if (smoothing) {
        strat_large <- strat_data %>% filter(nonzero >= min_nonzero)
        strat_small <- strat_data %>% filter(nonzero < min_nonzero)

        # For each row in strat_small, try to expand the window and do smoothing
        strat_smoothed <- strat_small %>%
          pmap_dfr(function(...) {
            row <- tibble(...)
            distbin <- row$distbin

            distbin_1 <- distbin_2 <- distbin
            expanded <- row
            # Expand until there are enough non-zero entries
            while (TRUE) {
              distbin_1 <- distbin_1 - 1
              distbin_2 <- distbin_2 + 1
              expanded <- strat_data %>%
                filter(between(distbin, distbin_1, distbin_2))
              if (sum(expanded$nonzero) > 10) {
                break
              } else if (distbin_1 < 0) {
                break
              }
            }
            row$total <-
              sum(expanded$total) / sum(expanded$nonzero) * row$nonzero
            # row$avg <- row$total / row$nonzero
            row
          })
        strat_data <- bind_rows(strat_smoothed, strat_large)
      }
      strat_data %>% mutate(chrom = chrom) %>% select(chrom, everything())
    }) %>%
      arrange(chrom, dist)
  }


#' Calculate oe normalized matrix using Juicer tools
#'
#' @param hic_matrix A Hi-C object
#' @export
oe_juicer <-
  function(hic_matrix,
           juicertools,
           java = "java",
           norm = "NONE",
           ref_genome = "hg19") {
    chroms <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2))
    resol <- attr(hic_matrix, "resol")

    hic_matrix <- as_tibble(hic_matrix)

    chroms %>% map_dfr(function(chrom) {
      temp_hic <- tempfile(fileext = ".hic")
      hic_matrix %>%
        filter(chrom1 == chrom & chrom2 == chrom) %>%
        write_juicer_hic(
          file_path = temp_hic,
          juicertools = juicertools,
          java = java,
          ref_genome = ref_genome
        )
      oe_matrix <-
        load_juicer_hic(
          file_path = temp_hic,
          chrom = chrom,
          resol = resol,
          matrix = "oe",
          norm = norm
        )
      unlink(temp_hic)
      oe_matrix
    })
  }


get_possible_dist <-
  function(resol,
           genome_wide = FALSE,
           chrom = NULL,
           ref_genome = c("hs37-1kg", "b37", "hs37d5", "hg19", "GRCh37")) {
    resol <- as.integer(resol)
    assert_that(is_scalar_integer(resol) && resol > 0)
    assert_that(is_scalar_logical(genome_wide))
    
    # Currently only support hs37-1kg and similar reference genomes
    ref_genome <- match.arg(ref_genome)
    assert_that(ref_genome %in% c("hs37-1kg", "hs37d5", "b37"))
    ref_genome <- "hs37-1kg"
    
    chrom_sizes <- local({
      env <- environment()
      dataset_name <- paste0("chrom_sizes.", ref_genome)
      data(list = dataset_name, package = "hictools", envir = env)
      env[[dataset_name]]
    })
    chrom_sizes[, n_bins := as.integer((size %/% resol) + 1)]
    chrom_sizes %<>% as_tibble()
    
    if (genome_wide)
      # Genome-wide
      chrom <- NULL

    if (is_null(chrom))
      chrom <- chrom_sizes$chrom

    result <- chrom %>% map_dfr(function(chrom) {
      # Chromosome-wide
      chrom_sizes %<>% filter(chrom == !!chrom)

      chrom_sizes %>% pmap_dfr(function(chrom, n_bins, ...) {
        distbin <- 0:(n_bins - 1)
        pd <- n_bins - distbin
        tibble(
          chrom = chrom,
          n_bins = n_bins,
          distbin = distbin,
          possible_dist = pd
        )
      }) %>%
        group_by(chrom, distbin) %>%
        summarise(possible_dist = sum(possible_dist),
                  .groups = "drop")
    })

    if (genome_wide) {
      result %<>%
        group_by(distbin) %>%
        summarize(possible_dist = sum(possible_dist),
                  .groups = "drop")
    }
    
    data.table::as.data.table(result)
  }



# Calculate the observed/expected matrix
#
#' @export
oe_ht <- function(hic_matrix,
                  method = c("lieberman", "obs_exp", "nonzero", "average"),
                  genome_wide = FALSE,
                  smoothing = TRUE,
                  min_nonzero = 4L) {
  assert_that(is(hic_matrix, "ht_table"))
  method <- match.arg(method)
  assert_that(is_scalar_logical(genome_wide))
  assert_that(is_scalar_logical(smoothing))
  min_nonzero <- as.integer(min_nonzero)
  assert_that(is_scalar_integer(min_nonzero) && min_nonzero >= 0)

  chroms <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2) %>% as.character())
  resol <- attr(hic_matrix, "resol")
  norm <- attr(hic_matrix, "norm")

  hic_matrix <- as_tibble(hic_matrix)

  chroms %>% map_dfr(function(chrom) {
    # Ignore inter-chromosomal contacts
    hic_matrix %<>%
      filter(chrom1 == chrom & chrom2 == chrom) %>%
      mutate(distbin = abs(pos1 - pos2) %/% resol)
    
    max_bin_idx <-
      max(c(hic_matrix$pos1, hic_matrix$pos2)) %/% resol + 1
    min_bin_idx <-
      min(c(hic_matrix$pos1, hic_matrix$pos2)) %/% resol + 1
    weight <-
      strat_dist(hic_matrix,
                 smoothing = smoothing,
                 min_nonzero = min_nonzero)
    
    joined <- inner_join(x = hic_matrix,
                         y = weight,
                         by = "distbin")
    
    (if (method == "obs_exp") {
      joined %>%
        mutate(observed = score, score = score / (total / (max_bin_idx - min_bin_idx + 1 - distbin)))
    } else if (method == "average") {
      joined %>%
        mutate(observed = score, score = score / avg)
    } else if (method == "nonzero") {
      joined %>%
        mutate(observed = score,
               score = score / (total / nonzero))
    } else if (method == "lieberman") {
      if (genome_wide) {
        possible_dist <-
          get_possible_dist(resol = resol, genome_wide = TRUE)
      } else {
        possible_dist <-
          get_possible_dist(resol = resol,
                            genome_wide = FALSE,
                            chrom = chrom) %>%
          filter(chrom == !!chrom)
      }
      inner_join(x = joined,
                 y = possible_dist,
                 by = "distbin") %>%
        mutate(observed = score,
               score = score / (total / possible_dist))
    } else {
      stop(str_interp("Invalid method ${method}"))
    }) %>%
      select(chrom1:score, observed) %>%
      set_attr(which = "resol", resol)
  }) %>%
    hictools::ht_table(
      resol = resol,
      type = "oe",
      norm = norm,
      genome = attr(hic_matrix, "genome")
    )
}


#' Calculate compartment scores using Juicer tools
#'
#' @param hic_matrix Either a `hic_matrix` object, or the path of the `.hic`
#'   name. If a file path, `resol` cannot be NULL since we can't infer
#'   resolution directly from the `.hic` file.
#' @param chrom A character vector indicating which chromosomes to run. If
#'   `NULL`, will calculate compartment scores for all chromosomes. However, in
#'   this case, `hic_matrx` cannot be the `.hic` file name.
#' @export
get_compartment <- function(hic_matrix,
                            method = c("juicer", "lieberman", "obs_exp", "nonzero", "average"),
                            chrom = NULL,
                            standard = NULL,
                            smooth = NULL,
                            # juicertools = get_juicer_tools(),
                            # java = "java",
                            # norm = c("NONE", "VC", "VC_SQRT", "KR", "SCALE"),
                            ...) {
  UseMethod("get_compartment")
}

#' @export
get_compartment.ht_table <- function(hic_matrix,
                                     method = c("juicer", "lieberman", "obs_exp", "nonzero", "average", "fanc"),
                                     chrom = NULL,
                                     standard = NULL,
                                     smooth = NULL,
                                     ...) {
  method <- match.arg(method)
  
  if (is_null(chrom))
    chrom <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2))
  else
    chrom <- as.character(chrom)
  assert_that(is_character(chrom) && length(chrom) >= 1)
  
  resol <- attr(hic_matrix, "resol")
  assert_that(is_scalar_integer(resol) &&
                resol %in% c(500e3L, 1000e3L, 2500e3L, 5000e3L))
  
  assert_that(is_null(standard) ||
                is(standard, "GRanges") ||
                is(standard, "data.frame"))
  smooth <-
    if (is_null(smooth))
      NULL
  else
    as.integer(smooth)
  assert_that(is_null(smooth) || smooth >= 1)
  
  genome <- attr(hic_matrix, "genome")
  assert_that(is_null(genome) || is_scalar_character(genome))
  
  if (method %in% c("juicer", "fanc")) {
    # Save hic to a temporary .hic file, and call juicer tools to get compartments
    hic_file <- tempfile(fileext = ".hic")
    on.exit(file.remove(hic_file), add = TRUE)
    
    args <- list(...)
    java <- args$java %||% "java"
    juicertools <- args$juicertools %||% get_juicer_tools()
    
    if (genome %in% c("hs37-1kg", "hg19", "GRCh37"))
      juicer_genome <- "hg19"
    else if (genome %in% c("hg38", "GRCh38"))
      juicer_genome <- "hg38"
    else
      stop(genome)
    
    hic_matrix %>% write_hic(
      file_path = hic_file,
      format = "juicer_hic",
      juicertools = juicertools,
      java = java,
      ref_genome = juicer_genome
    )
    
    comps <- get_compartment.character(
      hic_matrix = hic_file,
      method = method,
      chrom = chrom,
      standard = standard,
      smooth = 1L,
      resol = attr(hic_matrix, "resol") %||% guess_resol(hic_matrix),
      genome = attr(hic_matrix, "genome") %||% args$genome,
      juicertools = juicertools,
      java = java,
      norm = args$norm %||% "NONE"
    )
  } else {
    comps <- compartment_ht(
      hic_matrix = hic_matrix,
      method = method,
      chrom = chrom,
      npc = 2L,
      standard = standard
    )
  }
  
  if (isTRUE(smooth > 1))
    comps$score <- zoo::rollmean(
      comps$score,
      k = smooth,
      na.pad = TRUE,
      na.rm = TRUE,
      align = "center"
    )
  
  return(comps)
}


#' @export
get_compartment.character <- function(hic_matrix,
                                      method = c("juicer", "lieberman", "obs_exp", "nonzero", "average", "fanc"),
                                      chrom = NULL,
                                      standard = NULL,
                                      smooth = NULL,
                                      resol,
                                      genome = NULL,
                                      juicertools = get_juicer_tools(),
                                      java = "java",
                                      norm = c("NONE", "VC", "VC_SQRT", "KR", "SCALE")) {
  
  method <- match.arg(method)
  norm <- match.arg(norm)
  
  if (isTRUE(method == "juicer")) {
    # No need to load the .hic file. Instead, directly call compartments using Juicer tools
    comps <- compartment_juicer_file(
      hic_file = hic_matrix,
      juicertools = juicertools,
      java = java,
      norm = norm,
      chrom = chrom,
      resol = resol,
      standard = standard,
      genome = genome
    )
  } else if (method == "fanc") {
    # No need to load the .hic file. Instead, directly call compartments using Juicer tools
    comps <- compartment_fanc(
      hic_file = hic_matrix,
      norm = norm,
      resol = resol,
      genome = genome
    )
  } else {
    hic_matrix <-
      load_juicer_hic(
        file_path = hic_matrix,
        resol = resol,
        chrom = chrom,
        type = "observed",
        norm = norm, 
        genome = genome
      )
    comps <- get_compartment(
      hic_matrix = hic_matrix,
      method = method,
      chrom = chrom,
      standard = standard,
      smooth = 1L
    )
  }
  
  if (isTRUE(smooth > 1)) {
    comps$score <- zoo::rollmean(
      comps$score,
      k = smooth,
      na.pad = TRUE,
      na.rm = TRUE,
      align = "center"
    )
  }
  
  
  return(comps)
}


#' Get compartment using Juicer tools, from an existing .hic file
compartment_juicer_file <- function(hic_file,
                                    juicertools = get_juicer_tools(),
                                    java = "java",
                                    norm = c("NONE", "VC", "VC_SQRT", "KR", "SCALE"),
                                    chrom = NULL,
                                    resol = NULL,
                                    standard = NULL,
                                    genome = NULL) {
  assert_that(is_scalar_character(hic_file) && endsWith(hic_file, ".hic"))
  assert_that(is_scalar_character(java))
  assert_that(is_scalar_character(juicertools))
  norm <- match.arg(norm)
  chrom <- as.character(chrom)
  assert_that(length(chrom) >= 1)
  resol <- as.integer(resol)
  assert_that(is_scalar_integer(resol) && resol %in% c(500e3L, 1000e3L, 2500e3L, 5000e3L))
  assert_that(is_null(standard) || is(standard, "GRanges") || is(standard, "data.frame"))
  assert_that(is_null(genome) || is_scalar_character(genome))
  
  comps <- chrom %>% map(function(chrom) {
    ev_file <- tempfile()
    on.exit(unlink(ev_file), add = TRUE)
    cmd <-
      str_interp(
        "${java} -jar ${juicertools} eigenvector ${norm} ${hic_file} ${chrom} BP ${resol} ${ev_file}"
      )
    logging::loginfo(cmd)
    retcode <- system(cmd)
    assertthat::assert_that(retcode == 0,
                            msg = str_interp("Error in compartment calculation , RET: ${retcode}, CMD: ${cmd}"))
    
    comps <- data.table::fread(ev_file, col.names = "score", na.strings = c("", "NA", "NaN"))
    comps[, {
      start <- as.integer(0:(length(score) - 1) * resol)
      end <- as.integer(start + resol)
      list(chrom = chrom, start = start, end = end,
           score = ifelse(is.infinite(score) | is.nan(score), NA, score))
    }]
  }) %>%
    data.table::rbindlist()
  
  comps <- suppressWarnings(bedtorch::as.GenomicRanges(comps[!is.na(score)], genome = genome) %>% GenomicRanges::trim())
  
  if (!is_null(standard))
    comps %<>% normalize_compartment(standard = standard)
  
  return(comps)
}


#' Get compartment using FANC, from an existing .hic file
compartment_fanc <- function(hic_file,
                             norm = c("NONE", "VC", "VC_SQRT", "KR", "SCALE"),
                             resol,
                             ev = 1:2,
                             genome = c("hs37-1kg", "GRCh37", "GRCh38")) {
  assert_that(is_scalar_character(hic_file) && endsWith(hic_file, ".hic"))
  norm <- match.arg(norm)
  resol <- as.integer(resol)
  assert_that(is_scalar_integer(resol) && resol %in% c(250e3L, 500e3L, 1000e3L))
  assert_that(is_null(genome) || is_scalar_character(genome))
  
  ab_file <- tempfile(fileext = ".ab")
  on.exit(unlink(ab_file), add = TRUE)
  
  ev_tracks <- ev %>% map(function(ev) {
    ev_file <- tempfile(fileext = ".bed")
    on.exit(unlink(ev_file), add = TRUE)
    
    cmd <-
      str_interp("fanc compartments -v ${ev_file} -i ${ev} ${hic_file}@${resol%/%1000}kb@${norm} ${ab_file}")
    logging::loginfo(cmd)
    retcode <- system(cmd)
    assert_that(retcode == 0,
                msg = str_interp("Error in compartment calculation , RET: ${retcode}, CMD: ${cmd}"))
    
    
    gr <- bedtorch::read_bed(ev_file)
    # fanc output is 1-based bed
    GenomicRanges::start(gr) <- GenomicRanges::start(gr) - 1
    genome_info <- bedtorch::get_seqinfo(genome)
    seqlevels(gr, pruning.mode = "coarse") <- seqlevels(genome_info)
    seqinfo(gr) <- genome_info
    # GenomeInfoDb::keepSeqlevels(gr, value = bedtorch::get_seqinfo("GRCh38") %>% seqnames(), pruning.mode = "coarse")
    mcols(gr) <- data.frame(score = gr$V5)
    gr[gr$score != 0 & gr$score != 1]
  })
  
  # Combine eigenvectors
  comps <- purrr::reduce(seq_along(ev_tracks), .init = NULL, .f = function(x, idx) {
    y <- ev_tracks[[idx]]
    
    if (is_null(x)) {
      mcols(y) <- data.frame(ev_1 = y$score)
      return(y)
    }
    
    hits <- findOverlaps(x, y, type = "equal")
    x <- x[queryHits(hits)]
    mcols_x <- mcols(x)
    mcols_x[[paste0("ev_", idx)]] <- y[subjectHits(hits)]$score
    mcols(x) <- mcols_x
    
    return(x)
  })
  
  gc_track <- suppressWarnings({
    local({
      data_env <- env()
      data_name <- str_interp("gc.${genome}.${resol%/%1000}kbp")
      data(list = data_name, package = "hictools", envir = data_env)
      data_env[[data_name]]
    })
  })
  gc_track$score <- gc_track$gc
  gc_track$gc <- NULL
  
  if (!is_null(gc_track)) {
    comps <- normalize_compartment(comps, standard = gc_track, score_cols = colnames(mcols(comps)))

    score <- comps$score
    comps$score <- NULL
    mcols(comps) <- cbind(data.frame(score = score), mcols(comps))
  } else {
    score <- mcols(comps)[, 1]
    comps$score <- NULL
    mcols(comps) <- cbind(data.frame(score = score), mcols(comps))
  }
  
  return(comps)
}


#' Calculate O/E Pearson correlation matrix using Juicer tools
#'
#' @export
pearson_juicer <-
  function(hic_matrix,
           chrom,
           juicertools,
           java = "java",
           norm = "NONE") {
    stopifnot(length(chrom) == 1)
    resol <- attr(hic_matrix, "resol")
    pos_start <- min(c(hic_matrix$pos1, hic_matrix$pos2))
    
    temp_hic <- tempfile(fileext = ".hic")
    temp_matrix <- tempfile(fileext = ".txt")
    on.exit(file.remove(c(temp_hic, temp_matrix)))
    tryCatch({
      write_juicer_hic(
        hic_matrix = hic_matrix,
        file_path = temp_hic,
        juicertools = juicertools,
        java = java
      )
      cmd <- str_interp(
        paste0(
          "${java} -jar ${juicertools} pearsons ${norm} ",
          "${temp_hic} ${chrom} BP ${resol} ${temp_matrix}"
        )
      )
      cat(cmd, "\n")
      system(cmd)
      
      lines <- read_lines(temp_matrix) %>% str_trim()
      dim <- length(lines)
      mat <- lines %>%
        str_split(pattern = "[ ]+") %>%
        unlist() %>%
        as.numeric() %>%
        matrix(nrow = dim, byrow = TRUE)
      mat[is.nan(mat)] <- NA
      
      mat %>% convert_matrix_hic(chrom = chrom, resol = resol, pos_start = 0)
    }, finally = {
    })
  }


#' Calculate Pearson
#'
#' @export
pearson_ht <- function(hic_matrix, chrom, method = "lieberman") {
  stopifnot(length(chrom) == 1)
  pos_start <- min(c(hic_matrix$pos1, hic_matrix$pos2))
  resol <- attr(hic_matrix, "resol")
  oe <- oe_ht(hic_matrix, method = method)
  oe %>%
    convert_hic_matrix(chrom = chrom) %>%
    cor(use = "pairwise.complete.obs") %>%
    convert_matrix_hic(chrom = chrom, resol = resol, pos_start = pos_start)
}


#' Flip the compartment score signs according to `standard`
#'
#' @param compartment A `GRanges` or `data.frame` containing compartment scores.
#' @param standard A `GRanges` or `data.frame` of standard compartment scores.
#' @return The compartment. The signs of compartment scores may have been
#'   flipped based on `standard`.
#' @export
flip_compartment <- function(compartment, standard) {
  UseMethod("flip_compartment")
}

#' @export
flip_compartment.data.frame <- function(compartment, standard) {
  suppressWarnings(compartment %<>% bedtorch::as.GenomicRanges() %>% GenomicRanges::trim())
  compartment <- flip_compartment(compartment, standard)
  bedtorch::as.bedtorch_table(compartment)
}

#' @export
flip_compartment.GRanges <- function(compartment, standard) {
  if (inherits(standard, "data.frame"))
    suppressWarnings(standard %<>% bedtorch::as.GenomicRanges() %>% GenomicRanges::trim())
  
  compartment <- unique(GenomicRanges::seqnames(compartment)) %>%
    map(function(chrom) {
      compartment <- compartment[GenomicRanges::seqnames(compartment) == chrom]
      hits <- GenomicRanges::findOverlaps(compartment, standard)
      
      intervals_1 <- compartment[S4Vectors::queryHits(hits)]
      intervals_2 <- standard[S4Vectors::subjectHits(hits)]
      
      correlation <- cor(intervals_1$score, intervals_2$score, use = "complete.obs")
      if (correlation < 0)
        compartment$score <- compartment$score * (-1)
      
      compartment
    })
  
  rlang::exec(c, !!!compartment)
}


#' Normalize the compartment score signs according to `standard`
#'
#' @param compartment A `GRanges` or `data.frame` containing compartment scores.
#' @param standard A `GRanges` or `data.frame` of standard compartment scores.
#' @return The compartment. The signs of compartment scores may have been
#'   flipped based on `standard`.
#' @export
normalize_compartment <- function(compartment, standard, score_cols = c("score")) {
  UseMethod("normalize_compartment")
}

#' @export
normalize_compartment.data.frame <- function(compartment, standard, score_cols = c("score")) {
  suppressWarnings(compartment %<>% bedtorch::as.GenomicRanges() %>% GenomicRanges::trim())
  compartment <- normalize_compartment(compartment, standard, score_cols = score_cols)
  bedtorch::as.bedtorch_table(compartment)
}

#' @export
normalize_compartment.GRanges <- function(compartment, standard, score_cols = c("score")) {
  if (inherits(standard, "data.frame"))
    suppressWarnings(standard %<>% bedtorch::as.GenomicRanges() %>% GenomicRanges::trim())

  compartment <- unique(GenomicRanges::seqnames(compartment)) %>%
    map(function(chrom) {
      compartment <- compartment[GenomicRanges::seqnames(compartment) == chrom]
      
      # Calcualte the correlation with PC1, PC2, etc...
      cor_values <- map_dbl(score_cols, function(score_col) {
        compartment$score <- mcols(compartment)[[score_col]]
        hits <- GenomicRanges::findOverlaps(compartment, standard)
        intervals_1 <- compartment[S4Vectors::queryHits(hits)]
        intervals_2 <- standard[S4Vectors::subjectHits(hits)]
        cor(intervals_1$score, intervals_2$score, use = "complete.obs")
      })
      
      # Pick the one with the highest correlation coefficient, regardless of the sign
      # PC2 should be sufficiently larger than PC1 to be picked as compartment scores
      pick_index <- which.max(abs(cor_values))
      if (isTRUE(pick_index != 1) && abs(cor_values[pick_index]) - abs(cor_values[1]) > 0.5) {
        correlation_str <- paste(abs(cor_values), collapse = ", ")
        logging::logwarn(str_interp("Compartment for chromosome ${chrom}: choosing PC${pick_index}, correlation: ${correlation_str}"))
      } else
        pick_index <- 1
      
      if (cor_values[pick_index] < 0)
        compartment$score <- mcols(compartment)[[score_cols[pick_index]]] * (-1)
      else
        compartment$score <- mcols(compartment)[[score_cols[pick_index]]]
      
      return(compartment)
    })
  
  rlang::exec(c, !!!compartment)
}



#' Run hicPCA
#'
#' Run hicPCA and return the PCs, the O/E matrix and the Pearson correlation matrix
#'
#' @param hic_matrix A Hi-C matrix object
#' @export
hicexplorer_pca <- function(hic_matrix,
                            juicertools,
                            executable = "hicPCA",
                            executable_convert = "hicConvertFormat",
                            java = "java",
                            npc = 2,
                            chrom = NULL,
                            method = "lieberman") {
  if (is.null(chrom)) {
    # All chromosomes
    chrom <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2))
  }
  # Prepare temporary files
  hic_cool <- tempfile(fileext = ".cool")
  pc_bedgraph <-
    1:npc %>% map_chr(~ tempfile(fileext = ".bedGraph"))
  pm <- tempfile(fileext = ".cool")
  oem <- tempfile(fileext = ".cool")

  tryCatch({
    hic_matrix %>%
      write_cool(
        file_path = hic_cool,
        juicertools = juicertools,
        java = java,
        executable = executable_convert
      )

    pc_bedgraph_arg <- paste(pc_bedgraph, collapse = " ")
    cmd <- str_interp(
      paste0(
        "${executable} ",
        "-m ${hic_cool} ",
        "-noe ${npc} ",
        "-o ${pc_bedgraph_arg} ",
        "-f bedgraph --method ${method} ",
        "--chromosome ",
        paste(chrom, collapse = " "),
        " ",
        "-pm ${pm} ",
        "-oem ${oem} "
      )
    )
    cat(cmd, "\n")
    system(cmd)

    oe_matrix <- load_hic_cool(oem)
    pm_matrix <- load_hic_cool(pm)

    pc <-
      1:length(pc_bedgraph) %>%
      map_dfr(function(pc_idx) {
        bioessentials::load_genbed(pc_bedgraph[pc_idx]) %>%
          na.omit() %>%
          mutate(pc = paste0("PC", pc_idx)) %>%
          rename(score = X4)
      }) %>%
      pivot_wider(names_from = pc, values_from = score)
    list(pc = pc,
         oe = oe_matrix,
         pearson = pm_matrix)
  }, finally = {
    unlink(c(hic_cool, pc_bedgraph, pm, oem))
  })
}


# hictools's own implementation of compartment analysis, as opposed to juicer
compartment_ht <-
  function(hic_matrix,
           method = c("lieberman", "obs_exp", "nonzero", "average"),
           chrom = NULL,
           npc = 2L,
           standard = NULL) {
    assert_that(is(hic_matrix, "ht_table"))
    method <- match.arg(method)
    assert_that(is_scalar_integer(npc) && npc >= 1)
    
    assert_that(is_null(standard) ||
                  is(standard, "GRanges") || is(standard, "data.frame"))
    
    if (is_null(chrom))
      chrom <-
        unique(c(
          as.character(hic_matrix$chrom1),
          as.character(hic_matrix$chrom2)
        ))
    assert_that(is_character(chrom) && length(chrom) >= 1)
    
    resol <- attr(hic_matrix, "resol")
    assert_that(is_scalar_integer(resol) && resol > 0)
    type <- attr(hic_matrix, "type")
    assert_that(is_scalar_character(type))
    
    if (isTRUE(type == "observed"))
      hic_matrix <- oe_ht(hic_matrix, method = method)
    
    comps <- chrom %>% map(function(chrom) {
      hic_matrix <- hic_matrix[chrom1 == chrom & chrom2 == chrom]
      pos_start <- min(c(hic_matrix$pos1, hic_matrix$pos2))
      
      mc <- hic_matrix %>%
        convert_hic_matrix() %>%
        cor(use = "pairwise.complete.obs")
      
      # We need to exclude rows/columns which are all NAs
      all_na <-
        seq.int(ncol(mc)) %>% map_lgl(function(idx) {
          all(is.na(mc[, idx]))
        })
      
      # Sometimes, a single data point can be NA. This is different from the above case,
      # where the whole row/column is NA
      # We need to set these data points to 0 for PCA to work
      mc2 <- mc[!all_na,!all_na]
      mc2[which(is.na(mc2))] <- 0
      pc_subset <- prcomp(mc2, scale. = TRUE) %>% .$rotation %>% .[, 1:npc]
      pc <- matrix(rep(NA, ncol(mc) * npc), ncol = npc)
      pc[!all_na, ] <- pc_subset
      
      result <- data.table::data.table(
        chrom = chrom,
        start = as.integer((seq_along(pc[, 1]) - 1) * resol + pos_start)
      )
      result[, `:=`(end = as.integer(start + resol), score = pc[, 1])]
      result[, score := ifelse(is.infinite(score) | is.nan(score), NA, score)]
      data.table::setkey(result, "chrom", "start", "end")
      
      seq.int(npc) %>% walk(function(pc_idx) {
        result[, (paste0("PC", pc_idx)) := pc[, pc_idx]]
      })
      
      result
    }) %>%
      data.table::rbindlist() %>%
      bedtorch::as.bedtorch_table(attr(hic_matrix, "genome"))
    
    if (!is_null(standard))
      comps %<>% normalize_compartment(standard = standard, score_cols = paste0("PC", seq.int(npc)))
      # comps %<>% normalize_compartment(standard = standard)
      
    suppressWarnings(comps %>% bedtorch::as.GenomicRanges() %>% GenomicRanges::trim())
  }


#' Calculate correlation coefficients between two compartment tracks
#' 
#' @param x Compartment scores. Should be of bedtorch format.
#' @param y Compartment scores.
#' @return A named list of correlation coefficients
#' @export
comp_correlation <-
  function(x,
           y,
           correlation = c("pearson", "spearman", "kendall")) {
    x <- suppressWarnings(bedtorch::as.GenomicRanges(x) %>% GenomicRanges::trim())
    y <- suppressWarnings(bedtorch::as.GenomicRanges(y) %>% GenomicRanges::trim())
    
    correlation <- match.arg(correlation, several.ok = TRUE)
    
    hits <- GenomicRanges::findOverlaps(x, y)
    score_x <- x[S4Vectors::queryHits(hits)]$score
    score_y <- y[S4Vectors::subjectHits(hits)]$score
    
    result <- correlation %>% 
      map(~ cor(score_x, score_y, use = "complete.obs", method = .))
    names(result) <- correlation
    result
  }