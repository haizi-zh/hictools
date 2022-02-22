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
          norm = norm
        )
      oe_matrix <-
        load_juicer_hic(
          file_path = temp_hic,
          chrom = chrom,
          resol = resol,
          type = "oe",
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
    ht_table(
      resol = resol,
      type = "oe",
      norm = norm,
      genome = attr(hic_matrix, "genome")
    )
}


#' Calculate compartment scores using Juicer tools
#'
#' @param hic_matrix Either a `ht_table` object, or the path of the `.hic` name.
#'   If a file path, `resol` cannot be NULL since we can't infer resolution
#'   directly from the `.hic` file.
#' @param chrom A character vector indicating which chromosomes to run. If
#'   `NULL`, will calculate compartment scores for all chromosomes. However, in
#'   this case, `hic_matrx` cannot be the `.hic` file name.
#' @param resol Resolution of the matrix. Must be a positive integer when the
#'   input is a Hi-C file.
#' @param reference A reference compartment track which is used to "orient"
#'   principal components. If `NULL`, no orienting will be performed. Can be the
#'   name of a pre-computed reference track, e.g. gc.GRCh38.500kbp (see
#'   `data(package = "hictools")`), or path to the bed file for the reference
#'   track, or a `GRanges` object.
#' @param oe The method to compute observed/expected values. `"juicer"`
#'   indicates using using the values computed by Juicer tools. `"ht"` indicates
#'   using hictools for the calculation. `"as-is"` indicates observed values
#'   will be used as observed/expected ones, thus skipping the OE calculation.
#'   This may be useful when the input is cofrag matrix rather than standard
#'   Hi-C data. Won't take effect if `method` is `"juicer"`.
#' @param norm The norm used for calculating compartment scores. If `hic_matrix`
#'   is not a `.hic` file name, `norm` will be ignored.
#' @export
get_compartment <- function(hic_matrix,
                            method = c("juicer", "lieberman", "obs_exp", "nonzero", "average", "fanc"),
                            chrom = NULL,
                            reference = NULL,
                            oe = c("juicer", "ht", "as-is"),
                            norm = c("NONE", "VC", "VC_SQRT", "KR", "SCALE"),
                            ...) {
  UseMethod("get_compartment")
}


# Convert genome to hg19 or hg38
get_juicer_genome <- function(genome) {
  if (genome %in% c("hs37-1kg", "hs37d5", "hg19") || grepl(pattern = "^GRCh37(\\.p[0-9]+)?$", x = genome))
    juicer_genome <- "hg19"
  else if (genome == "hg38" || grepl(pattern = "^GRCh38(\\.p[0-9]+)?$", x = genome))
    juicer_genome <- "hg38"
  else
    stop(pastep("Unknown genome: ", genome))
  return(juicer_genome)
}


#' Calculate compartment scores
#' 
#' Calculate compartment scores for Hi-C matrix in ht_table format.
#' 
#' @param hic_matrix A `ht_table` object for the Hi-C matrix.
#' 
#' @export
get_compartment.ht_table <- function(hic_matrix,
                                     method = c("juicer", "lieberman", "obs_exp", "nonzero", "average", "fanc"),
                                     chrom = NULL,
                                     reference = NULL,
                                     oe = c("juicer", "ht", "as-is"),
                                     juicertools = get_juicer_tools(),
                                     java = "java",
                                     ...) {
  assert_that(
    hic_type(hic_matrix) %in% c("observed", "cofrag") || (hic_type(hic_matrix) == "oe" && oe == "as-is"),
    msg = paste0("Incorrect HiC data type: ", hic_type(hic_matrix))
  )
  method <- match.arg(method)
  oe <- match.arg(oe)
    
  if (is_null(chrom))
    chrom <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2))
  else
    chrom <- as.character(chrom)
  assert_that(is_character(chrom) && length(chrom) >= 1)
  
  resol <- hic_resol(hic_matrix)
  assert_that(is_valid_resol(resol) && resol >= 100e3L)
  
  genome <- hic_genome(hic_matrix)
  
  if (method %in% c("juicer", "fanc") || oe == "juicer") {
    # Need to convert the input to .hic files if we need to calculate Juicer
    # eigenvectors, or use Juicer tools to get oe scores.
    
    # Save hic to a temporary .hic file, and call juicer tools to get compartments
    hic_file <- tempfile(fileext = ".hic")
    on.exit(file.remove(hic_file), add = TRUE)
    
    hic_matrix %>% write_hic(
      file_path = hic_file,
      format = "juicer_hic",
      juicertools = juicertools,
      java = java
    )
    
    if (method == "juicer")
      comps <- get_compartment.character(
        hic_matrix = hic_file,
        method = method,
        chrom = chrom,
        reference = reference,
        resol = resol,
        genome = genome,
        juicertools = juicertools,
        java = java,
        norm = "NONE",
        ...
      )
    else if (method == "fanc")
      comps <- get_compartment.character(
        hic_matrix = hic_file,
        method = method,
        chrom = chrom,
        reference = reference,
        resol = resol,
        genome = genome,
        norm = "NONE",
        ...
      )
    else {
      comps <- chrom %>%
        map(function(chrom) {
            load_juicer_hic(
              file_path = hic_file,
              chrom = chrom,
              resol = resol,
              type = "oe",
              norm = "NONE",
              genome = genome
            )
        }) %>% 
        concat_hic() %>%
        get_compartment.ht_table(
          method = method,
          chrom = chrom,
          reference = reference,
          oe = "as-is",
          juicertools = juicertools,
          java = java,
          ...
        )
    }
  } else {
    comps <- compartment_ht(
      hic_matrix = hic_matrix,
      method = method,
      chrom = chrom,
      npc = 2L,
      oe = oe
    )
  }
  
  if (!is_null(reference)) {
    # PC1/2/3, ..., or simply score
    score_cols <- colnames(mcols(comps)) %>% grep(pattern = "^PC[0-9]+", value = TRUE)
    if (length(score_cols) == 0)
      score_cols <- colnames(mcols(comps))[1]
    
    comps <- orient_compartment(comps, reference, score_cols = score_cols)
  }

  return(comps)
}


#' Calculate compartment scores
#' 
#' Load Hi-C matrix from .hic file and calculate compartment scores.
#' 
#' @param genome Name of the reference genome.
#' @param resol A positive integer for resolution.
#' @export
get_compartment.character <- function(hic_matrix,
                                      method = c("juicer", "lieberman", "obs_exp", "nonzero", "average", "fanc"),
                                      chrom = NULL,
                                      resol,
                                      reference = NULL,
                                      oe = c("juicer", "ht", "as-is"),
                                      genome = NULL,
                                      norm = c("NONE", "VC", "VC_SQRT", "KR", "SCALE"),
                                      ...) {
  assert_that(grepl("\\.hic$", x = hic_matrix) && assertthat::is.readable(hic_matrix))
  method <- match.arg(method)
  oe <- match.arg(oe)
  norm <- match.arg(norm)
  
  # Must provide chromosome names unless work with FAN-C
  if (method != "fanc")
    assert_that(is_character(chrom) && length(chrom) >= 1)
  
  assert_that(is_valid_resol(resol) && resol >= 100e3L)
  assert_that(!is_null(bedtorch::get_seqinfo(genome = genome)))
  
  if (method == "juicer") {
    # No need to load the .hic file. Instead, directly call compartments using Juicer tools
    comps <- compartment_juicer_file(
      hic_file = hic_matrix,
      norm = norm,
      chrom = chrom,
      resol = resol,
      genome = genome,
      ...
    )
  } else if (method == "fanc") {
    # Check fanc availability
    fanc_version <-
      tryCatch({
        system2(
          command = "fanc",
          args = "--version",
          stdout = TRUE,
          stderr = FALSE
        )
      }, error = function(e)
        NULL)
    assert_that(!is.null(fanc_version), msg = "FAN-C is not available")
    # stop("FAN-C compartment is not supported at this point")
    # No need to load the .hic file. Instead, directly call compartments using FAN-C
    comps <- compartment_fanc(
      hic_file = hic_matrix,
      norm = norm,
      resol = resol,
      genome = genome,
      ...
    )
  } else {
    if (oe == "juicer") {
      juicer_type <- "oe"
      oe <- "as-is"
    } else {
      juicer_type <- "observed"
    }
    
    comps <- chrom %>%
      map(function(chrom) {
        load_juicer_hic(
          file_path = hic_matrix,
          chrom = chrom,
          resol = resol,
          # Which type to load?
          type = juicer_type, 
          norm = "NONE",
          genome = genome
        )
      }) %>% 
      concat_hic() %>%
      get_compartment.ht_table(
        method = method,
        chrom = chrom,
        reference = reference,
        oe = oe,
        ...
      )
    
    
    # 
    # 
    # hic_matrix <-
    #   load_juicer_hic(
    #     file_path = hic_matrix,
    #     resol = resol,
    #     chrom = chrom,
    #     # Which type to load?
    #     type = switch(oe, juicer = "oe", "observed"),
    #     norm = norm, 
    #     genome = genome
    #   )
    # comps <- get_compartment(
    #   hic_matrix = hic_matrix,
    #   method = method,
    #   chrom = chrom,
    #   reference = reference,
    #   oe = switch(oe, juicer = "as-is", oe)
    # )
  }
  
  
  if (!is_null(reference)) {
    # PC1/2/3, ..., or simply score
    score_cols <- colnames(mcols(comps)) %>% grep(pattern = "^PC[0-9]+", value = TRUE)
    if (length(score_cols) == 0)
      score_cols <- colnames(mcols(comps))[1]
    
    comps <- orient_compartment(comps, reference, score_cols = score_cols)
  }
  
  return(comps)
}


#' Compute compartment scores using Juicer tools, from an existing .hic file
#' 
#' @param bpparam BiocParallelParam for parallel operation
compartment_juicer_file <- function(hic_file,
                                    juicertools = get_juicer_tools(),
                                    java = "java",
                                    norm = c("NONE", "VC", "VC_SQRT", "KR", "SCALE"),
                                    chrom = NULL,
                                    resol = NULL,
                                    genome = NULL,
                                    bpparam = BiocParallel::bpparam("SerialParam")) {
  assert_that(grepl(pattern = "\\.hic$", x = hic_file) && assertthat::is.readable(hic_file))
  assert_that(is_scalar_character(java))
  assert_that(is_scalar_character(juicertools))
  norm <- match.arg(norm)
  chrom <- as.character(chrom)
  assert_that(length(chrom) >= 1)
  resol <- as.integer(resol)
  assert_that(is_valid_resol(resol) && resol >= 100e3L)
  assert_that(is_null(genome) || (is_scalar_character(genome) && !is_null(bedtorch::get_seqinfo(genome))))
  
  comps <- BiocParallel::bplapply(chrom, BPPARAM = bpparam, FUN = function(chrom) {
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
  
  # Need to trim to fit the genome
  comps <- suppressWarnings(bedtorch::as.GenomicRanges(comps[!is.na(score)], genome = genome) %>% GenomicRanges::trim())
  
  return(comps)
}


#' Compute compartment scores using hictools's own implementation of compartment analysis, as opposed to juicer
#' 
#' @param npc A positive integer indicating how many principal components to keep.
compartment_ht <-
  function(hic_matrix,
           method = c("lieberman", "obs_exp", "nonzero", "average"),
           oe = c("ht", "as-is"),
           chrom = NULL,
           npc = 2L) {
    assert_that(is(hic_matrix, "ht_table"))
    method <- match.arg(method)
    oe <- match.arg(oe)
    assert_that(is_scalar_integer(npc) && npc >= 1)
    
    if (is_null(chrom))
      chrom <-
      unique(c(
        as.character(hic_matrix$chrom1),
        as.character(hic_matrix$chrom2)
      ))
    assert_that(is_character(chrom) && length(chrom) >= 1)
    
    resol <- attr(hic_matrix, "resol")
    assert_that(is_valid_resol(resol) && resol >= 100e3L)
    type <- attr(hic_matrix, "type")
    assert_that(is_scalar_character(type))
    
    genome <- attr(hic_matrix, "genome")
    
    if (oe != "as-is")
      hic_matrix <- oe_ht(hic_matrix, method = method)
    
    comps <- chrom %>% map(function(chrom) {
      hic_matrix <- hic_matrix[chrom1 == chrom & chrom2 == chrom]
      pos_start <- min(c(hic_matrix$pos1, hic_matrix$pos2))
      
      mc <- hic_matrix %>%
        convert_hic_matrix() %>%
        cor(use = "pairwise.complete.obs")
      
      # We need to exclude rows/columns which are all NAs
      all_na <- apply(mc, MARGIN = 1, FUN = function(x) all(is.na(x)))
      
      # Sometimes, a single data point can be NA. This is different from the above case,
      # where the whole row/column is NA
      # We need to set these data points to 0 for PCA to work
      mc2 <- mc[!all_na, !all_na]
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
      data.table::rbindlist()
    
    return(suppressWarnings(comps %>% 
                                bedtorch::as.GenomicRanges(genome = genome) %>% 
                                GenomicRanges::trim()))
    
    # gc_track <- suppressWarnings({
    #   local({
    #     data_env <- env()
    #     data_name <- str_interp("gc.${genome}.${resol%/%1000}kbp")
    #     data(list = data_name, package = "hictools", envir = data_env)
    #     data_env[[data_name]]
    #   })
    # })
    # gc_track$score <- gc_track$gc
    # gc_track$gc <- NULL
    # 
    # comps <- suppressWarnings(comps %>% bedtorch::as.GenomicRanges() %>% GenomicRanges::trim())
    # 
    # if (!is_null(gc_track)) {
    #   comps <- orient_compartment(comps, reference = gc_track, score_cols = colnames(mcols(comps)))
    #   
    #   score <- comps$score
    #   comps$score <- NULL
    #   mcols(comps) <- cbind(data.frame(score = score), mcols(comps))
    # } else {
    #   score <- mcols(comps)[, 1]
    #   comps$score <- NULL
    #   mcols(comps) <- cbind(data.frame(score = score), mcols(comps))
    # }
    # 
    # 
    # if (!is_null(reference))
    #   comps %<>% orient_compartment(reference = reference, score_cols = paste0("PC", seq.int(npc)))
    # # comps %<>% orient_compartment(reference = reference)
    
    # suppressWarnings(comps %>% bedtorch::as.GenomicRanges() %>% GenomicRanges::trim())
  }


#' Get compartment using FANC, from an existing .hic file
#' 
#' @param ev an integer vector specifying which eigenvectors to calculate
#' @param ab_file path to a pre-calculated AB matrix file
#' @param keep_ab if `TRUE`, keep the AB matrix file. This requires `ab_file =
#'   NULL`
compartment_fanc <- function(hic_file,
                             norm = c("NONE", "VC", "VC_SQRT", "KR", "SCALE"),
                             resol,
                             ev = 1:2,
                             genome = NULL,
                             ab_file = NULL,
                             keep_ab = !is_null(ab_file)) {
  assert_that(is_scalar_character(hic_file) && endsWith(hic_file, ".hic"))
  norm <- match.arg(norm)
  
  resol <- as.integer(resol)
  assert_that(is_scalar_integer(resol) &&
                resol %in% strawr::readHicBpResolutions(hic_file),
              msg = str_interp("Resolution ${resol} does not exist in ${hic_file}"))
  assert_that(is_valid_resol(resol) && resol >= 100e3L)

  # Test genome validity
  assert_that(is_null(genome) || !is_null(bedtorch::get_seqinfo(genome)))
  
  assert_that(is_scalar_logical(keep_ab))
  
  if (is_null(ab_file)) {
    ab_file <- tempfile(fileext = ".ab")
    if (!keep_ab) {
      on.exit(unlink(ab_file), add = TRUE)
    }
  } else {
    assert_that(is_scalar_character(ab_file))
  }
  
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
    # ev_file may contains seqnames not available in genome
    # trim them
    genome_info <- bedtorch::get_seqinfo(genome)
    seqlevels(gr, pruning.mode = "coarse") <- seqlevels(genome_info)
    seqinfo(gr) <- genome_info

    # fanc output is 1-based bed
    GenomicRanges::start(gr) <- GenomicRanges::start(gr) - 1
    mcols(gr) <- data.frame(score = gr$V5)
    gr[!is.na(gr$score) & gr$score != 0 & gr$score != 1]
  })
  
  # Combine eigenvectors
  comps <- purrr::reduce(seq_along(ev_tracks), .init = NULL, .f = function(x, idx) {
    y <- ev_tracks[[idx]]
    
    if (is_null(x)) {
      mcols(y) <- data.frame(PC1 = y$score)
      return(y)
    }
    
    hits <- findOverlaps(x, y, type = "equal")
    x <- x[queryHits(hits)]
    mcols_x <- mcols(x)
    mcols_x[[paste0("PC", idx)]] <- y[subjectHits(hits)]$score
    mcols(x) <- mcols_x
    
    return(x)
  })
  
  return(comps)
}


#' Calculate O/E Pearson correlation matrix using Juicer tools
#'
#' @export
pearson_juicer <-
  function(hic_matrix,
           chrom = NULL,
           juicertools = get_juicer_tools(),
           java = "java",
           norm = "NONE") {
    assert_that(is(hic_matrix, "ht_table"))
    
    all_chroms <- hic_matrix[, unique(chrom1)]
    if (is.null(chrom))
      chrom <- all_chroms
    else {
      assert_that(all(chrom %in% all_chroms))
    }
    
    resol <- hic_resol(hic_matrix)
    pos_start <- min(c(hic_matrix$pos1, hic_matrix$pos2))
    
    result <- map(chrom, function(chrom) {
      hic_matrix <- hic_matrix[chrom1 == chrom & chrom1 == chrom2]
      
      temp_hic <- tempfile(fileext = ".hic")
      temp_matrix <- tempfile(fileext = ".txt")
      on.exit(file.remove(c(temp_hic, temp_matrix)))
      tryCatch({
        write_juicer_hic(
          hic_matrix = hic_matrix,
          file_path = temp_hic,
          juicertools = juicertools,
          java = java,
          norm = "NONE"
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
        
        hic_pearson <- mat %>% convert_matrix_hic(
          chrom = chrom,
          resol = resol,
          pos_start = 0,
          copy_from = hic_matrix
        )
        hic_type(hic_pearson) <- "pearson"
        return(hic_pearson)
      }, finally = {
        
      })
    })
    
    data.table::rbindlist(result) %>% ht_table(copy_from = result[[1]])
  }


#' Calculate Pearson
#'
#' @param method `lieberman` indicates to calculate Pearson correlation matrix
#'   based on observed/expected values, using the method described in Lieberman
#'   2009. `none` indicates directly calculate Pearson correlation from the Hi-C
#'   matrix.
#' @export
pearson_ht <- function(hic_matrix, chrom, method = c("lieberman", "none")) {
  method = match.arg(method)
  stopifnot(length(chrom) == 1)
  pos_start <- min(c(hic_matrix$pos1, hic_matrix$pos2))
  resol <- attr(hic_matrix, "resol")
  if (method == "none")
    oe <- hic_matrix
  else
    oe <- oe_ht(hic_matrix, method = method)
  hic_pearson <- oe %>%
    convert_hic_matrix(chrom = chrom) %>%
    cor(use = "pairwise.complete.obs") %>%
    convert_matrix_hic(
      chrom = chrom,
      resol = resol,
      pos_start = pos_start,
      copy_from = hic_matrix
    )
  hic_type(hic_pearson) <- "pearson"
  return(hic_pearson)
}


#' Orient the compartment score signs according to `reference`
#'
#' @param compartment A `GRanges` object for compartment scores.
#' @param reference A `GRanges` for the reference compartment track.
#' @param score_cols A character vector indicating which columns are candidates
#'   for compartment scores. For example, `c("PC1", "PC2")` means from `PC1` and
#'   `PC2`, picking the one most likely to be associated with A/B compartments,
#'   and orient it if necessary.
#' @return The compartment. The signs of compartment scores may have been
#'   flipped based on `reference`.
#' @export
orient_compartment <- function(compartment, reference, score_cols = c("score")) {
  assert_that(is_character(score_cols) && length(score_cols) >= 1)
  assert_that(is(compartment, "GRanges"))
  assert_that(is(reference, "GRanges"))
  
  # Check for compatibility
  resol1 <- local({
    dt <- compartment %>% bedtorch::as.bedtorch_table()
    dt[, dplyr::lead(start) - start, by = "chrom"][, min(V1, na.rm = TRUE)]
  })
  resol2 <- local({
    dt <- reference %>% bedtorch::as.bedtorch_table()
    dt[, dplyr::lead(start) - start, by = "chrom"][, min(V1, na.rm = TRUE)]
  })
  assert_that(is_valid_resol(resol1) && is_valid_resol(resol2) && resol1 == resol2)
  
  # Make compatible genomes
  original_seqinfo <- seqinfo(compartment)
  seqlevels(compartment) <- seqlevels(reference)
  seqinfo(compartment) <- seqinfo(reference)
  
  # By default, the first column of reference is used for determining the orientation
  reference$score <- mcols(reference)[[1]]
  
  compartment <- unique(seqnames(compartment)) %>%
    map(function(chrom) {
      compartment <- compartment[seqnames(compartment) == chrom]
      
      hits <- findOverlaps(compartment, reference)
      if (length(hits) == 0) {
        warning(paste0("Chromosome ", chrom, " does not exist in the reference track"),
                call. = FALSE)
        return(compartment)
      }
      
      # Calcualte the correlation with PC1, PC2, etc...
      for (score_name in score_cols) {
        original_score <- mcols(compartment)[[score_name]]
        cor_value <- cor(original_score[queryHits(hits)],
                         reference[subjectHits(hits)]$score, 
                         use = "complete.obs")
        mcols(compartment)[[score_name]] <-
          sign(cor_value) * mcols(compartment)[[score_name]]
      }
      
      return(compartment)
    })
  
  compartment <- rlang::exec(c, !!!compartment)
  seqlevels(compartment) <- seqlevels(original_seqinfo)
  seqinfo(compartment) <- original_seqinfo
  return(compartment)
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


#' Calculate correlation coefficients between two compartment tracks
#'
#' @param x Compartment scores.
#' @param y Compartment scores.
#' @param score_x name of the score column in `x`
#' @param score_y name of the score column in `y`
#' @param method Type of correlation scores. Must be one of `"pearson"`,
#'   `"spearman"`, and `"kendall"`.
#' @param overall if FALSE, calculate correlation scores for each chromosome
#'   separately
#' @return A data frame for correlation.
#' @export
comp_correlation <-
  function(x,
           y,
           score_x = 1,
           score_y = 1,
           method = c("pearson", "spearman", "kendall"),
           overall = FALSE) {
    assert_that(is(x, "GRanges"))
    assert_that(is(y, "GRanges"))
    assert_that(is_scalar_logical(overall))
    method <- match.arg(method, several.ok = TRUE)
    
    chroms <- intersect(seqnames(x), seqnames(y))
    seqlevels(x, pruning.mode = "coarse") <- chroms
    seqlevels(y, pruning.mode = "coarse") <- chroms
    seqinfo(y) <- seqinfo(x)
    
    
    cor_helper <- function(x, y, method) {
      hits <- findOverlaps(x, y)
      score_x <- mcols(x[queryHits(hits)])[[score_x]]
      score_y <- mcols(y[subjectHits(hits)])[[score_y]]
      return(cor(score_x, score_y, use = "complete.obs", method = method))
    }
    
    if (overall) {
      map_dfr(method, function(method) {
        tibble(method = method, cor = cor_helper(x, y, method))
      })
    } else {
      expand_grid(method = method,
                  chrom = chroms) %>%
        pmap_dfr(function(method, chrom) {
          x <- x[seqnames(x) == chrom]
          y <- y[seqnames(y) == chrom]
          
          tibble(chrom = chrom,
                 method = method,
                 cor = cor_helper(x, y, method))
        })
    }
  }


#' Normalize Hi-C matrix
#' 
#' @param norm Currently only KR is supported
#' @export
normalize_hic <- function(hic_matrix,
                          norm = "KR") {
  assert_that(is(hic_matrix, "ht_table"))
  norm <- match.arg(norm)
  
  chrom <- hic_matrix[, unique(c(chrom1, chrom2))]
  
  result <- map(chrom, function(chrom) {
    hic_matrix <- hic_matrix[chrom1 == chrom & chrom2 == chrom]
    pos_start <- hic_matrix[, min(c(pos1, pos2))]
    
    convert_hic_matrix(hic_matrix = hic_matrix,
                                 chrom = chrom,
                                 missing_score = 0) %>%
      HiCcompare::KRnorm() %>%
      convert_matrix_hic(
        chrom = chrom,
        resol = hic_resol(hic_matrix),
        pos_start = pos_start,
        copy_from = hic_matrix
      )
  }) %>%
    data.table::rbindlist() %>%
    ht_table(copy_from = hic_matrix)
  
  hic_norm(result) <- "KR"
  return(result)
}