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

    chroms %>% map_dfr(function(chrom) {
      strat_data <- hic_matrix %>%
        mutate(dist = abs(pos1 - pos2)) %>%
        group_by(dist) %>%
        summarize(
          avg = mean(score, na.rm = TRUE),
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
            # browser()
            row$total <-
              sum(expanded$total) / sum(expanded$nonzero) * row$nonzero
            row$avg <- row$total / row$nonzero
            row
          })
        strat_data <- bind_rows(strat_smoothed, strat_large)
      }
      strat_data %>% mutate(chrom = chrom) %>% select(chrom, everything())
    })
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

    chroms %>% map_dfr(function(chrom) {
      temp_hic <- tempfile(fileext = ".hic")
      hic_matrix %>%
        filter(chrom1 == chrom & chrom2 == chrom) %>%
        dump_juicer_hic(
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
           ref_genome = "hg19") {
    chrom_sizes %<>%
      filter(ref_genome == !!ref_genome) %>%
      mutate(n_bins = (size %/% resol) + 1)

    if (genome_wide) {
      # Genome-wide
      chrom <- NULL
    }

    if (is.null(chrom)) {
      chrom <- chrom_sizes %>%
        filter(ref_genome == !!ref_genome) %>%
        .$chrom
    }

    result <- chrom %>% map_dfr(function(chrom) {
      # if (!is.null(chrom)) {
      # Chromosome-wide
      chrom_sizes %<>% filter(chrom == !!chrom)
      # }

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
      result %>%
        group_by(distbin) %>%
        summarize(possible_dist = sum(possible_dist),
                  .groups = "drop")
    } else {
      result
    }
  }



# Calculate the observed/expected matrix
#
#' @export
oe_ht <- function(hic_matrix,
                  method,
                  genome_wide = FALSE,
                  smoothing = TRUE,
                  min_nonzero = 4) {
  if (method %in% c("lieberman", "obs_exp", "nonzero", "average")) {
    chroms <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2))
    resol <- attr(hic_matrix, "resol")

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
          mutate(observed = score, score = score / (total / (
            max_bin_idx - min_bin_idx + 1 - distbin
          )))
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
    })
  } else {
    stop(str_interp("Invalid method ${method}"))
  }
}


#' Calculate compartment scores using Juicer tools
#'
#' @export
compartment_juicer <-
  function(hic_matrix,
           juicertools,
           java = "java",
           norm = "NONE",
           resol = NULL) {
    if (is.null(resol)) {
      resol = attr(hic_matrix, "resol")
    }
    chroms <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2))

    # Dump the matrix to .hic
    hic_file <- tempfile(fileext = ".hic")
    hic_matrix %>% dump_hic(
      file_path = hic_file,
      format = "juicer_hic",
      juicertools = juicertools,
      java = java
    )

    comps <- chroms %>% map_dfr(function(chrom) {
      ev_file <- tempfile()
      cmd <-
        str_interp(
          "${java} -jar ${juicertools} eigenvector ${norm} ${hic_file} ${chrom} BP ${resol} ${ev_file}"
        )
      retcode <- system(cmd)
      if (retcode != 0) {
        unlink(hic_file)
        unlink(ev_file)
        stop(str_interp(
          "Error in compartment calculation , RET: ${retcode}, CMD: ${cmd}"
        ))
      }
      comps <-
        read_tsv(
          ev_file,
          col_names = "score",
          col_types = "n",
          na = c("", "NA", "NaN")
        ) %>%
        mutate(
          chrom = chrom,
          start = 0:(length(score) - 1) * resol,
          end = start + resol
        ) %>%
        select(chrom, start, end, score)
      unlink(ev_file)
      comps
    })

    unlink(hic_file)
    comps
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
    tryCatch({
      dump_juicer_hic(
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
      unlink(c(temp_hic, temp_matrix))
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


#' @export
compartment_ht <-
  function(hic_matrix,
           method = "lieberman",
           matrix = "observed",
           npc = 2,
           ...) {
    stopifnot(matrix %in% c("observed", "oe"))
    # stopifnot(method %in% c("juicer", "pca"))

    chroms <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2))
    resol <- attr(hic_matrix, "resol")

    if (matrix == "observed") {
      hic_matrix <- oe_ht(hic_matrix, method = method, ...)
    }

    chroms %>% map_dfr(function(chrom) {
      # hic_matrix %<>%
      #   filter(chrom1 == chrom & chrom2 == chrom) %>%
      #   mutate(x = pos1 %/% resol + 1,
      #          y = pos2 %/% resol + 1)
      # d <- max(c(hic_matrix$x, hic_matrix$y))
      # m <-
      #   matrix(
      #     data = 0,
      #     nrow = d,
      #     ncol = d,
      #     dimnames = list(1:d, 1:d)
      #   )
      # m[hic_matrix %>% select(x, y) %>% as.matrix()] <-
      #   hic_matrix$score
      # m[hic_matrix %>% select(y, x) %>% as.matrix()] <-
      #   hic_matrix$score
      # mc <- m %>% cor(use = "pairwise.complete.obs")

      hic_matrix %<>%
        filter(chrom1 == chrom & chrom2 == chrom)
      pos_start <- min(c(hic_matrix$pos1, hic_matrix$pos2))
      mc <- hic_matrix %>%
        convert_hic_matrix() %>%
        cor(use = "pairwise.complete.obs")

      all_na <-
        1:ncol(mc) %>% map_lgl(function(idx) {
          all(is.na(mc[, idx]))
        })

      pc_subset <-
        mc[!all_na,!all_na] %>% prcomp(scale. = TRUE) %>% .$rotation %>% .[, 1:npc]
      pc <- matrix(rep(NA, ncol(mc) * npc), ncol = npc)
      pc[!all_na, ] <- pc_subset

      result <- tibble(
        chrom = chrom,
        start = as.integer((seq_along(pc[, 1]) - 1) * resol + pos_start),
        end = start + resol,
        score = pc[, 1]
      )
      #
      1:npc %>% walk(function(pc_idx) {
        result[[paste0("PC", pc_idx)]] <<- pc[, pc_idx]
      })
      result %>%
        mutate(score = PC1) %>%
        select(chrom, start, end, score, everything())
    })
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
      dump_cool(
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
