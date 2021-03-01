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
# The above copyright notice and this permission notice shall be included in all
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
# Tools for analyzing and manipulating Hi-C dataset


supported_resol <- c(2.5e6L,
                     1e6L,
                     500e3L,
                     250e3L,
                     100e3L,
                     50e3L,
                     25e3L,
                     10e3L,
                     5e3L,
                     2.5e3L,
                     1e3L)


#' Load Hi-C dataset from a .hic file
#'
#' This function invokes Juicer tools for the dumping
#' @param file_path Path to the .hic file
#' @param chrom A character vector to indicate what chromosomes to load
#' @param resol Resolution in base pair (BP)
#' @param matrix Should be \code{observed} or \code{oe}
#' @param norm See the manual of Juicer tools.
#' @return A Hi-C object
#' @export
load_juicer_hic <- function(file_path, chrom, resol, matrix = "observed", norm = "NONE") {
  stopifnot(!is.null(resol) && resol %in% supported_resol)
  stopifnot(matrix %in% c("observed", "oe"))
  stopifnot(norm %in% c("NONE", "KR", "VC", "VC_SQRT"))


  chrom %>% map_dfr(function(chrom) {
    tryCatch({
      strawr::straw(
        matrix = matrix,
        norm = norm,
        fname = file_path,
        chr1loc = chrom,
        chr2loc = chrom,
        unit = "BP",
        binsize = resol
      ) %>%
        as_tibble() %>%
        mutate(
          chrom1 = chrom,
          chrom2 = chrom,
          pos1 = x,
          pos2 = y,
          score = counts
        ) %>%
        select(chrom1, pos1, chrom2, pos2, score)
    }, error = function(e) {
      warning(str_interp("File doesn't have data for chromosome: ${chrom}"))
      NULL
    })
  }) %>%
    set_attr("resol", as.integer(resol)) %>%
    set_attr("type", matrix) %>%
    set_attr("norm", norm) %>%
    set_attr("chrom", sort(chrom)) %>%
    arrange("chrom1", "chrom2", "pos1", "pos2")
}


#' @export
load_juicer_short <- function(file_path, chrom = NULL, matrix = "unknown", norm = "unknown") {
  # 0 22 16000000 0 0 22 16000000 1 95
  data <- read_delim(
    file = file_path,
    delim = " ",
    col_names = F,
    col_types = "iciiiciin"
  ) %>%
    rename(
      chrom1 = X2,
      pos1 = X3,
      chrom2 = X6,
      pos2 = X7,
      score = X9
    ) %>%
    select(chrom1, pos1, chrom2, pos2, score)

  if (!is.null(chrom)) {
    data %<>% filter(chrom1 %in% chrom & chrom2 %in% chrom)
  } else {
    chrom <- c(data$chrom1, data$chrom2) %>% unique()
  }

  attr(data, "resol") <- guess_resol(data)
  data %>%
    set_attr("type", matrix) %>%
    set_attr("norm", norm) %>%
    set_attr("chrom", sort(chrom)) %>%
    arrange("chrom1", "chrom2", "pos1", "pos2")
}

#' @export
load_juicer_dump <- function(file_path, chrom, matrix = "unknown", norm = "unknown") {
  stopifnot(length(chrom) == 1)
  # 16000000        16000000        95.0
  data <- read_tsv(
    file = file_path,
    col_names = c("pos1", "pos2", "score"),
    col_types = "iin"
  ) %>%
    mutate(chrom1 = chrom, chrom2 = chrom) %>%
    select(chrom1, pos1, chrom2, pos2, score)
  attr(data, "resol") <- guess_resol(data)
  data %>%
    set_attr("type", matrix) %>%
    set_attr("norm", norm) %>%
    set_attr("chrom", sort(chrom)) %>%
    arrange("chrom1", "chrom2", "pos1", "pos2")
}


#' @export
load_hic_genbed <- function(file_path, chrom = NULL, matrix = "unknown", norm = "unknown") {
  data <- bioessentials::load_genbed(file_path = file_path) %>%
    modify_at(4, as.character) %>%
    modify_at(5:6, as.integer) %>%
    modify_at(7, as.numeric) %>%
    select(-3, -6)
  colnames(data)[1:5] <- c("chrom1", "pos1", "chrom2", "pos2", "score")

  if (!is.null(chrom)) {
    data %<>% filter(chrom1 %in% chrom & chrom2 %in% chrom)
  } else {
    chrom <- c(data$chrom1, data$chrom2) %>% unique()
  }

  attr(data, "resol") <- guess_resol(data)
  data %>%
    set_attr("type", matrix) %>%
    set_attr("norm", norm) %>%
    set_attr("chrom", sort(chrom)) %>%
    arrange("chrom1", "chrom2", "pos1", "pos2")
}


#' Load Hi-C data in cool format
#' @export
load_hic_cool <- function(file_path, chrom = NULL, matrix = "unknown", norm = "unknown", hdf5 = TRUE, cooler = "cooler") {
  stopifnot(hdf5)

  # Load the bin table
  if (hdf5) {
    bin_table <- h5read(file_path, "bins") %>%
      as_tibble() %>%
      mutate(chrom = as.character(chrom), bin_id = row_number() - 1)
  } else {
    cmd <-
      str_interp("cooler dump ${file_path} -t bins -H -c chrom,start,end")
    bin_table <-
      read.delim2(file = textConnection(system(cmd, intern = TRUE)), sep = "\t") %>%
      as_tibble() %>%
      mutate(bin_id = row_number() - 1)
  }

  helper_cool <- function(chrom) {
    cmd <-
      str_interp("cooler dump ${file_path} -H")
    if (!is.null(chrom)) {
      cmd <- paste0(cmd, " -r ", chrom)
    }
    observed <-
      read.delim2(file = textConnection(system(cmd, intern = TRUE)),
                  sep = "\t") %>%
      as_tibble() %>%
      modify_at(1:2, as.integer) %>%
      modify_at(3, as.numeric) %>%
      inner_join(x = .,
                 y = bin_table,
                 by = c(bin1_id = "bin_id")) %>%
      inner_join(x = .,
                 y = bin_table,
                 by = c(bin2_id = "bin_id")) %>%
      rename(
        chrom1 = chrom.x,
        chrom2 = chrom.y,
        pos1 = start.x,
        pos2 = start.y,
        score = count
      ) %>%
      select(chrom1, pos1, chrom2, pos2, score)
  }

  helper_h5 <- function() {
    h5read(file_path, "pixels") %>% as_tibble() %>%
      inner_join(x = .,
                 y = bin_table,
                 by = c(bin1_id = "bin_id")) %>%
      inner_join(x = .,
                 y = bin_table,
                 by = c(bin2_id = "bin_id")) %>%
      rename(
        chrom1 = chrom.x,
        chrom2 = chrom.y,
        pos1 = start.x,
        pos2 = start.y,
        score = count
      ) %>%
      select(chrom1, pos1, chrom2, pos2, score)
  }

  # Load the contents
  observed <- if (is.null(chrom)) {
    # Genome-wide dump
    if (hdf5) {
      helper_h5()
    } else {
      helper_cool(chrom = chrom)
    }
  } else {
    if (hdf5) {
      helper_h5() %>% filter(chrom1 %in% chrom & chrom2 %in% chrom)
    } else {
      chrom %>% map_dfr(helper_cool)
    }
  }

  if (!is.null(chrom)) {
    observed %<>% filter(chrom1 %in% chrom & chrom2 %in% chrom)
  } else {
    chrom <- c(observed$chrom1, observed$chrom2) %>% unique()
  }

  attr(observed, "resol") <- guess_resol(observed)
  observed %>%
    set_attr("type", matrix) %>%
    set_attr("norm", norm) %>%
    set_attr("chrom", sort(chrom)) %>%
    arrange("chrom1", "chrom2", "pos1", "pos2")
}


# Guess Hi-C format from file name
guess_format <- function(file_path) {
  if (endsWith(file_path, ".hic")) {
    return("juicer_hic")
  } else if (endsWith(file_path, ".bed") || endsWith(file_path, ".bed.gz")) {
    return("genbed")
  } else if (endsWith(file_path, ".short") || endsWith(file_path, ".short.gz")) {
    return("juicer_short")
  } else if(endsWith(file_path, ".cool")) {
    return("cool")
  } else {
    stop(str_interp("Unknown input format: ${file_path}"))
  }
}


# Guess the resolution from data
guess_resol <- function(data) {
  sorted_pos <- c(data$pos1, data$pos2) %>% unique() %>% sort()
  sorted_pos[2] - sorted_pos[1]
}


#' Load Hi-C dataset from file
#'
#' @export
load_hic <- function(file_path, format = NULL, resol = NULL, ...) {
  if (is.null(format)) {
    format <- guess_format(file_path)
  }
  if (format == "juicer_short") {
    data <- load_juicer_short(file_path, ...)
  } else if (format == "juicer_dump") {
    data <- load_juicer_dump(file_path, ...)
  } else if (format == "juicer_hic") {
    data <- load_juicer_hic(file_path, resol = resol, ...)
  } else if (format == "genbed") {
    data <- load_hic_genbed(file_path = file_path)
  } else if (format == "cool") {
    data <- load_hic_cool(file_path = file_path, ...)
  } else {
    stop(str_interp("Invalid format ${format}"))
  }
  data
}



#' Dump a Hi-C object as .hic format
#'
#' This function invokes Juicer tools for the dumping
#' @param hic_matrix A Hi-C object
#' @param file_path Path to the output .hic file
#' @param java Path to JVM. Default is \code{java}
#' @param ref_genome Reference genome \code{hic_matrix} is using. Default is \code{hg19}
#' @export
dump_juicer_hic <- function(hic_matrix, file_path, juicertools, java = "java", ref_genome = "hg19") {
  stopifnot(ref_genome == "hg19")

  juicer_short_path <- tempfile(fileext = ".short")
  tryCatch({
    dump_juicer_short(hic_matrix, file_path = juicer_short_path)

    cmd <- str_interp("${java} -jar ${juicertools} pre ${juicer_short_path} ${file_path} ${ref_genome}")
    retcode <- system(cmd)
    if (retcode != 0) {
      stop(str_interp("Error in creating .hic , RET: ${retcode}, CMD: ${cmd}"))
    }
  }, finally = {
    unlink(juicer_short_path)
  })
}


#' @export
dump_juicer_short <- function(hic_matrix, file_path) {
  # 0 22 16000000 0 0 22 16000000 1 95
  hic_matrix %>%
    mutate(str1 = 0, frag1 = 0, str2 = 0, frag2 = 1) %>%
    select(str1, chrom1, pos1, frag1, str2, chrom2, pos2, frag2, score) %>%
    na.omit() %>%
    write_delim(file = file_path, col_names = FALSE, delim = " ")
}


#' @export
dump_hic_genbed <- function(hic_matrix, file_path) {
  resol <- attr(hic_matrix, "resol")
  stopifnot(is.integer(resol) && resol > 0)

  hic_matrix %>%
    rename(start1 = pos1, start2 = pos2) %>%
    mutate(end1 = start1 + resol, end2 = start2 + resol) %>%
    select(chrom1, start1, end1, chrom2, start2, end2, score, everything()) %>%
    bioessentials::write_genbed(file = file_path)
}


convert_matrix_hic <- function(mat, chrom, resol, pos_start) {
  stopifnot(nrow(mat) == ncol(mat))
  dim <- nrow(mat)
  stopifnot(length(chrom) == 1)
  stopifnot(length(resol) == 1 & resol > 0)

  all_pos <- (1:dim - 1) * resol + pos_start
  expand_grid(pos1 = all_pos,
              pos2 = all_pos) %>%
    mutate(
      chrom1 = chrom,
      chrom2 = chrom,
      pos1 = as.integer(pos1),
      pos2 = as.integer(pos2),
      score = as.numeric(mat)
    ) %>%
    filter(pos1 <= pos2) %>%
    na.omit() %>%
    set_attr("resol", resol) %>%
    select(chrom1, pos1, chrom2, pos2, score)
}

#' Convert a Hi-C map to a naive matrix
#' @export
convert_hic_matrix <- function(hic_matrix, chrom = NULL) {
  if (is.null(chrom)) {
    # Infer chrom from input
    chrom <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2))
    stopifnot(length(chrom) == 1)
  }

  resol <- attr(hic_matrix, "resol")
  stopifnot(resol > 0)

  # Single-chromosome matrix
  hic_matrix %<>%
    filter(chrom1 == chrom & chrom2 == chrom)

  min_pos <- min(c(hic_matrix$pos1, hic_matrix$pos2))
  max_pos <- max(c(hic_matrix$pos1, hic_matrix$pos2))
  hic_matrix %<>%
    mutate(x = (pos1 - min_pos) %/% resol + 1,
           y = (pos2 - min_pos) %/% resol + 1)

  mat_dim <- (max_pos - min_pos) %/% resol + 1
  mat <- matrix(rep(NA, mat_dim * mat_dim), nrow = mat_dim)
  mat[hic_matrix %>% select(x, y) %>% as.matrix()] <- hic_matrix$score
  mat[hic_matrix %>% select(y, x) %>% as.matrix()] <- hic_matrix$score

  mat_pos <- (1:mat_dim - 1) * resol + min_pos
  mat_header <- mat_pos %>% map_chr(~ str_interp("chr${chrom}-$[d]{.}"))
  colnames(mat) <- mat_header
  rownames(mat) <- mat_header
  mat
}


#' @export
dump_hic_matrix <- function(mat, file_path, chrom = NULL) {
  # Detect all-NA rows/columns
  na_band <- apply(
    mat,
    MARGIN = 1,
    FUN = function(v)
      all(is.na(v))
  )
  mat <- mat[!na_band,!na_band]

  header <- paste(c("", colnames(mat)), collapse = "\t")
  body <- 1:nrow(mat) %>% map_chr(function(idx) {
    pos_str <- colnames(mat)[idx]
    row_str <- paste(mat[idx,], collapse = "\t")
    paste0(pos_str, "\t", row_str)
  })
  c(header, body) %>% write_lines(file = file_path)
}


#' Dump Hi-C data as .cool
#' @export
dump_cool <- function(hic_matrix, file_path, juicertools, java = "java", executable = "hicConvertFormat") {
  # First dump as as a temporary genbed file, then call hicConvertFormat
  temp_hic <- tempfile(fileext = ".hic")

  tryCatch({
    resol <- attr(hic_matrix, "resol")
    hic_matrix %>%
      dump_juicer_hic(file_path = temp_hic,
                      juicertools = juicertools,
                      java = java)

    # Output name
    # For single resolution .cool dataset, hicConvertFormat will implicit change the output
    # file name. For example: foo.cool -> foo_500000.cool
    output_name <- str_interp("${str_sub(file_path, end = -6L)}_${resol}.cool")

    cmd <- str_interp(
      paste0(
        "${executable} -m ${temp_hic} --inputFormat hic ",
        "-r ${resol} ",
        "-o ${file_path} --outputFormat cool"
      )
    )
    cat(cmd, "\n")
    system(cmd)
    file.rename(output_name, file_path)
  }, finally = {
    unlink(temp_hic)
  })
}


#' Dump Hi-C data to file
#'
#' @export
dump_hic <- function(hic_matrix, file_path, format = NULL, ...) {
  if (is.null(format)) {
    format <- guess_format(file_path)
  }
  if (format == "juicer_short") {
    dump_juicer_short(hic_matrix, file_path)
  } else if (format == "juicer_hic") {
    dump_juicer_hic(hic_matrix, file_path, ...)
  } else if (format == "genbed") {
    dump_hic_genbed(hic_matrix, file_path)
  } else {
    stop(str_interp("Invalid format ${format}"))
  }
}
