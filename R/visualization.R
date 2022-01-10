#' Plot the A/B compartment profile
#'
#' @param score_col name of the column for compartment scores
#' @param full_scale whether to plot the track for the full extent of the chromosome
#' @export
#' @return A ggplot2 object showing the compartments
plot_compartment <-
  function(comps,
           chrom = NULL,
           score_col = "score",
           full_scale = FALSE) {
    assert_that(is(comps, "GRanges") || is(comps, "data.frame"))
    if (!is(comps, "GRanges"))
      comps %<>% bedtorch::as.GenomicRanges()
    
    assert_that(score_col %in% colnames(mcols(comps)))
    comps$score <- mcols(comps)[[score_col]]
    comps <- comps[!is.na(comps$score)]
    
    if (is.null(chrom)) {
      chrom <- as.character(unique(seqnames(comps)))
    }
    # Only deal with single-chromosome track
    assert_that(is_scalar_character(chrom))
    comps <- sort(comps[seqnames(comps) == chrom])
    
    # If the last and the first bin have different sizes, it means the last bin has been trimmed
    resol <- as.integer(unique(GenomicRanges::width(comps)))
    if (length(resol) == 2)
      resol <- resol[1]
    assert_that(is_scalar_integer(resol))
    
    # if (type == "bar") {
    #   p <- comps %>%
    #     bedtorch::as.bedtorch_table() %>%
    #     mutate(compartment = factor(ifelse(score > 0, "open", "close"), levels = c("open", "close"))) %>%
    #     ggplot(aes(x = start, y = score, fill = compartment)) +
    #     geom_col(width = resol * .75) +
    #     xlab("coordinate") +
    #     theme(legend.position = "top")
    #
    #   if (full_scale) {
    #     p <- p + scale_x_continuous(labels = scales::comma, limits = c(0, max(comps$start)))
    #   } else {
    #     p <- p + scale_x_continuous(labels = scales::comma)
    #   }
    #   return(p)
    # }
    # Fill the gap with 0s, this helps in interpolating scores
    comps_gapless <-
      GenomicRanges::tileGenome(
        seqlengths = GenomeInfoDb::keepSeqlevels(seqinfo(comps), chrom),
        tilewidth = resol,
        cut.last.tile.in.chrom = TRUE
      )
    comps_gapless$score <- 0
    hits <- GenomicRanges::findOverlaps(comps_gapless, comps)
    comps_gapless[S4Vectors::queryHits(hits)]$score <-
      comps[S4Vectors::subjectHits(hits)]$score
    
    # Interpolate data
    lin_interp <- function(x, y, length.out = 1000) {
      approx(x, y, xout = seq(min(x), max(x), length.out = length.out))$y
    }
    pos <- GenomicRanges::start(comps_gapless) + 0.5 * resol
    pos_interp <- lin_interp(pos, pos)
    score_interp <- lin_interp(pos, comps_gapless$score)
    df_interp <-
      data.table::data.table(pos = pos_interp, score = score_interp)
    
    # Make a grouping variable for each pos/neg segment
    cat_rle <- rle(sign(df_interp$score))
    df_interp$group <- rep.int(1:length(cat_rle$lengths), times = cat_rle$lengths)
    df_interp[, compartment := factor(ifelse(score > 0, "A", "B"))]
    df_interp <- df_interp[score != 0]
    
    p <- df_interp %>%
      ggplot(aes(
        x = pos,
        y = score,
        fill = compartment,
        group = group
      )) + geom_area() +
      xlab(NULL) + ylab(NULL) +
      theme(legend.position = "none")
    
    if (full_scale) {
      p + scale_x_continuous(labels = scales::comma, limits = c(0, max(df_interp$pos)))
    } else {
      p + scale_x_continuous(labels = scales::comma)
    }
  }


#' Visualize Hi-C data
#'
#' Show a heatmap of the Hi-C data
#'
#' @param chrom Indicate which chromosome to process. If `NULL`, `hic`
#'   should contain only one chromosome, which will be used in the
#'   visualization.
#' @export
plot_hic_matrix <- function(hic,
                            chrom = NULL,
                            control_hic = NULL,
                            # control_gm = NULL,
                            scale_factor = c(1, 1),
                            norm = c(0, 1),
                            transform = c("log10", "linear"),
                            # symfill = TRUE,
                            missing_value = NULL,
                            n.breaks = NULL,
                            color_palette = "viridis",
                            gamma = "auto",
                            tile_outline = NULL,
                            verbose = FALSE) {
  
  # Scale a numeric vector to range [a, b]
  scale_to_range <- function(x, a = 0, b = 1) {
    return((x - min(x, na.rm = TRUE)) /
             (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) * (b - a) + a)
  }
  
  # Calculate the best gamma value for the plot
  auto_gamma <- function(score_diag, score_off_diag) {
    billboard <- purrr::map_dfr(seq(1, 10, by = 0.1), function(gamma) {
      score_diag <- score_diag ** gamma
      score_off_diag <- score_off_diag ** gamma
      list(
        gamma = gamma,
        delta_m = abs(
          median(score_diag, na.rm = TRUE) - median(score_off_diag, na.rm = TRUE)
        ),
        sd_off_diag = sd(score_off_diag, na.rm = TRUE)
      )
    }) %>%
      mutate(
        delta_m = scale_to_range(delta_m),
        sd_off_diag = scale_to_range(sd_off_diag)
      ) %>%
      mutate(order_value = delta_m + sd_off_diag) %>%
      arrange(desc(order_value))
    
    return(billboard$gamma[1])
  }
  
  # Transform & scale to [0, 1]
  transcale <- function(score, transform, scale_min = 0, scale_max = 1) {
    if (transform == "log10") {
      # Ensure all log(score) values are non-negative
      # Exclude all non-positive values
      score[score <= 0] <- NA
      score[is.na(score)] <- NA
      factor <- 1 / min(score, na.rm = TRUE)
      score <- score * factor
      score[is.na(score)] <- 1
      score <- log10(score)
    } else if (transform != "linear") {
      stop(paste0("Invalid transform: ", transform), call. = FALSE)
    }
    
    m1 <- min(score)
    m2 <- max(score)
    return((score - m1) / (m2 - m1) * (scale_max - scale_min) + scale_min)
  }
  
  # Build the full matrix using the upper upper triangular
  build_full_matrix <- function(gm1, gm2) {
    rbind(
      gm1,
      gm2 %>% filter(pos1 != pos2) %>%
        mutate(
          s = pos1 + pos2,
          pos1 = (s - pos1),
          pos2 = (s - pos2)
        ) %>% select(-s)
    )
  }
  
  if (is_null(chrom))
    chrom <- unique(c(hic$chrom1, hic$chrom2) %>% as.character())
  assert_that(is_scalar_character(chrom))

  resol <- attr(hic, "resol")
  stopifnot(resol > 0)
  
  assert_that(is(hic, "ht_table"))
  assert_that(is.null(control_hic) || is(control_hic, "ht_table"))

  transform <- match.arg(transform)

  # Ensure the input is upper trangular
  hic %<>% filter(chrom1 == chrom & chrom2 == chrom)
  stopifnot(with(hic, all(pos1 <= pos2)))
  if (!is.null(control_hic)) {
    control_hic %<>% filter(chrom1 == chrom & chrom2 == chrom)
    stopifnot(with(control_hic, all(pos1 <= pos2)))
    stopifnot(attr(control_hic, "resol") == resol)
  }

  if (hic_type(hic) != "pearson")
    hic[, score := transcale(score, transform)]

  full_matrix <- if (is.null(control_hic)) {
    full_matrix <- build_full_matrix(hic, hic)
  } else {
    if (hic_type(control_hic) != "pearson")
      control_hic[, score := transcale(score, transform)]
    full_matrix <- build_full_matrix(hic, control_hic)
  }
  
  if (hic_type(hic) != "pearson") {
    if (gamma == "auto") {
      gamma <- auto_gamma(
        score_diag = filter(hic, pos1 == pos2)$score,
        score_off_diag = filter(hic, pos1 != pos2)$score
      )
      if (verbose)
        cat(paste0("Automatically estimate the best gamma: ", gamma))
    }
    full_matrix %<>% mutate(score = score ** gamma)
  }

  gr <- local({
    pos <- with(full_matrix, c(pos1, pos2))
    c(min(pos), max(pos))
  })
  gr_str <-
    str_interp("${chrom}:$[d]{gr[1] + 1}-$[d]{gr[2] + resol}")

  hic_colors <- list(
    colors = c(
      "#FFFFFF",
      "#FFF2F2",
      "#FFE8E8",
      "#FFCBCB",
      "#FFB3B3",
      "#FFA4A4",
      "#FF6565",
      "#FF0402"
    ),
    values = (c(0, 56, 95, 218, 265, 369, 603, 1033) / 1033)# ** 0.45
  )

  p <- ggplot(
    full_matrix %>% mutate(
      pos1 = pos1 / scale_factor[1],
      pos2 = pos2 / scale_factor[2]
    ),
    aes(x = pos1, y = pos2, fill = score)
  ) +
    (if (is.null(tile_outline))
      geom_tile()
     else
       geom_tile(color = tile_outline)) +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 1),
      n.breaks = n.breaks,
      position = "top"
    ) +
    scale_y_reverse(labels = scales::number_format(accuracy = 1),
                    n.breaks = n.breaks) +
    xlab(paste0(ifelse(
      is.null(control_hic), "", "Control: "
    ), gr_str)) +
    ylab(gr_str) +
    theme(legend.position = "none")

  if (hic_type(hic) == "pearson")
    p + scale_fill_gradientn(colors = c("blue", "black", "red"), values = c(0, 0.5, 1))
  else
    p + scale_fill_gradientn(colors = hic_colors$colors, values = hic_colors$values)
}


#' Get the scatter plot between two compartment tracks
#' @export
plot_compartment_scatter <-
  function(x,
           y,
           size = 0.1,
           alpha = 0.5,
           add = "reg.line",
           add.params = list(fill = "lightgray"),
           conf.int = TRUE,
           stat_cor_method = c("pearson", "spearman", "kendall"),
           ...) {
    x <- bedtorch::as.GenomicRanges(x)
    y <- bedtorch::as.GenomicRanges(y)
    assert_that(is_null(stat_cor_method) || is_character(stat_cor_method))
    if (!is_null(stat_cor_method))
      stat_cor_method <- match.arg(stat_cor_method[1], choices = c("pearson", "spearman", "kendall"))
    
    hits <- GenomicRanges::findOverlaps(x, y)
    score_x <- x[S4Vectors::queryHits(hits)]$score
    score_y <- y[S4Vectors::subjectHits(hits)]$score
    
    p <- tibble(x = x[S4Vectors::queryHits(hits)]$score,
           y = y[S4Vectors::subjectHits(hits)]$score) %>%
      ggpubr::ggscatter(
        x = "x",
        y = "y",
        size = size,
        alpha = alpha,
        add = add,
        # Add regressin line
        add.params = add.params,
        # Customize reg. line
        conf.int = conf.int,
        ...
      ) 
    
    if (!is_null(stat_cor_method))
      p <- p + ggpubr::stat_cor(method = stat_cor_method)
    
    p
  }