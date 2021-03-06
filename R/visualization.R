#' Plot the A/B compartment profile
#'
#' @export
#' @return A ggplot2 object showing the compartments
plot_compartment <-
  function(comps,
           chrom = NULL,
           full_scale = FALSE) {
    assert_that(is(comps, "GRanges") || is(comps, "data.frame"))
    if (is(comps, "GRanges"))
      comps %<>% bedtorch::as.bedtorch_table()

    if (is.null(chrom)) {
      chrom <- unique(comps$chrom %>% as.character())
    }
    
    # Only deal with single-chromosome track
    assert_that(is_scalar_character(chrom))
    
    pick_idx <- comps$chrom == chrom
    comps <- comps[pick_idx]
    
    resol <- comps$start %>% diff() %>% min()
    p <- comps[!is.na(score)] %>%
      mutate(compartment = factor(ifelse(score > 0, "open", "close"), levels = c("open", "close"))) %>%
      ggplot(aes(x = start, y = score, fill = compartment)) +
      geom_col(width = resol * .75) +
      xlab("coordinate")

    if (full_scale) {
      p + scale_x_continuous(labels = scales::comma, limits = c(0, max(comps$start)))
    } else {
      p + scale_x_continuous(labels = scales::comma)
    }
  }


#' Visualize Hi-C data
#'
#' Show a heatmap of the Hi-C data
#'
#' @param chrom Indicate which chromosome to process. If `NULL`, `hic_matrix`
#'   should contain only one chromosome, which will be used in the
#'   visualization.
#' @export
plot_hic_matrix <- function(hic_matrix,
                            chrom = NULL,
                            control_hic_matrix = NULL,
                            # control_gm = NULL,
                            scale_factor = c(1, 1),
                            norm = c(0, 1),
                            transform = NULL,
                            # symfill = TRUE,
                            missing_value = NULL,
                            n.breaks = NULL,
                            color_palette = "viridis",
                            matrix = "observed",
                            gamma = 2.3,
                            tile_outline = NULL) {
  if (is_null(chrom))
    chrom <- unique(c(hic_matrix$chrom1, hic_matrix$chrom2) %>% as.character())
  assert_that(is_scalar_character(chrom))

  resol <- attr(hic_matrix, "resol")
  stopifnot(resol > 0)

  if (is(hic_matrix, "data.table"))
    hic_matrix <- as_tibble(hic_matrix)

  # Ensure the input is upper trangular
  hic_matrix %<>% filter(chrom1 == chrom & chrom2 == chrom)
  stopifnot(with(hic_matrix, all(pos1 <= pos2)))
  if (!is.null(control_hic_matrix)) {
    control_hic_matrix %<>% filter(chrom1 == chrom & chrom2 == chrom)
    stopifnot(with(control_hic_matrix, all(pos1 <= pos2)))
    stopifnot(attr(control_hic_matrix, "resol") == resol)
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

  full_matrix <- if (is.null(control_hic_matrix)) {
    full_matrix <- build_full_matrix(hic_matrix, hic_matrix)
  } else {
    full_matrix <- build_full_matrix(control_hic_matrix, hic_matrix)
  }

  gr <- local({
    pos <- with(full_matrix, c(pos1, pos2))
    c(min(pos), max(pos))
  })
  gr_str <-
    str_interp("${chrom}:$[d]{gr[1] + 1}-$[d]{gr[2] + resol}")

  # Transform & normalization
  local({
    if (is.null(transform))
      score <- full_matrix$score
    else if (transform == "log10") {
      score <- full_matrix$score
      score[score < 1] <- 1
      score <- log10(score)
    } else {
      stop(transform)
    }

    # browser()
    m1 <- min(score)
    m2 <- max(score)
    score <- (score - m1) / (m2 - m1)
    full_matrix$score <<- score * (norm[2] - norm[1]) + norm[1]
  })
  full_matrix %<>% mutate(score = score ** gamma)

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
    ylab(paste0(ifelse(
      is.null(control_hic_matrix), "", "Control: "
    ), gr_str)) +
    xlab(gr_str) +
    theme(legend.position = "none")

  if (matrix == "pearson")
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