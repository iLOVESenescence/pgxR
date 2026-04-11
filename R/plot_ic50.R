# plot_ic50.R — IC50 and AUC visualizations

#' Plot IC50 estimates as a dot plot with confidence intervals
#'
#' Displays point estimates and CI error bars for all cell lines ordered by
#' ascending IC50. The recommended default for overall sensitivity comparison.
#'
#' @param ic50_df IC50 data frame from [extractALL()].
#' @param title Character. Plot title. Default `""`.
#' @param colors Optional ggplot2 color scale. Default `NULL`.
#' @param units Character. Concentration units for y-axis label. Default `"nM"`.
#' 
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plotDot(ic50, title = "Drug IC50")
#' }
plotDot <- function(ic50_df, title = "",units = "nM", colors = NULL) {
  validateCols(ic50_df,
                   c("cell_line", "Estimate", "Lower", "Upper"))
  
  p <- ggplot2::ggplot(
    ic50_df,
    ggplot2::aes(
      x     = reorder(cell_line, Estimate),
      y     = Estimate,
      color = cell_line
    )
  ) +
    ggplot2::geom_point(size = 4, alpha = 0.9) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Lower, ymax = Upper),
      width = 0.2, alpha = 0.8
    ) +
    ggplot2::scale_y_continuous(name = "IC50 [%s]") +
    ggplot2::scale_x_discrete(name = "") +
    drc_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::ggtitle(title)
  
  if (!is.null(colors)) p <- p + colors
  p
}


#' Plot IC50 estimates faceted by genomic feature
#'
#' Point estimates and CI error bars per cell line, faceted by feature group
#' (e.g. translocation status, driver mutation). Ancestry encoded by linetype.
#'
#' @param ic50_df IC50 data frame from [extractALL()].
#' @param title Character. Plot title. Default `""`.
#' @param colors Optional ggplot2 fill scale. Default `NULL`.
#' @param ancestry_linetypes Named character vector mapping ancestry to
#'   linetypes. Default [ANCESTRY_LINETYPES].
#' @param units Character. Concentration units for y-axis label. Default `"nM"`.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plotFeature(ic50, title = "IC50 by Feature")
#' }
plotFeature <- function(ic50_df,
                         title = "",
                         units = "nM",
                         colors = NULL,
                         ancestry_linetypes = ANCESTRY_LINETYPES) {
  validateCols(ic50_df,
                   c("cell_line", "Estimate", "Lower", "Upper",
                     "ancestry", "feature"))
  
  p <- ggplot2::ggplot(
    ic50_df,
    ggplot2::aes(
      x        = reorder(cell_line, Estimate),
      y        = Estimate,
      color    = cell_line,
      linetype = ancestry
    )
  ) +
    ggplot2::geom_point(size = 4, alpha = 0.9) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Lower, ymax = Upper),
      width = 0.2, alpha = 0.8
    ) +
    ggplot2::facet_wrap(~feature) +
    ggplot2::scale_y_continuous(name = "IC50 [%s]") +
    ggplot2::scale_x_discrete(name = "") +
    ggplot2::scale_linetype_manual(values = ancestry_linetypes) +
    drc_theme() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::ggtitle(title)
  
  if (!is.null(colors)) p <- p + colors
  p
}


#' Plot IC50 estimates grouped by ancestry
#'
#' Displays IC50 point estimates and CI error bars colored by ancestry group.
#' Cell lines are ordered by ascending IC50. Directly addresses whether
#' ancestry-of-origin associates with drug sensitivity.
#'
#' @param ic50_df IC50 data frame from [extractALL()].
#' @param title Character. Plot title. Default `""`.
#' @param ancestry_colors Named character vector of hex colors per ancestry
#'   group. Default [ANCESTRY_COLORS].
#' @param units Character. Concentration units for y-axis label. Default `"nM"`.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plotAnc(ic50, title = "IC50 by Ancestry")
#' }
plotAnc <- function(ic50_df,
                     title = "",
                     units = "nM",
                    ancestry_colors = ANCESTRY_COLORS) {
  validateCols(ic50_df,
                   c("cell_line", "Estimate", "Lower", "Upper", "ancestry"))
  
  ggplot2::ggplot(
    ic50_df,
    ggplot2::aes(
      x     = reorder(cell_line, Estimate),
      y     = Estimate,
      color = ancestry
    )
  ) +
    ggplot2::geom_point(size = 4, alpha = 0.9) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Lower, ymax = Upper),
      width = 0.2, alpha = 0.8
    ) +
    ggplot2::scale_y_continuous(name = "IC50 [%s]") +
    ggplot2::scale_x_discrete(name = "") +
    ggplot2::scale_color_manual(values = ancestry_colors) +
    drc_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::ggtitle(title)
}


#' Plot AUC estimates across cell lines
#'
#' Lollipop plot of normalized AUC values ordered by ascending sensitivity.
#' Higher AUC indicates greater drug sensitivity. Complements IC50 for cell
#' lines where curves do not reach 50% inhibition.
#'
#' @param auc_df AUC data frame from [calcAUC()].
#' @param title Character. Plot title. Default `""`.
#' @param ancestry_colors Named character vector of hex colors per ancestry
#'   group. Default [ANCESTRY_COLORS].
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plotAUC(auc, title = "Melphalan AUC")
#' }
plotAUC <- function(auc_df, title = "", ancestry_colors = ANCESTRY_COLORS) {
  validateCols(auc_df, c("cell_line", "AUC", "ancestry"))
  
  ggplot2::ggplot(
    auc_df,
    ggplot2::aes(
      x     = reorder(cell_line, AUC),
      y     = AUC,
      color = ancestry
    )
  ) +
    ggplot2::geom_point(size = 4, alpha = 0.9) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x    = reorder(cell_line, AUC),
        xend = reorder(cell_line, AUC),
        y    = 0,
        yend = AUC
      ),
      linewidth = 0.8, alpha = 0.6
    ) +
    ggplot2::scale_y_continuous(name = "AUC (normalized)", limits = c(0, 1)) +
    ggplot2::scale_x_discrete(name = "") +
    ggplot2::scale_color_manual(values = ancestry_colors) +
    drc_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::ggtitle(title)
}

#' Plot sensitivity metric by group with statistical comparisons
#'
#' Displays individual cell line values as jittered points with mean and
#' confidence interval bars. Pairwise p-values are annotated between groups
#' using [ggpubr::stat_compare_means()]. Works for any grouping variable
#' (ancestry, feature, mutation status etc.) and any sensitivity metric
#' (IC50, AUC, hill_slope).
#'
#' @param df A data frame containing the sensitivity metric and grouping
#'   column. Typically output from [summarizeDRC()] or [extractALL()].
#' @param group_col Character. Name of the grouping column (e.g. `"ancestry"`,
#'   `"feature"`).
#' @param metric Character. Name of the sensitivity metric column to plot.
#'   Default `"IC50"`.
#' @param colors Optional ggplot2 fill/color scale. Default `NULL`.
#' @param comparisons Optional list of group pairs for p-value annotation.
#'   If `NULL` (default) all pairwise comparisons are shown.
#' @param y_label Character. Y-axis label. Default uses `metric` value.
#' @param title Character. Plot title. Default `""`.
#' @param units Character. Concentration units for y-axis label when
#'   `metric = "IC50"`. Default `"nM"`.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' # IC50 by ancestry
#' plotSensitivity(metrics, group_col = "ancestry", metric = "IC50")
#'
#' # AUC by feature
#' plotSensitivity(metrics, group_col = "feature", metric = "AUC")
#' }
plotSensitivity <- function(df,
                             group_col,
                             metric = "IC50",
                             units = "nM",
                             colors = NULL,
                             comparisons = NULL,
                             y_label = NULL,
                             title = "") {
  validateCols(df, c(group_col, metric))
  
  if (is.null(y_label)) {
    y_label <- switch(metric,
                      "IC50" = sprintf("IC50 [%s]", units),
                      "AUC" = "AUC (normalized)",
                      "hill_slope" = "Hill Slope",
                      metric        
    )
  }
  
  # build all pairwise comparisons if not supplied
  if (is.null(comparisons)) {
    groups      <- unique(df[[group_col]])
    comparisons <- utils::combn(as.character(groups), 2,
                                simplify = FALSE)
  }
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x    = .data[[group_col]],
      y    = .data[[metric]],
      fill = .data[[group_col]]
    )
  ) +
    ggplot2::geom_jitter(
      width  = 0.15,
      size   = 3,
      alpha  = 0.7,
      ggplot2::aes(color = .data[[group_col]])
    ) +
    ggplot2::stat_summary(
      fun      = mean,
      fun.min  = function(x) mean(x) - stats::sd(x),
      fun.max  = function(x) mean(x) + stats::sd(x),
      geom     = "errorbar",
      width    = 0.2,
      linewidth = 0.8
    ) +
    ggplot2::stat_summary(
      fun  = mean,
      geom = "crossbar",
      width = 0.4,
      linewidth = 0.6,
      fill = "white",
      alpha = 0.6
    ) +
    ggpubr::stat_compare_means(
      comparisons = comparisons,
      method      = "wilcox.test",
      label       = "p.format",
      tip.length  = 0.01
    ) +
    ggplot2::scale_y_continuous(name = y_label) +
    ggplot2::scale_x_discrete(name = "") +
    drc_theme() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::ggtitle(title)
  
  if (!is.null(colors)) p <- p + colors
  p
}