# plot_drc.R — dose-response curve visualizations

#' Plot dose-response curves colored by cell line
#'
#' Overlays observed data points (with SD error bars) and smooth fitted curves
#' for all cell lines on a single log10-scaled dose axis.
#'
#' @param agg_data Aggregated data frame from [combine_reps()].
#' @param pred_data Prediction data frame from [predict_drc()].
#' @param title Character. Plot title. Default `""`.
#' @param colors Optional ggplot2 color scale (
#'   `ggsci::scale_colour_nejm()`). Default `NULL`.
#' @param y_label Character. Y-axis label. Default `"Effect (%)"`.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plot_drc(agg, preds, title = "Drug Response")
#' }
plot_drc <- function(agg_data, pred_data, title = "", 
                     y_label = "Effect (%)" , colors = NULL) {
  validate_columns(agg_data, c("dose", "mean_response", "sd", "cell_line"))
  validate_columns(pred_data, c("dose", "predicted_response", "cell_line"))
  
  p <- ggplot2::ggplot(
    agg_data,
    ggplot2::aes(x = dose, y = mean_response, color = cell_line)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_response - sd, ymax = mean_response + sd),
      width = 0.1, alpha = 0.6
    ) +
    ggplot2::geom_line(
      data = pred_data,
      ggplot2::aes(x = dose, y = predicted_response, color = cell_line),
      linewidth = 1
    ) +
    ggplot2::scale_x_log10(name = "Dose [\u00b5M]") +
    ggplot2::scale_y_continuous(name = y_label, limits = c(0, 110)) +
    drc_theme() +
    ggplot2::ggtitle(title)
  
  if (!is.null(colors)) p <- p + colors
  p
}


#' Plot dose-response curves grouped by ancestry and faceted by feature
#'
#' Color and linetype both encode ancestry. Individual cell lines are labeled
#' at the highest dose point using [ggrepel::geom_text_repel()]. Panels are
#' faceted by feature group.
#'
#' @param agg_data Aggregated data frame from [combine_reps()].
#' @param pred_data Prediction data frame from [predict_drc()].
#' @param title Character. Plot title. Default `""`.
#' @param ancestry_colors Named character vector of hex colors per ancestry
#'   group. Default [ANCESTRY_COLORS].
#' @param ancestry_linetypes Named character vector of linetypes per ancestry
#'   group. Default [ANCESTRY_LINETYPES].
#' @param label_size Numeric. Font size for cell line labels. Default `3`.
#' @param y_label Character. Y-axis label. Default `"Effect (%)"`.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plot_drc_anc(agg, preds, title = "Drug Response by Ancestry")
#' }
plot_drc_anc <- function(agg_data,
                         pred_data,
                         title = "",
                         y_label = "Effect (%)",
                         ancestry_colors    = ANCESTRY_COLORS,
                         ancestry_linetypes = ANCESTRY_LINETYPES,
                         label_size         = 3) {
  validate_columns(agg_data,
                   c("dose", "mean_response", "sd", "ancestry",
                     "cell_line", "feature"))
  validate_columns(pred_data,
                   c("dose", "predicted_response", "ancestry",
                     "cell_line", "feature"))
  
  label_data <- pred_data |>
    dplyr::group_by(cell_line, feature) |>
    dplyr::slice_max(dose, n = 1) |>
    dplyr::ungroup()
  
  ggplot2::ggplot(
    agg_data,
    ggplot2::aes(
      x        = dose,
      y        = mean_response,
      color    = ancestry,
      linetype = ancestry,
      group    = cell_line
    )
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_response - sd, ymax = mean_response + sd),
      width = 0.1, alpha = 0.6
    ) +
    ggplot2::geom_line(
      data = pred_data,
      ggplot2::aes(
        x        = dose,
        y        = predicted_response,
        color    = ancestry,
        linetype = ancestry
      ),
      linewidth = 1
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(
        x     = dose,
        y     = predicted_response,
        label = cell_line,
        color = ancestry
      ),
      size         = label_size,
      fontface     = "bold",
      show.legend  = FALSE,
      box.padding  = 0.4,
      max.overlaps = Inf,
      nudge_y      = 5
    ) +
    ggplot2::facet_wrap(~feature, scales = "free") +
    ggplot2::scale_x_log10(name = "Dose [\u00b5M]") +
    ggplot2::scale_y_continuous(name = y_label, limits = c(0, 110)) +
    ggplot2::scale_color_manual(values = ancestry_colors) +
    ggplot2::scale_linetype_manual(values = ancestry_linetypes) +
    drc_theme(axis_title_size = 16, bold_axis_titles = TRUE) +
    ggplot2::ggtitle(title)
}