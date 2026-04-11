# predict.R — Generate smooth predictions for plotting

#' Generate smooth model predictions across the dose range
#'
#' Produces a tidy data frame of predicted responses for plotting
#' fitted curves. Cell line metadata (`ancestry`, `feature`) is joined
#' from `agg_data` so it propagates automatically to all downstream plot
#' functions.
#'
#' @param fits Named list of `drc` model objects from [fitALL()].
#' @param agg_data Aggregated data frame from [combineReps()].
#' @param n_points Integer. Number of evenly-spaced dose points per cell line
#'   over the observed dose range. Default `200`.
#'
#' @return A data frame with columns `dose`, `cell_line`,
#'   `predicted_response`, `ancestry`, and `feature`.
#' @export
#'
#' @examples
#' \dontrun{
#' agg   <- combineReps(load_data("data.csv"))
#' fits  <- fitALL(agg, unique(agg$cell_line))
#' preds <- predictDRC(fits, agg)
#' }
predictDRC <- function(fits, agg_data, n_points = 200) {
  validateCols(agg_data,
                   c("cell_line", "dose", "ancestry", "feature"),
                   arg_name = "agg_data")
  
  meta <- agg_data |>
    dplyr::select(cell_line, ancestry, feature) |>
    dplyr::distinct()
  
  dose_seq <- seq(
    min(agg_data$dose),
    max(agg_data$dose),
    length.out = n_points
  )
  
  tidyr::expand_grid(dose = dose_seq, cell_line = names(fits)) |>
    dplyr::mutate(
      predicted_response = purrr::map2_dbl(
        cell_line, dose,
        ~as.numeric(
          stats::predict(fits[[.x]], newdata = data.frame(dose = .y))
        )
      )
    ) |>
    dplyr::left_join(meta, by = "cell_line")
}