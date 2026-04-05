
# data.R — load and preprocess data


#' Load and preprocess dose-response data
#'
#' Reads a CSV and coerces columns to expected types. Negative response values
#' are clamped to 0 before any downstream fitting.
#'
#' @details
#' The input CSV must contain the following columns:
#' - `dose` — numeric drug concentration
#' - `response` — numeric viability or death response (0–100 scale expected)
#' - `cell_line` — cell line identifier
#' - `ancestry` — population/ancestry label
#' - `feature` — translocation or mutation status (use `"none"` if absent)
#'
#' @param filepath Character. Path to the input CSV file.
#'
#' @return A data frame with `dose` and `response` as numeric, and
#'   `cell_line`, `ancestry`, `feature` as factors.
#' @export
#'
#' @examples
#' \dontrun{
#' raw <- load_data("drug_response_data.csv")
#' }
load_data <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("File not found: '%s'", filepath), call. = FALSE)
  }
  
  data <- utils::read.csv(filepath)
  
  validate_columns(
    data,
    c("dose", "response", "cell_line", "ancestry", "feature"),
    arg_name = "filepath"
  )
  
  data$dose              <- as.numeric(data$dose)
  data$response          <- as.numeric(data$response)
  data$response[data$response < 0] <- 0
  data$cell_line         <- factor(data$cell_line)
  data$ancestry          <- factor(data$ancestry)
  data$feature     <- factor(data$feature)
  
  data
}


#' Aggregate replicate measurements by cell line and dose
#'
#' Collapses technical or biological replicates to per-group means and standard
#' deviations. This aggregated output is the expected input for all downstream
#' fitting and plotting functions.
#'
#' @param data A data frame from [load_data()].
#'
#' @return A data frame with columns `cell_line`, `dose`, `ancestry`,
#'   `feature`, `mean_response`, `sd`, and `n` (replicate count).
#' @export
#'
#' @examples
#' \dontrun{
#' raw <- load_data("drug_response_data.csv")
#' agg <- combine_reps(raw)
#' }
combine_reps <- function(data) {
  validate_columns(
    data,
    c("cell_line", "dose", "ancestry", "feature", "response")
  )
  
  data |>
    dplyr::group_by(cell_line, dose, ancestry, feature) |>
    dplyr::summarise(
      sd            = stats::sd(response, na.rm = TRUE),
      mean_response = mean(response, na.rm = TRUE),
      n             = dplyr::n(),
      .groups       = "drop"
    ) |>
    dplyr::mutate(cell_line = factor(cell_line))
}
