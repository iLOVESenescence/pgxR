
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
#' @param col_map Named list mapping pgxR's expected column names to the 
#' actual column names in your CSV. Only specify columns that differ.
#' For example, if your CSV has `translocation` instead of `feature`:
#'   `col_map = list(feature = "translocation")`.
#'   Or if ancestry is stored as `population`:
#'   `col_map = list(ancestry = "population", feature = "translocation")`.
#'   
#' @return A data frame with `dose` and `response` as numeric, and
#'   `cell_line`, `ancestry`, `feature` as factors.
#' @export
#'
#' @examples
#' \dontrun{
#' raw <- loadData("drug_response_data.csv")
#' raw <- loadData("drug_response_data.csv", col_map = list(condition = "treatment_group"))
#' }
loadData <- function(filepath, col_map = NULL) {
  if (!file.exists(filepath)) {
    stop(sprintf("File not found: '%s'", filepath), call. = FALSE)
  }
  
  data <- utils::read.csv(filepath)
  
  #rename user columns to pgxR expected names
  if (!is.null(col_map)) {
    for (pgxr_name in names(col_map)) {
      user_name <- col_map[[pgxr_name]]
      if (!user_name %in% colnames(data)) {
        stop(
          sprintf(
            "col_map error: column '%s' not found in data. Available columns: %s",
            user_name,
            paste(colnames(data), collapse = ", ")
          ),
          call. = FALSE
        )
      }
      colnames(data)[colnames(data) == user_name] <- pgxr_name
    }
  }
  
  validateCols(
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
  
  #additional optional columns that will be factored if present, but ignore if not
  optional_factor_cols <- c("condition")
  for (col in optional_factor_cols) {
    if (col %in% colnames(data)) {
      data[[col]] <- factor(data[[col]])
    }
  }
  data
}


#' Aggregate replicate measurements by cell line and dose
#'
#' Collapses technical or biological replicates to per-group means and standard
#' deviations. This aggregated output is the expected input for all downstream
#' fitting and plotting functions.
#'
#' Extra columns (e.g. replicate identifiers, experiment labels) are ignored
#' during aggregation — only `cell_line`, `dose`, `ancestry`, `feature`, and
#' `response` are used.
#'
#' @param data A data frame from [loadData()].
#' @param replicate_col Character or `NULL`. Name of the column identifying
#'   replicates (e.g. `"experiment"`, `"replicate"`, `"run"`). If provided,
#'   validates the column exists before aggregating. If `NULL` (default),
#'   replicates are collapsed implicitly across all rows sharing the same
#'   `cell_line`, `dose`, `ancestry`, and `feature`.
#'
#' @return A data frame with columns `cell_line`, `dose`, `ancestry`,
#'   `feature`, `mean_response`, `sd`, and `n` (replicate count).
#' @export
#'
#' @examples
#' \dontrun{
#' # implicit replicate collapsing
#' agg <- combineReps(raw)
#'
#' # explicit replicate column
#' agg <- combineReps(raw, replicate_col = "experiment")
#' }
combineReps <- function(data, replicate_col = NULL) {
  validateCols(
    data,
    c("cell_line", "dose", "ancestry", "feature", "response")
  )
  
  if (!is.null(replicate_col)) {
    if (!replicate_col %in% colnames(data)) {
      stop(
        sprintf(
          "replicate_col '%s' not found in data. Available columns: %s",
          replicate_col,
          paste(colnames(data), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    message(sprintf(
      "Aggregating across %d replicates identified by '%s'.",
      length(unique(data[[replicate_col]])),
      replicate_col
    ))
  }
  #if condition is present include it in the grouping variable
  group_vars <- c("cell_line", "dose", "ancestry", "feature")
  if ("condition" %in% colnames(data)) {
    group_vars <- c(group_vars, "condition")
  }
  
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