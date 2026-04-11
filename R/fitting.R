
# fitting.R — dose-response model fitting

#' Fit a dose-response model for a single cell line
#'
#' Wraps [drc::drm()] for a single cell line subset of the aggregated data.
#' The default 4-parameter log-logistic model (`LL.4`) is appropriate for
#' sigmoidal dose-response relationships.
#'
#' @param data Aggregated data frame from [combineReps()].
#' @param cell_line_name Character. Name of the cell line to fit.
#' @param model A `drc` model function. Default `drc::LL.4()`. Other options
#'   include `drc::LL.3()` (fixed lower asymptote at 0) or `drc::W1.4()`.
#'
#' @return A `drc` model object.
#' @export
#'
#' @examples
#' \dontrun{
#' agg  <- combineReps(load_data("data.csv"))
#' fit  <- fitCL(agg, "SKMM-1")
#' }
fitCL <- function(data, cell_line_name, model = drc::LL.4()) {
  validateCols(data, c("cell_line", "dose", "mean_response"))
  
  subset_data <- subset(data, cell_line == cell_line_name)
  
  if (nrow(subset_data) == 0) {
    stop(
      sprintf("No data found for cell line: '%s'", cell_line_name),
      call. = FALSE
    )
  }
  
  drc::drm(mean_response ~ dose, data = subset_data, fct = model)
}


#' Fit dose-response models for all cell lines
#'
#' Iterates over a character vector of cell line names and calls
#' [fitCL()] for each, returning a named list of `drc` model objects.
#'
#' @param data Aggregated data frame from [combineReps()].
#' @param cell_lines Character vector of cell line names. Typically
#'   `unique(agg$cell_line)`.
#' @param model A `drc` model function passed to [fitCL()].
#'   Default `drc::LL.4()`.
#'
#' @return A named list of `drc` model objects, one entry per cell line.
#' @export
#'
#' @examples
#' \dontrun{
#' agg  <- combineReps(load_data("data.csv"))
#' fits <- fitALL(agg, unique(agg$cell_line))
#' }
fitALL <- function(data, cell_lines, model = drc::LL.4()) {
  if (length(cell_lines) == 0) {
    stop("'cell_lines' must contain at least one cell line name.", call. = FALSE)
  }
  
  fits <- purrr::map(cell_lines, ~fitCL(data, .x, model = model))
  stats::setNames(fits, as.character(cell_lines))
}
