# metrics.R — pharmacological metrics: IC50, AUC, hill slope, summary

#' Extract IC50 from a single fitted model
#'
#' Uses [drc::ED()] with delta-method confidence intervals to estimate the
#' dose producing 50% effect.
#'
#' @param fit A `drc` model object from [fitCL()].
#' @param cell_line Character. Cell line name appended as a column.
#' @param ancestry Character or `NA`. Ancestry label appended as a column.
#' @param feature Character or `NA`. Feature label appended as a column.
#'
#' @return A one-row data frame with columns `Estimate`, `SE`, `Lower`,
#'   `Upper`, `cell_line`, `ancestry`, `feature`.
#' @export
extractIC50 <- function(fit, cell_line,
                         ancestry = NA, feature = NA) {
  ic50_mat <- suppressWarnings(drc::ED(fit, 50, interval = "delta"))
  ic50_df  <- as.data.frame(ic50_mat)
  colnames(ic50_df)[1:4] <- c("Estimate", "SE", "Lower", "Upper")
  
  ic50_df$cell_line <- cell_line
  ic50_df$ancestry  <- ancestry
  ic50_df$feature   <- feature
  
  ic50_df
}


#' Extract IC50 estimates for all fitted cell lines
#'
#' Joins cell line metadata from `agg_data` and calls [extractIC50()] for
#' each model in the list.
#'
#' @param fits Named list of `drc` model objects from [fitALL()].
#' @param agg_data Aggregated data frame from [combineReps()], used to recover
#'   `ancestry` and `feature` metadata per cell line.
#'
#' @return A data frame with columns `Estimate`, `SE`, `Lower`, `Upper`,
#'   `cell_line`, `ancestry`, `feature`.
#' @export
extractALL <- function(fits, agg_data) {
  validateCols(agg_data, c("cell_line", "ancestry", "feature"),
                   arg_name = "agg_data")
  
  cell_line_info <- agg_data |>
    dplyr::select(cell_line, ancestry, feature) |>
    dplyr::distinct()
  
  purrr::map_df(names(fits), function(cl) {
    info <- dplyr::filter(cell_line_info, cell_line == cl)
    extractIC50(fits[[cl]], cl,
                 ancestry = info$ancestry[1],
                 feature  = info$feature[1])
  })
}


#' calc area under the dose-response curve (AUC)
#'
#' computes AUC using the trapezoidal rule over the log10-transformed dose
#' axis. AUC complements IC50 by summarizing the full curve shape — it is
#' defined even when curves do not reach 50% inhibition.
#'
#' A higher AUC indicates greater drug sensitivity. AUC is normalized to
#' \[0, 1\] by default by dividing by the total log-dose range.
#'
#' @param pred_data Prediction data frame from [predictDRC()].
#' @param normalize Logical. If `TRUE` (default), AUC is divided by the total
#'   log10 dose range so values fall in \[0, 1\].
#'
#' @return A data frame with columns `cell_line`, `ancestry`, `feature`,
#'   and `AUC`.
#' @export
#'
#' @examples
#' \dontrun{
#' auc <- calcAUC(preds)
#' }
calcAUC <- function(pred_data, normalize = TRUE) {
  validateCols(pred_data,
                   c("cell_line", "dose", "predicted_response",
                     "ancestry", "feature"))
  
  meta <- pred_data |>
    dplyr::select(cell_line, ancestry, feature) |>
    dplyr::distinct()
  
  auc_vals <- pred_data |>
    dplyr::group_by(cell_line) |>
    dplyr::arrange(dose, .by_group = TRUE) |>
    dplyr::summarise(
      AUC = {
        log_dose <- log10(dose)
        resp     <- predicted_response / 100
        
        n       <- length(log_dose)
        widths  <- diff(log_dose)
        heights <- (resp[-n] + resp[-1]) / 2
        raw_auc <- sum(widths * heights)
        
        if (normalize) {
          raw_auc / (max(log_dose) - min(log_dose))
        } else {
          raw_auc
        }
      },
      .groups = "drop"
    )
  
  dplyr::left_join(auc_vals, meta, by = "cell_line")
}


#' Extract Hill slope from fitted dose-response models
#'
#' The Hill slope (the `b` parameter of the LL.4 model) describes the
#' steepness of the dose-response transition. A steeper slope indicates
#' a sharper switch between no-effect and full-effect doses.
#'
#' @param fits Named list of `drc` model objects from [fitALL()].
#' @param agg_data Aggregated data frame from [combineReps()].
#'
#' @return A data frame with columns `cell_line`, `ancestry`, `feature`,
#'   `hill_slope`, and `hill_slope_se`.
#' @export
#'
#' @examples
#' \dontrun{
#' hills <- extractHillSlope(fits, agg)
#' }
extractHillSlope <- function(fits, agg_data) {
  validateCols(agg_data, c("cell_line", "ancestry", "feature"),
                   arg_name = "agg_data")
  
  cell_line_info <- agg_data |>
    dplyr::select(cell_line, ancestry, feature) |>
    dplyr::distinct()
  
  purrr::map_df(names(fits), function(cl) {
    coef_mat <- stats::coef(summary(fits[[cl]]))
    hill_row <- coef_mat[grepl("^b", rownames(coef_mat)), , drop = FALSE]
    info     <- dplyr::filter(cell_line_info, cell_line == cl)
    
    data.frame(
      cell_line     = cl,
      ancestry      = info$ancestry[1],
      feature       = info$feature[1],
      hill_slope    = hill_row[1, "Estimate"],
      hill_slope_se = hill_row[1, "Std. Error"],
      stringsAsFactors = FALSE
    )
  })
}


#' Summarize all pharmacological metrics for fitted models
#'
#' Convenience wrapper combining IC50, AUC, and Hill slope into a single
#' tidy data frame. Recommended starting point for downstream statistical
#' analysis and omics correlation.
#'
#' @param fits Named list of `drc` model objects from [fitALL()].
#' @param agg_data Aggregated data frame from [combineReps()].
#' @param pred_data Prediction data frame from [predictDRC()].
#' @param normalize_auc Logical. Passed to [calcAUC()]. Default `TRUE`.
#'
#' @return A data frame with columns `cell_line`, `ancestry`, `feature`,
#'   `IC50`, `IC50_lower`, `IC50_upper`, `AUC`, `hill_slope`.
#' @export
#'
#' @examples
#' \dontrun{
#' metrics <- summarizeDRC(fits, agg, preds)
#' }
summarizeDRC <- function(fits, agg_data, pred_data, normalize_auc = TRUE) {
  ic50 <- extractALL(fits, agg_data) |>
    dplyr::select(cell_line, ancestry, feature,
                  IC50 = Estimate, IC50_lower = Lower, IC50_upper = Upper)
  
  auc  <- calcAUC(pred_data, normalize = normalize_auc) |>
    dplyr::select(cell_line, AUC)
  
  hill <- extractHillSlope(fits, agg_data) |>
    dplyr::select(cell_line, hill_slope)
  
  ic50 |>
    dplyr::left_join(auc,  by = "cell_line") |>
    dplyr::left_join(hill, by = "cell_line")
}