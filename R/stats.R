# stats.R — statistical comparisons across groups

#' Compare sensitivity metric across any grouping variable
#'
#' Runs a one-way ANOVA with Tukey HSD post-hoc test comparing a sensitivity
#' metric across levels of any grouping column. Works for ancestry, feature,
#' or any other categorical variable in your data.
#'
#' @param metrics_df A data frame from [summarizeDRC()] or [extractALL()].
#' @param group_col Character. Name of the grouping column (e.g. `"ancestry"`,
#'   `"feature"`).
#' @param metric Character. Name of the numeric column to compare.
#'   Default `"IC50"`.
#'
#' @return Invisibly, a named list with elements:
#'   - `anova` — an [stats::aov()] object
#'   - `tukey` — a [stats::TukeyHSD()] object
#' @export
#'
#' @examples
#' \dontrun{
#' # by ancestry
#' res <- runAOV(metrics, group_col = "ancestry", metric = "IC50")
#'
#' # by feature
#' res <- runAOV(metrics, group_col = "feature", metric = "AUC")
#' }
runAOV <- function(metrics_df, group_col, metric = "IC50") {
  validateCols(metrics_df, c(metric, group_col))
  
  formula <- stats::as.formula(paste(metric, "~", group_col))
  aov_res <- stats::aov(formula, data = metrics_df)
  tukey   <- stats::TukeyHSD(aov_res)
  
  cat(sprintf("\n--- ANOVA: %s ~ %s ---\n", metric, group_col))
  print(summary(aov_res))
  cat(sprintf("\n--- Tukey HSD post-hoc ---\n"))
  print(tukey)
  
  invisible(list(anova = aov_res, tukey = tukey))
}


#' Test sensitivity differences between two groups
#'
#' Runs a Wilcoxon rank-sum test comparing IC50 or AUC between two levels of
#' a grouping variable. Useful for mutant vs. wild-type or resistant vs.
#' sensitive comparisons.
#'
#' @param metrics_df A data frame from [summarizeDRC()] or [extractALL()].
#' @param group_col Character. Name of the binary grouping column.
#' @param metric Character. Name of the sensitivity column. Default `"IC50"`.
#'
#' @return A tidy data frame with columns `group_col`, `group_a`, `group_b`,
#'   `metric`, `statistic`, `p_value`, `median_a`, `median_b`.
#' @export
#'
#' @examples
#' \dontrun{
#' testGroups(metrics, group_col = "feature", metric = "AUC")
#' }
testGroups <- function(metrics_df, group_col, metric = "IC50") {
  validate_columns(metrics_df, c(group_col, metric))
  
  groups <- unique(metrics_df[[group_col]])
  if (length(groups) != 2) {
    stop(
      sprintf(
        "`%s` must have exactly 2 levels for a pairwise test; found: %s",
        group_col,
        paste(groups, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  vals_a <- metrics_df[[metric]][metrics_df[[group_col]] == groups[1]]
  vals_b <- metrics_df[[metric]][metrics_df[[group_col]] == groups[2]]
  
  wt <- stats::wilcox.test(vals_a, vals_b)
  
  data.frame(
    group_col  = group_col,
    group_a    = as.character(groups[1]),
    group_b    = as.character(groups[2]),
    metric     = metric,
    statistic  = wt$statistic,
    p_value    = wt$p.value,
    median_a   = stats::median(vals_a, na.rm = TRUE),
    median_b   = stats::median(vals_b, na.rm = TRUE),
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
}