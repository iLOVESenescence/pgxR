#' Synthetic paclitaxel dose-response dataset
#'
#' A synthetic dataset simulating paclitaxel dose-response measurements
#' across 8 lung cancer cell lines with varying genomic backgrounds and
#' ancestry groups. Generated using 4-parameter log-logistic (LL.4) curves
#' with added gaussian noise, reflecting known paclitaxel biology:
#' NF1-deleted lines are sensitized, KRAS-mutant lines are resistant.
#'
#' @format A data frame with 192 rows and 5 columns:
#' \describe{
#'   \item{cell_line}{Cell line identifier (pLC1-pLC8)}
#'   \item{dose}{Paclitaxel concentration in nM}
#'   \item{response}{Percent cell death (0-100)}
#'   \item{ancestry}{1KGP/HGDP superpopulation abbreviation
#'     (AFR, EUR, EAS, AMR)}
#'   \item{feature}{Genomic feature group (NF1-del, KRAS-mut,
#'     TP53-mut, WT)}
#' }
#' @source Synthetically generated via data-raw/pgxr_example.R
#' @references
#'   Tomas M, et al. (2015). NF1 inactivation and paclitaxel sensitivity.
#'   KRAS mutations and taxane resistance in NSCLC.
"pgxr_example"