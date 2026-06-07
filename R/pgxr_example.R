#' Synthetic paclitaxel dose-response dataset
#'
#' A synthetic dataset for illustrating pgxR functionality. Dose-response
#' parameters were chosen to reflect the range of paclitaxel sensitivity
#' observed in human lung cancer cell lines (Georgiadis et al., 1997).
#' Ancestry and genomic features included to highlight faceting parameters.
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
#'   Georgiadis MS, Russell EK, Gazdar AF, Johnson BE. Paclitaxel cytotoxicity
#'   against human lung cancer cell lines increases with prolonged exposure
#'   durations. \emph{Clin Cancer Res.} 1997;3(3):449-454. PMID: 9815704.
"pgxr_example"