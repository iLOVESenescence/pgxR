# 
# theme.R — Shared ggplot2 theme and package-level aesthetic constants
# 

#' Default ancestry color palette
#'
#' A named character vector mapping HGDP/1KG superpopulation abbrevs to hex colors.
#' Pass a custom named vector to any plot function's
#' `ancestry_colors` argument to override.
#'
#'Superpopulastions follow the HGDP + 1KGP callset standard:
#' \itemize{
#' \item AFR - African
#' \item AMR - Admixed American
#' \item EAS - East Asian
#' \item EUR - European
#' \item SAS - South Asian 
#' }
#' 
#' @export
ANCESTRY_COLORS <- c(
  "AFR" = "#EE4C97FF",
  "AMR" = "#E18727FF",
  "EAS" = "#6F99ADFF",
  "EUR" = "#7876B1FF",
  "SAS" = "#20854EFF"
)

#' Default ancestry linetype mapping
#'
#' A named character vector mapping HGDP+1KGP labels to ggplot2 linetype strings. 
#' Pass a custom named vector to any plot function's
#' `ancestry_linetypes` argument to override.
#'
#' @export
ANCESTRY_LINETYPES <- c(
  "AFR" = "solid",
  "AMR" = "dotdash",
  "EAS" = "dotted",
  "EUR" = "dashed",
  "SAS" = "longdash"
)

#' Shared ggplot2 theme for pgxR plots
#'
#'
#' @param base_size Numeric. Base font size in points. Default `12`.
#' @param axis_title_size Numeric. Axis title size. Default `base_size + 2`.
#' @param bold_axis_titles Logical. Whether axis titles are bold. Default
#'   `FALSE`.
#'
#' @return A ggplot2 `theme` object.
#' @export
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   drc_theme()
drc_theme <- function(base_size        = 12,
                      axis_title_size  = base_size + 2,
                      bold_axis_titles = FALSE) {
  title_face <- if (bold_axis_titles) "bold" else "plain"
  
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.title  = ggplot2::element_blank(),
      legend.text   = ggplot2::element_text(size = base_size),
      axis.title.x  = ggplot2::element_text(size = axis_title_size,
                                            face = title_face),
      axis.title.y  = ggplot2::element_text(size = axis_title_size,
                                            face = title_face),
      axis.text.x   = ggplot2::element_text(size = base_size),
      axis.text.y   = ggplot2::element_text(size = base_size),
      strip.text    = ggplot2::element_text(size = base_size),
      plot.title    = ggplot2::element_text(size = base_size + 4,
                                            face = "bold")
    )
}

#' Map ancestry labels to HGDP+1KP super pops
#' 
#' Helpfer function for users whose data uses full pop names rather than standard abbrevs. 
#' Unrecognized labels are passed through unchanged with a warning.
#' 
#' @param x Character vectors of ancestry label
#' 
#' @return Character vector with standardized abbreviations
#' @export
#' 
#' @examples
#' standardize_ancestry(c("African", "European", "East Asian"))
 standardize_ancestry <- function(x) {
   lookup <- c(
     "African" = "AFR",
     "Admixed American" = "AMR",
     "Latino" = "AMR",
     "East Asian" = "EAS",
     "European" = "EUR",
     "South Asian" = "SAS"
 )

  result <- lookup[x]
  unrecognized <- is.na(result)
  
  if(any(unrecognized)) {
    warning(
      sprintf(
        "Unrecognized ancestry label(s) passed through unchanged: %s",
        paste(unique(x[unrecognized]), collapse = ",")
      ),
      call. = FALSE
    )
    result[unrecognized] <- x[unrecognized]
  }
  
  unname(result)
  }



#' Validate required columns in a data frame
#'
#' @param data A data frame.
#' @param required Character vector of required column names.
#' @param arg_name Name of the argument being checked, used in error messages.
#' @keywords internal
 # Internal helpers (not exported)
validate_columns <- function(data, required, arg_name = "data") {
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "`%s` is missing required column(s): %s",
        arg_name,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

utils::globalVariables(c(
  "cell_line", "dose", "response", "ancestry", "feature",
  "mean_response", "sd", "predicted_response",
  "Estimate", "Lower", "Upper", "AUC", "hill_slope",
  ".data", "reorder"
))