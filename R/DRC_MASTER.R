#dose response pipeline
library(drc) #citation Ritz, C., Baty, F., Streibig, J. C., Gerhard, D. (2015) Dose-Response Analysis Using R PLOS ONE, 10(12), e0146021
library(tidyverse)
library(RColorBrewer)
library(ggrepel)

#data preprocessing

load_data <- function(filepath) {
  data <- read.csv(filepath)
  data$dose <- as.numeric(data$dose)
  data$response[data$response < 0] <- 0 #set neg val to 0
  data$cell_line <- factor(data$cell_line)
  data$ancestry <- factor(data$ancestry)
  data$translocation <- factor(data$translocation)
  return(data)
}

#combine reps by cell line and dose
combine_reps <- function(data) {
  data_avg <- data|>
    group_by(cell_line, dose, ancestry, translocation)|>
    summarise(sd = sd(response, na.rm = TRUE),
              mean_response = mean(response,na.rm = TRUE),
              n = n(),
              .groups = "drop")|>
    mutate(cell_line = factor(cell_line))
  return(data_avg)
}

#now fit the models
fit_cellline <- function(data, cell_line_name, model = LL.4()) {
  drm(mean_response ~ dose,
                 data = subset(data,cell_line == cell_line_name),
                 fct = model)
}

#fit all cell line models
fit_all <- function(data, cell_lines) {
  setNames(
    map(cell_lines, ~fit_cellline(data, .x)),
    cell_lines
  )
  
}

#extract IC50/ED50
extract_ic50 <- function(fit, cell_line, ancestry = NA, translocation = NA) {
  ic50 <- suppressWarnings(ED(fit, 50, interval = "delta"))
  ic50_df <- as.data.frame(ic50)
  colnames(ic50_df)[1:4] <- c("Estimate", "SE", "Lower", "Upper")
  ic50_df$cell_line <- cell_line
  ic50_df$ancestry <- ancestry
  ic50_df$translocation <- translocation
  ic50_df
}


#extract all ic50's
extract_all <- function(fits, agg_data) {
  cell_line_info <- agg_data|>
    dplyr::select(cell_line, ancestry, translocation)|>
    dplyr::distinct()
  
  purrr::map_df(names(fits), function(cl) {
    info <- dplyr::filter(cell_line_info, cell_line == cl)
    extract_ic50(fits[[cl]], cl, info$ancestry[1], info$translocation[1])
  })
}

pred <- function(fits, agg_data, n_points = 200) {
  #want my metadata to end up in my final plot instead of getting lost
  meta <- agg_data |>
    dplyr::select(cell_line, ancestry, translocation)|>
    dplyr::distinct()
  dose_seq <- seq(min(agg_data$dose), max(agg_data$dose), length.out = n_points)
  
  tidyr::expand_grid(dose = dose_seq, cell_line = names(fits))|>
    dplyr::mutate(
      predicted_response = purrr::map2_dbl(
        cell_line, dose, 
                                         ~as.numeric(predict(fits[[.x]], newdata = data.frame(dose = .y))
                                         ))
      )|>
    dplyr::left_join(meta, by ="cell_line")
}


#' Plot dose-response curves
#' @param agg_data Aggregated data from combine_reps()
#' @param pred_data Predictions from pred()
#' @param title Plot title
#' @export
plot_drc <- function(agg_data, pred_data, title = "", colors = NULL) {
  p <- ggplot(agg_data, aes(x = dose, y = mean_response, color = cell_line)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_response - sd, ymax = mean_response + sd),
                  width = 0.1, alpha = 0.6) +
    geom_line(data = pred_data,
              aes(x = dose, y = predicted_response, color = cell_line),
              size = 1) +
    scale_x_log10(name = "Dose [µM]") +
    scale_y_continuous(name = "Cell Death (%)", limits = c(0, 110)) +
    theme_classic() +
    theme(legend.title = element_blank(),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x =  element_text(size = 12),
          axis.text.y =  element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 16, face = "bold")) +
    ggtitle(title)
  
  if (!is.null(colors)) {
    p <- p + colors
  }
  return(p)
}

#' IC50 bar plot
#' @param ic50_df IC50 data frame from extract_all()
#' @param title Plot title
#' @export
plot_ic50 <- function(ic50_df, title = "", colors = NULL) {
  p <- ggplot(ic50_df, aes(x = reorder(cell_line, Estimate), y = Estimate, fill = cell_line, linetype = ancestry)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
    facet_wrap(~translocation) +
    scale_y_continuous(name = "IC50 [µM]") +
    scale_x_discrete(name = "") +
    scale_linetype_manual(values = c("African" = "solid",
                                     "European" = "dashed",
                                     "East Asian" = "dotted")) +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "right",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 16, face = "bold")) +
    ggtitle(title)
  
  if(!is.null(colors)) {
    p <- p + colors
  }
  return(p)
}

plot_drc_by_ancestry <- function(agg_data, pred_data, title = "", colors = NULL) {
  #to get a plot with ancestry
  #pred_data <- pred_data|>
   # left_join(agg_data|>
    #            dplyr::select(cell_line, ancestry, translocation)|>
     #           dplyr::distinct(),
      #        by = "cell_line")
  label_data <- pred_data |>
    dplyr::group_by(cell_line, translocation) |>
    dplyr::slice_max(dose, n = 1) |>
    dplyr::ungroup()
  
  p <- ggplot(agg_data, aes(x = dose, y = mean_response, 
                            color = ancestry, linetype = ancestry, group = cell_line)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_response - sd, ymax = mean_response + sd),
                  width = 0.1, alpha = 0.6) +
    geom_line(data = pred_data,
              aes(x = dose, y = predicted_response, 
                  color = ancestry, linetype = ancestry),
              size = 1) +
    facet_wrap(~translocation, scales = "free") +
    scale_x_log10(name = "Dose [µM]") +
    scale_y_continuous(name = "Effect (%)", limits = c(0, 110)) +
    scale_linetype_manual(values = c("African" = "solid",
                                     "European" = "dashed",
                                     "East Asian" = "dotted")) +
    scale_color_manual(values = c("African" = "#EE4C97FF",
                                  "European" = "#7876B1FF",
                                  "East Asian" = "#6F99ADFF")) +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(x = dose,
          y = predicted_response,
          label = cell_line,
          color = ancestry),
      size = 3,
      fontface = "bold",
      show.legend = FALSE,
      box.padding = 0.4,
      max.overlaps = Inf,
      nudge_y = 5
    ) +
    theme_classic() +
    theme(legend.title = element_blank(),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.x =  element_text(size = 14),
          axis.text.y =  element_text(size = 14),
          legend.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold")
          ) +
    ggtitle(title)
  
  if (!is.null(colors)) {
    p <- p + colors
  }
  return(p)
}

#pal_nejm()(8)
#[1] "#BC3C29FF" "#0072B5FF" "#E18727FF" "#20854EFF" "#7876B1FF" "#6F99ADFF" "#FFDC91FF"
#[8] "#EE4C97FF"

#boxplot
plot_box <- function(ic50_df, title = "", colors = NULL) {
  p <- ggplot(ic50_df, aes(x = reorder(cell_line, Estimate), y = Estimate, fill = cell_line)) +
    geom_boxplot(alpha = 0.8, outlier.size = 2) +
    geom_point(size = 3, alpha = 0.6, position = position_jitter(width = 0.1)) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, alpha = 0.6) +
    scale_y_continuous(name = "IC50 [µM]") +
    scale_x_discrete(name = "") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title)
  
  if(!is.null(colors)) {
    p <- p + colors
  }
  return(p)
}

#violin
plot_violin <- function(ic50_df, title = "", colors = NULL) {
  p <- ggplot(ic50_df, aes(x = reorder(cell_line, Estimate), y = Estimate, fill = cell_line)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.8, fill = "white") +
    geom_point(size = 3, alpha = 0.6, position = position_jitter(width = 0.1)) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, alpha = 0.6) +
    scale_y_continuous(name = "IC50 [µM]") +
    scale_x_discrete(name = "") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title)
  
  if(!is.null(colors)) {
    p <- p + colors
  }
  return(p)
}

plot_box_expand <- function(ic50_df, title = "", colors = NULL) {
  p <- ggplot(ic50_df, aes(x = reorder(cell_line, Estimate), y = Estimate, 
                           fill = cell_line, linetype = ancestry)) +
    geom_boxplot(alpha = 0.8, outlier.size = 2) +
    geom_point(size = 3, alpha = 0.6, position = position_jitter(width = 0.1)) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, alpha = 0.6) +
    facet_wrap(~translocation) +
    scale_y_continuous(name = "IC50 [µM]") +
    scale_x_discrete(name = "") +
    scale_linetype_manual(values = c("African" = "solid",
                                     "European" = "dashed",
                                     "East Asian" = "dotted")) +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = "right",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 16, face = "bold")) +
    ggtitle(title)
  
  if(!is.null(colors)) {
    p <- p + colors
  }
  return(p)
}
