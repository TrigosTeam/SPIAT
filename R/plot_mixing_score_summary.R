#' plot_mixing_score_summary
#'
#' @description Plots mixing scores against target cells from the input mixing_score_summary data.frame.
#' It produces multiple scatterplots for each reference cell observed, with line of best fit. 
#' It gives a Rho statistic with p-value based off Spearman's correlation.
#' It produces different sets of plots based off mixing scores and normalised mixing scores.

#' @param mixing_score_summary SingleCellExperiment object in the form of the output of format_image_to_sce
#' @param individual.plots (OPTIONAL) Individual plots for each reference cell are given, rather than being arranged in a grid.
#' @import ggpubr
#' @import gridExtra
#' @return Various plots of log2(mixing score) vs log2(number of target cells) for each reference cell.
#' @examples
#' plot_mixing_score <- plot_mixing_score_summary(mixing_score_summary) 
#' @export

plot_mixing_score_summary <- function(mixing_score_summary_df, individual.plots=FALSE) 
{
  mixing.df <- data.frame(mixing_score_dataframe)
  mixing.df$log_ref_no <- log2(mixing.df$Number_of_reference_cells)
  mixing.df$log_tar_no <- log2(mixing.df$Number_of_target_cells)
  mixing.df$log_mix_score <- log2(mixing.df$Mixing_score)
  mixing.df$log_normalised_score <- log2(mixing.df$Normalised_mixing_score)
  if (isFALSE(individual.plots)) {
    p1 <- ggscatter(mixing.df,
                    x="log_mix_score", y="log_tar_no", 
                    add = "reg.line",
                    conf.int=TRUE,
                    add.params=list(color="blue", fill="lightgrey")) +  
      stat_cor(method = "spearman", label.y=3) + 
      facet_grid(~Reference, scales="free") +
      ggtitle("Mixing scores for all reference cells")
    p2 <- ggscatter(mixing.df,
                    x="log_normalised_score", y="log_tar_no", 
                    add = "reg.line",
                    conf.int=TRUE,
                    add.params=list(color="blue", fill="lightgrey")) +  
      stat_cor(method = "spearman", label.y=3) + 
      facet_grid(~Reference, scales="free") +
      ggtitle("Mixing scores for all reference cells")
    grid.arrange(p1, p2, nrow = 2)
  }
  else {
    for (reference in unique(mixing.df$Reference)) {
      reference.df <- mixing.df[mixing.df$Reference == reference,]
      p1 <- ggscatter(reference.df,
                      x="log_mix_score", y="log_tar_no") 
      p2 <- ggscatter(reference.df,
                      x="log_normalised_score", y="log_tar_no") 
      grid.arrange(p1, p2, nrow = 1, top=paste("Mixing scores for", reference, "as reference cell"))
    }
  }
}