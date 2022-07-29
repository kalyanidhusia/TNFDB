#' Make a p-value histogram plot
#'
#' Internal function for plotting p-value histograms of raw- and adjusted p-values
#'
#' @param data Per-contrast DE results to be plotted, as prepared by
#'   \code{\link{prep_plot_model_data}}.
#' @param contrast The contrast being plotted. Used for generating the plot title.
#'
#' @return A ggplot object
#'
#' @examples
#' # No examples yet
static_pval_histogram <- function(data, contrast) {
   output <- with(data, {
      ggplot(data) +
         geom_histogram(aes(x = P.Value), binwidth = 0.025, color = "black", na.rm = T) +
         xlim(c(-0.05,1.05)) +
         theme_bw() +
         xlab("raw P value") +
         theme(panel.border = element_rect(fill = NA, color = "grey30")) +
         ggplot(data) +
         geom_histogram(aes(x = adj.P.Val), binwidth = 0.025, color = "black", na.rm = T) +
         xlim(c(-0.05,1.05)) +
         theme_bw() +
         xlab("adjusted P value") +
         theme(panel.border = element_rect(fill = NA, color = "grey30")) +
         patchwork::plot_annotation(
            title = stringr::str_replace(contrast, "_vs_", " vs ")
         )
   })
   output
}
