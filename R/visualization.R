#' Volcano plot for DEG results
#'
#' @importFrom magrittr %>%
#' @param df A data.frame with log2FoldChange, significance, expression columns.
#' @param title Character plot title.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal
#' @export
plot_volcano <- function(df, title) {
  ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = significance)) +
    ggplot2::geom_point(ggplot2::aes(color = expression), alpha = 0.6) +
    ggplot2::labs(
      title = title,
      x     = "Log2 Fold Change",
      y     = "-Log10(padj)"
    ) +
    ggplot2::theme_minimal()
}

#' PCA plot for a DESeq2 dataset
#'
#' @param dds A DESeqDataSet object.
#' @param intgroup Character vector of metadata columns for grouping.
#' @return A ggplot object.
#' @importFrom DESeq2 vst plotPCA
#' @export
plot_PCA <- function(dds, intgroup = c("group_time", "batch")) {
  vsd <- DESeq2::vst(dds, blind = FALSE)
  ggplot2::autoplot(vsd, intgroup = intgroup)
}
