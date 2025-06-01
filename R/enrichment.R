# R/enrichment.R

#' GO and KEGG enrichment analysis
#'
#' @importFrom magrittr %>%
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom enrichplot dotplot
#' @importFrom ggplot2 ggtitle
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @param genes Character vector of Entrez gene IDs.
#' @param OrgDb AnnotationDb object, e.g. org.Hs.eg.db.
#' @return List with GO and KEGG results and ggplot objects.
#' @export
enrich_analysis <- function(genes, OrgDb) {
  ego <- clusterProfiler::enrichGO(
    gene          = genes,
    OrgDb         = OrgDb,
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  ek <- clusterProfiler::enrichKEGG(
    gene          = genes,
    organism      = "hsa",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )
  plots <- list(
    go_dot   = enrichplot::dotplot(ego) + ggplot2::ggtitle("GO BP Enrichment")
    # 如果你还想做 KEGG 条形图，可用:
    # , kegg_bar = enrichplot::barplot(ek) + ggplot2::ggtitle("KEGG Pathway Enrichment")
  )
  return(list(GO = ego, KEGG = ek, plots = plots))
}
