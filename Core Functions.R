# ---------------------------------------------
# 3. Core Functions (R/ directory)
# ---------------------------------------------
## R/simulation.R
#' Simulate sample metadata for two organisms
#'
#' @param org Character: organism label ("EC" or "PP").
#' @param times Character vector of time points.
#' @param groups Character vector of group labels.
#' @param batches Character vector of batch labels.
#' @param reps Integer: number of replicates per condition.
#' @return A data.frame of sample metadata.
#' @export
simulate_samples <- function(org,
                             times   = c("0h", "2h", "4h", "8h"),
                             groups  = if (org == "EC") c("less", "equal", "more", "all") else c("less", "equal", "more", "none"),
                             batches = paste0("B", 1:2),
                             reps    = 4) {
  expand.grid(
    organism = org,
    time     = times,
    group    = groups,
    batch    = batches,
    rep      = seq_len(reps),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(batch, time, group, rep) %>%
    dplyr::mutate(
      sample_id = paste(organism, group, time, batch, rep, sep = "_"),
      ratio0    = group
    )
}

#' Generate count data with negative binomial model
#'
#' @param meta Data.frame of sample metadata.
#' @param genes Character vector of gene IDs.
#' @param base_means Numeric named vector of base expression means per gene.
#' @param dispersions Numeric named vector of dispersions per gene.
#' @param size_factors Numeric named vector of library size factors per sample.
#' @return Tibble of counts for each gene-sample.
#' @export
generate_counts <- function(meta, genes, base_means, dispersions, size_factors) {
  tidyr::expand_grid(meta, gene = genes) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      mu0       = base_means[gene],
      fc        = dplyr::case_when(
        group != ifelse(organism == "EC", "all", "none") & runif(1) < 0.25 ~ sample(c(0.5, 1.5, 2, 0.25), 1),
        TRUE ~ 1
      ),
      batch_eff = stats::rlnorm(1, meanlog = 0, sdlog = 0.1),
      sf        = size_factors[sample_id],
      mu        = mu0 * fc * batch_eff * sf,
      count     = stats::rnbinom(1, mu = mu, size = 1 / dispersions[gene])
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample_id, organism, time, group, batch, gene, count)
}

## R/deseq_utils.R
#' Build DESeq2 dataset with batch correction
#'
#' @param counts_df Data.frame of counts with columns sample_id, organism, time, group, batch, gene, count.
#' @param org Character: organism subset.
#' @param size_factors Numeric named vector of size factors.
#' @return DESeqDataSet object.
#' @export
build_deseq_objects <- function(counts_df, org, size_factors) {
  df   <- dplyr::filter(counts_df, organism == org)
  ctrl <- ifelse(org == "EC", "all", "none")
  df   <- dplyr::filter(df, group != ctrl)
  
  count_mat <- df %>%
    dplyr::select(gene, sample_id, count) %>%
    tidyr::pivot_wider(names_from = sample_id, values_from = count) %>%
    tibble::column_to_rownames("gene")
  
  col_meta <- df %>%
    dplyr::distinct(sample_id, group, time, batch) %>%
    dplyr::mutate(
      group_time = factor(paste(group, time, sep = "_")),
      batch      = factor(batch)
    ) %>%
    dplyr::arrange(sample_id)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = col_meta,
    design    = ~ batch + group_time
  )
  DESeq2::sizeFactors(dds) <- size_factors[colnames(dds)]
  return(dds)
}

#' Extract differential expression results
#'
#' @param dds DESeqDataSet object after running DESeq().
#' @param contrasts List of contrast vectors for DESeq2 results().
#' @param lfcThresh Numeric threshold for |log2FC|.
#' @param padj Numeric adjusted p-value cutoff.
#' @return Named list of DEG data.frames.
#' @export
get_DEGs <- function(dds, contrasts, lfcThresh = 1, padj = 0.05) {
  res_list <- lapply(contrasts, function(ct) {
    res <- DESeq2::results(dds, contrast = ct, lfcThreshold = lfcThresh, alpha = padj)
    as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::filter(!is.na(padj), padj < padj) %>%
      dplyr::mutate(
        expression   = ifelse(log2FoldChange > 0, "up", "down"),
        significance = -log10(padj)
      )
  })
  names(res_list) <- sapply(contrasts, function(ct) paste(ct, collapse = ":"))
  return(res_list)
}

## R/visualization.R
#' Volcano plot for DEG results
#'
#' @param df DEG data.frame with log2FoldChange and significance columns.
#' @param title Character plot title.
#' @return ggplot object.
#' @export
plot_volcano <- function(df, title) {
  ggplot2::ggplot(df, aes(x = log2FoldChange, y = significance)) +
    ggplot2::geom_point(aes(color = expression), alpha = 0.6) +
    ggplot2::labs(title = title, x = "Log2 Fold Change", y = "-Log10(padj)") +
    ggplot2::theme_minimal()
}

#' PCA plot for DESeq2 dataset
#'
#' @param dds DESeqDataSet object.
#' @param intgroup Character vector of metadata columns for grouping.
#' @return ggplot object.
#' @export
plot_PCA <- function(dds, intgroup = c("group_time", "batch")) {
  vsd <- DESeq2::vst(dds, blind = FALSE)
  p   <- DESeq2::plotPCA(vsd, intgroup = intgroup)
  return(p)
}

## R/enrichment.R
#' GO and KEGG enrichment analysis
#'
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
    go_dot   = enrichplot::dotplot(ego) + ggplot2::ggtitle("GO BP Enrichment"),
    kegg_bar = enrichplot::barplot(ek) + ggplot2::ggtitle("KEGG Pathway Enrichment")
  )
  return(list(GO = ego, KEGG = ek, plots = plots))
}

#' Export DEG results to CSV files
#'
#' @param deg_list Named list of DEG data.frames.
#' @param outdir Character output directory path.
#' @export
export_DEG <- function(deg_list, outdir = "DEG_Results") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  invisible(lapply(names(deg_list), function(nm) {
    utils::write.csv(
      deg_list[[nm]],
      file = file.path(outdir, paste0(nm, ".csv")),
      row.names = FALSE
    )
  }))
}

# ---------------------------------------------
## R/simulation.R
#' Simulate sample metadata
simulate_samples <- function(org, times, groups, batches, reps) { ... }

#' Generate count matrix
generate_counts <- function(meta, genes, base_means, dispersions, size_factors) { ... }

## R/deseq_utils.R
#' Build DESeq2 dataset with batch correction
build_deseq_objects <- function(counts_df, org, size_factors) { ... }

#' Extract DEGs
get_DEGs <- function(dds, contrasts, lfcThresh, padj) { ... }

## R/visualization.R
#' Volcano plot
plot_volcano <- function(df, title) { ... }

#' PCA plot
plot_PCA <- function(dds, intgroup) { ... }

## R/enrichment.R
#' GO/KEGG enrichment
enrich_analysis <- function(genes, OrgDb) { ... }

#' Export DEG to CSV
export_DEG <- function(deg_list, outdir) { ... }