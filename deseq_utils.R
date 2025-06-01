#' Build a DESeqDataSet with batch correction
#'
#' @importFrom magrittr %>%
#' @param counts_df Data.frame with columns sample_id, organism, time, group, batch, gene, count.
#' @param org Character: organism subset ("EC" or "PP").
#' @param size_factors Named numeric vector of library size factors.
#' @return A DESeqDataSet object.
#' @importFrom dplyr filter distinct mutate arrange select
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom DESeq2 DESeqDataSetFromMatrix sizeFactors
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
#' @param dds A DESeqDataSet object after running DESeq().
#' @param contrasts List of contrast vectors for DESeq2::results().
#' @param lfcThresh Numeric threshold for |log2 fold change|.
#' @param padj Numeric adjusted p‑value cutoff.
#' @return Named list of DEG data.frames.
#' @importFrom DESeq2 results
#' @export
get_DEGs <- function(dds, contrasts, lfcThresh = 1, padj = 0.05) {
  res_list <- lapply(contrasts, function(ct) {
    res <- DESeq2::results(
      dds,
      contrast     = ct,
      lfcThreshold = lfcThresh,
      alpha        = padj
    )
    df <- as.data.frame(res)
    df <- tibble::rownames_to_column(df, "gene")
    df <- dplyr::filter(df, !is.na(padj), padj < padj)
    df <- dplyr::mutate(df,
                        expression   = ifelse(log2FoldChange > 0, "up", "down"),
                        significance = -log10(padj))
    df
  })
  names(res_list) <- sapply(contrasts, function(ct) paste(ct, collapse = ":"))
  return(res_list)
}
