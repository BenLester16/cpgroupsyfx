# R/plumber.R
#* @apiTitle RNAseqSimPkg API
#* @apiDescription Simulate & analyze RNA-seq data via RESTful endpoints.

#* 模拟样本元信息
#* @param org      string  Organism code, e.g. "EC"
#* @param times    list    时间点向量
#* @param groups   list    组别向量
#* @param batches  list    批次向量
#* @param reps     int     重复次数
#* @post /simulate
function(org, times, groups, batches, reps) {
  meta <- simulate_samples(
    org      = org,
    times    = strsplit(times, ",")[[1]],
    groups   = strsplit(groups, ",")[[1]],
    batches  = strsplit(batches, ",")[[1]],
    reps     = as.integer(reps)
  )
  jsonlite::toJSON(meta, pretty = TRUE, auto_unbox = TRUE)
}

#* 差异表达分析
#* @param counts_csv string  原始计数矩阵（CSV文本，首列为 gene）
#* @param contrast   string 对比标签，如 "less_0h_vs_all_0h"
#* @post /analyze
function(counts_csv, contrast) {
  # 1. 读取计数矩阵
  df <- read.csv(text = counts_csv, row.names = 1, check.names = FALSE)
  # 2. 构造 colData
  col_meta <- do.call(rbind, strsplit(colnames(df), "_")) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("group","time","batch","rep")) %>%
    tibble::rownames_to_column("sample_id")
  col_meta <- col_meta %>%
    dplyr::mutate(
      group_time = factor(paste(group, time, sep = "_")),
      batch      = factor(batch)
    )
  rownames(col_meta) <- col_meta$sample_id
  
  # 3. DESeq2 分析
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = df,
    colData   = col_meta,
    design    = if (nlevels(col_meta$batch)>1) ~ batch + group_time else ~ group_time
  )
  dds <- DESeq2::DESeq(dds, fitType = "local")
  
  # 4. 提取结果
  parts <- strsplit(contrast, "_vs_")[[1]]
  res   <- DESeq2::results(dds,
                           contrast     = c("group_time", parts[1], parts[2]),
                           alpha        = 0.05,
                           lfcThreshold = 0)
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(
      expression   = ifelse(log2FoldChange > 0, "up", "down"),
      significance = -log10(padj)
    )
  jsonlite::toJSON(res_df, pretty = TRUE, auto_unbox = TRUE)
}
