# inst/shiny/app.R

library(shiny)
library(shinythemes)
library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(RNAseqSimPkg)

# ---- UI ----
ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(tags$style(HTML("
    .shiny-busy {
      position: fixed; top: 0; left: 0;
      width: 100%; height: 100%;
      background: rgba(0,0,0,0.3);
      z-index: 10000;
    }
    .shiny-busy-panel {
      position: absolute; top: 50%; left: 50%;
      transform: translate(-50%, -50%);
      background: #fff; padding: 20px; border-radius: 8px;
      font-size: 1.2em;
    }
  "))),
  titlePanel("RNAseqSimPkg Explorer"),
  sidebarLayout(
    sidebarPanel(
      fileInput("counts", "Upload Count CSV", accept = c(".csv","text/csv")),
      selectInput("contrast", "Select Contrast", choices = NULL),
      actionButton("run", "Run Analysis")
    ),
    mainPanel(
      conditionalPanel(
        condition = "$('html').hasClass('shiny-busy')",
        div(class="shiny-busy-panel","Processing...")
      ),
      tabsetPanel(
        tabPanel("Volcano",
                 plotOutput("volcanoPlot"),
                 downloadButton("downloadVolcanoPDF","Download PDF"),
                 downloadButton("downloadVolcanoPNG","Download PNG")
        ),
        tabPanel("PCA",
                 plotOutput("pcaPlot"),
                 downloadButton("downloadPCAPDF","Download PDF"),
                 downloadButton("downloadPCAPNG","Download PNG")
        ),
        tabPanel("Enrichment (GO & KEGG)",
                 plotOutput("goPlot"),
                 downloadButton("downloadGOPDF","Download GO PDF"),
                 downloadButton("downloadGOPNG","Download GO PNG"),
                 hr(),
                 plotOutput("keggPlot"),
                 downloadButton("downloadKEGGPDF","Download KEGG PDF"),
                 downloadButton("downloadKEGGPng","Download KEGG PNG")
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  # 1. Read count matrix
  counts_mat <- reactive({
    req(input$counts)
    df <- read.csv(input$counts$datapath,
                   row.names       = 1,
                   check.names     = FALSE,
                   stringsAsFactors = FALSE)
    validate(need(ncol(df) >= 2, "Please upload at least two samples"))
    df
  })

  # 2. Build sample metadata
  sample_meta <- reactive({
    req(counts_mat())
    df <- counts_mat()
    parts <- strsplit(colnames(df), "_")
    meta  <- do.call(rbind, lapply(parts, function(x){
      data.frame(
        sample_id = paste(x, collapse = "_"),
        group     = x[1],
        time      = x[2],
        batch     = x[3],
        rep       = x[4],
        stringsAsFactors = FALSE
      )
    }))
    rownames(meta) <- colnames(df)

    # <<< fully qualify select() and everything() >>>
    meta %>%
      dplyr::mutate(
        group_time = factor(paste(group, time, sep = "_")),
        batch      = factor(batch)
      ) %>%
      dplyr::select(sample_id, dplyr::everything())
  })

  # 3. Populate contrasts
  observeEvent(input$counts, {
    req(sample_meta())
    meta  <- sample_meta()
    glvls <- unique(meta$group)
    if (all(c("less","all") %in% glvls)) {
      base <- "all"; cmp <- "less"
    } else {
      base <- "none"; cmp <- "less"
    }
    times   <- unique(meta$time)
    choices <- setNames(
      paste0(cmp, "_", times, "_vs_", base, "_", times),
      paste0(cmp, " vs ", base, " @ ", times)
    )
    updateSelectInput(session, "contrast",
                      choices  = choices,
                      selected = choices[1])
  })

  # 4. Run DESeq2 analysis
  analysis <- eventReactive(input$run, {
    req(counts_mat(), sample_meta(), input$contrast)
    df   <- counts_mat()
    meta <- sample_meta()
    design <- if (nlevels(meta$batch) > 1) ~ batch + group_time else ~ group_time

    dds <- DESeqDataSetFromMatrix(countData = df,
                                  colData   = meta,
                                  design    = design)
    dds <- DESeq(dds, fitType = "local")

    parts    <- strsplit(input$contrast, "_vs_")[[1]]
    contrast <- c("group_time", parts[1], parts[2])
    res      <- results(dds, contrast = contrast,
                        alpha = 0.05, lfcThreshold = 0)
    as.data.frame(res) %>%
      rownames_to_column("gene") %>%
      dplyr::filter(!is.na(padj)) %>%
      dplyr::mutate(
        expression   = ifelse(log2FoldChange > 0, "up", "down"),
        significance = -log10(padj)
      ) -> deg_df

    list(dds = dds, deg = deg_df)
  })

  # 5. Volcano plot + downloads (unchanged) …
  output$volcanoPlot <- renderPlot({
    req(analysis())
    df <- analysis()$deg
    validate(need(nrow(df) > 0, "No DE genes"))
    RNAseqSimPkg::plot_volcano(df, input$contrast)
  })
  output$downloadVolcanoPDF <- downloadHandler(
    filename = function() paste0("volcano_", input$contrast, ".pdf"),
    content  = function(file) {
      df <- analysis()$deg
      ggsave(file,
             plot   = RNAseqSimPkg::plot_volcano(df, input$contrast),
             device = "pdf", width = 8, height = 6)
    }
  )
  output$downloadVolcanoPNG <- downloadHandler(
    filename = function() paste0("volcano_", input$contrast, ".png"),
    content  = function(file) {
      df <- analysis()$deg
      ggsave(file,
             plot   = RNAseqSimPkg::plot_volcano(df, input$contrast),
             device = "png", width = 8, height = 6)
    }
  )

  # 6. PCA + downloads (unchanged) …
  output$pcaPlot <- renderPlot({
    req(analysis())
    dds <- analysis()$dds
    validate(need(ncol(dds) > 1, "Insufficient samples"))
    vsd <- tryCatch(
      vst(dds, blind = FALSE),
      error = function(e) varianceStabilizingTransformation(dds, blind = FALSE)
    )
    plotPCA(vsd, intgroup = c("group_time", "batch"))
  })
  output$downloadPCAPDF <- downloadHandler(
    filename = function() paste0("PCA_", input$contrast, ".pdf"),
    content  = function(file) {
      vsd <- tryCatch(
        vst(analysis()$dds, blind = FALSE),
        error = function(e) varianceStabilizingTransformation(analysis()$dds, blind = FALSE)
      )
      ggsave(file,
             plot   = plotPCA(vsd, intgroup = c("group_time","batch")),
             device = "pdf", width = 8, height = 6)
    }
  )
  output$downloadPCAPNG <- downloadHandler(
    filename = function() paste0("PCA_", input$contrast, ".png"),
    content  = function(file) {
      vsd <- tryCatch(
        vst(analysis()$dds, blind = FALSE),
        error = function(e) varianceStabilizingTransformation(analysis()$dds, blind = FALSE)
      )
      ggsave(file,
             plot   = plotPCA(vsd, intgroup = c("group_time","batch")),
             device = "png", width = 8, height = 6)
    }
  )

  # 7. GO enrichment barplot + downloads
  output$goPlot <- renderPlot({
    req(analysis())
    up  <- analysis()$deg$gene[analysis()$deg$expression == "up"]
    validate(need(length(up) > 0, "No up-regulated genes"))
    ego <- enrichGO(gene          = up,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    df_go <- as.data.frame(ego)
    if (nrow(df_go) == 0) {
      ggplot() + theme_void() +
        annotate("text", x=0.5, y=0.5, label="No significant GO terms", size=6)
    } else {
      topn <- df_go %>%
        arrange(p.adjust) %>%
        slice_head(n = 10)
      ggplot(topn, aes(x = reorder(Description, Count), y = Count, fill = -log10(p.adjust))) +
        geom_col() + coord_flip() +
        labs(title="GO BP Enrichment", x=NULL, fill="-log10(p.adjust)") +
        theme_minimal()
    }
  })
  output$downloadGOPDF <- downloadHandler(
    filename = function() paste0("GO_", input$contrast, ".pdf"),
    content  = function(file) {
      # repeat code to build the plot
      up  <- analysis()$deg$gene[analysis()$deg$expression == "up"]
      ego <- enrichGO(gene=up, OrgDb=org.Hs.eg.db, ont="BP",
                      pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
      df_go <- as.data.frame(ego)
      if (nrow(df_go) == 0) {
        p <- ggplot() + theme_void() +
          annotate("text", x=0.5,y=0.5,label="No significant GO terms", size=6)
      } else {
        topn <- df_go %>% arrange(p.adjust) %>% slice_head(n=10)
        p <- ggplot(topn, aes(x=reorder(Description, Count), y=Count, fill=-log10(p.adjust))) +
          geom_col() + coord_flip() +
          labs(title="GO BP Enrichment", x=NULL, fill="-log10(p.adjust)") +
          theme_minimal()
      }
      ggsave(file, plot = p, device = "pdf", width = 8, height = 6)
    }
  )
  output$downloadGOPNG <- downloadHandler(
    filename = function() paste0("GO_", input$contrast, ".png"),
    content  = function(file) {
      up  <- analysis()$deg$gene[analysis()$deg$expression == "up"]
      ego <- enrichGO(gene=up, OrgDb=org.Hs.eg.db, ont="BP",
                      pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
      df_go <- as.data.frame(ego)
      if (nrow(df_go) == 0) {
        p <- ggplot() + theme_void() +
          annotate("text", x=0.5,y=0.5,label="No significant GO terms", size=6)
      } else {
        topn <- df_go %>% arrange(p.adjust) %>% slice_head(n=10)
        p <- ggplot(topn, aes(x=reorder(Description, Count), y=Count, fill=-log10(p.adjust))) +
          geom_col() + coord_flip() +
          labs(title="GO BP Enrichment", x=NULL, fill="-log10(p.adjust)") +
          theme_minimal()
      }
      ggsave(file, plot = p, device = "png", width = 8, height = 6)
    }
  )

  # 8. KEGG enrichment barplot + downloads
  output$keggPlot <- renderPlot({
    req(analysis())
    up  <- analysis()$deg$gene[analysis()$deg$expression == "up"]
    validate(need(length(up) > 0, "No up-regulated genes"))
    ekegg <- enrichKEGG(gene=up, organism="hsa",
                        pAdjustMethod="BH", qvalueCutoff=0.05)
    df_kegg <- as.data.frame(ekegg)
    if (nrow(df_kegg) == 0) {
      ggplot() + theme_void() +
        annotate("text", x=0.5,y=0.5,label="No significant KEGG pathways", size=6)
    } else {
      topk <- df_kegg %>% arrange(p.adjust) %>% slice_head(n=10)
      ggplot(topk, aes(x=reorder(Description, Count), y=Count, fill=-log10(p.adjust))) +
        geom_col() + coord_flip() +
        labs(title="KEGG Pathway Enrichment", x=NULL, fill="-log10(p.adjust)") +
        theme_minimal()
    }
  })
  output$downloadKEGGPDF <- downloadHandler(
    filename = function() paste0("KEGG_", input$contrast, ".pdf"),
    content  = function(file) {
      up  <- analysis()$deg$gene[analysis()$deg$expression == "up"]
      ekegg <- enrichKEGG(gene=up, organism="hsa",
                          pAdjustMethod="BH", qvalueCutoff=0.05)
      df_kegg <- as.data.frame(ekegg)
      if (nrow(df_kegg) == 0) {
        p <- ggplot() + theme_void() +
          annotate("text", x=0.5,y=0.5,label="No significant KEGG pathways", size=6)
      } else {
        topk <- df_kegg %>% arrange(p.adjust) %>% slice_head(n=10)
        p <- ggplot(topk, aes(x=reorder(Description, Count), y=Count, fill=-log10(p.adjust))) +
          geom_col() + coord_flip() +
          labs(title="KEGG Pathway Enrichment", x=NULL, fill="-log10(p.adjust)") +
          theme_minimal()
      }
      ggsave(file, plot = p, device = "pdf", width = 8, height = 6)
    }
  )
  output$downloadKEGGPng <- downloadHandler(
    filename = function() paste0("KEGG_", input$contrast, ".png"),
    content  = function(file) {
      up  <- analysis()$deg$gene[analysis()$deg$expression == "up"]
      ekegg <- enrichKEGG(gene=up, organism="hsa",
                          pAdjustMethod="BH", qvalueCutoff=0.05)
      df_kegg <- as.data.frame(ekegg)
      if (nrow(df_kegg) == 0) {
        p <- ggplot() + theme_void() +
          annotate("text", x=0.5,y=0.5,label="No significant KEGG pathways", size=6)
      } else {
        topk <- df_kegg %>% arrange(p.adjust) %>% slice_head(n=10)
        p <- ggplot(topk, aes(x=reorder(Description, Count), y=Count, fill=-log10(p.adjust))) +
          geom_col() + coord_flip() +
          labs(title="KEGG Pathway Enrichment", x=NULL, fill="-log10(p.adjust)") +
          theme_minimal()
      }
      ggsave(file, plot = p, device = "png", width = 8, height = 6)
    }
  )
}

# ---- Launch ----
shinyApp(ui, server)
