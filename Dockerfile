FROM rocker/r-ver:4.2.2

# 安装运行时依赖
RUN R -e "install.packages(c(
    'shiny','shinythemes','DESeq2','dplyr','tidyr','tibble',
    'ggplot2','clusterProfiler','org.Hs.eg.db','enrichplot',
    'plumber','callr'
  ), repos='https://cloud.r-project.org')"

# 复制包源码与运行脚本
COPY . /srv/RNAseqSimPkg
WORKDIR /srv/RNAseqSimPkg

# 安装包
RUN R -e "devtools::install()"

# Expose Shiny & Plumber ports
EXPOSE 3838 8000

# 运行 Shiny 与 Plumber
CMD ["Rscript", "run_app.R"]
