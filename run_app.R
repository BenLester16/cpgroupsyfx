# run_app.R — 同时启动 Shiny (3838) 与 Plumber API (8000)

suppressPackageStartupMessages({
  library(callr)
  library(shiny)
  library(plumber)
})

# Shiny App 目录（相对路径）
app_dir <- "inst/shiny"
# Plumber 路由文件
plumber_file <- "R/plumber.R"

# 检查路径
if (!dir.exists(app_dir)) {
  stop("找不到 Shiny 目录：", app_dir, "\n请确认当前工作目录为包根目录。")
}
if (!file.exists(plumber_file)) {
  stop("找不到 Plumber 路由文件：", plumber_file)
}

# 1. 后台启动 Shiny
shiny_proc <- callr::r_bg(
  func = function(path) {
    shiny::runApp(path, host = "0.0.0.0", port = 3838)
  },
  args   = list(app_dir),
  stdout = "|", stderr = "|"
)
message("✅ Shiny 已在后台启动：http://localhost:3838")

# 2. 前台启动 Plumber API
pr <- plumber::plumb(plumber_file)
message("🔌 Plumber API 启动：http://localhost:8000")
pr$run(host = "0.0.0.0", port = 8000)
