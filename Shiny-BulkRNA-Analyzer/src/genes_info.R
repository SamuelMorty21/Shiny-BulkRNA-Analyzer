library(shiny)

gene_ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
      .blue-button {
        background-color: #007bff; /* 蓝色背景颜色 */
        color: white; /* 文本颜色为白色 */
      }
    ")
    )
  ),
  
  titlePanel("基因输入界面"),
  fluidRow(
    column(6, textInput("gene_input", "基因名：", placeholder = "输入基因名")),
    column(6, textInput("category_input", "基因类别：", placeholder = "输入基因类别"))
  ),
  actionButton("add_row", "添加基因"),
  actionButton("remove_row", "减少最后一行"),
  actionButton("store_data", "存储基因数据", class = "blue-button"),
  dataTableOutput("gene_data_table")
)

gene_server <- function(input, output, session) {
  gene_data <- reactiveValues(data = data.frame(GeneName = character(0), Category = character(0)))
  
  
  observeEvent(input$add_row, {
    gene_names <- unlist(strsplit(input$gene_input, "\\s+"))  # 使用空格作为分隔符
    gene_category <- isolate(input$category_input)
    
    for (gene_name in gene_names) {
      if (gene_name != "") {
        gene_data$data <- rbind(gene_data$data, data.frame(GeneName = gene_name, Category = gene_category))
      }
    }
  })
  
  observeEvent(input$remove_row, {
    gene_data$data <- gene_data$data[-nrow(gene_data$data), , drop = FALSE]
  })
  
  observeEvent(input$store_data, {
    # 在这里可以将基因数据保存到文件或者数据库中
    assign("genes_info", gene_data$data, envir = .GlobalEnv)
    saveRDS(gene_data$data, paste0(file_path, "genes_info.rds"))
    # 显示存储成功的通知弹窗
    showNotification("基因数据已成功存储！", type = "message")
  })
  
  output$gene_data_table <- renderDataTable({
    gene_data$data
  })
}


genes_info <- function(){
  shinyApp(gene_ui, gene_server)
}
