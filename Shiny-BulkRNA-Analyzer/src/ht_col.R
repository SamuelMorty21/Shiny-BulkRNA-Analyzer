

# 定义UI界面
ui <- fluidPage(
  titlePanel("选择HT的列名"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput("columns", "选择列名:", choices = NULL),
      actionButton("submitBtn", "确定")
    ),
    mainPanel(
      verbatimTextOutput("selectedColsText")
    )
  )
)

# 定义服务器逻辑
server <- function(input, output, session) {
  
  observe({
    # 检查Heatmap_result是否存在
    if (!is.null(Heatmap_result)) {
      # 更新checkboxGroupInput的选项
      updateCheckboxGroupInput(session, "columns", choices = colnames(Heatmap_result[[1]]))
    }
  })
  
  # 创建一个空的reactiveValues对象来存储选择的列名
  selected_cols <- reactiveValues(value = NULL)
  
  # 当确定按钮被点击时，更新选择的列名
  observeEvent(input$submitBtn, {
    selected_cols$value <- input$columns
    showModal(
      modalDialog(
        title = "提示",
        "选好了！",
        footer = actionButton("closeModal", "关闭")
      )
    )
  })
  
  # 关闭模态框
  observeEvent(input$closeModal, {
    removeModal()
  })
  
  # 在UI中显示已选列名
  output$selectedColsText <- renderPrint({
    if (!is.null(selected_cols$value)) {
      paste("已选列名:", paste(selected_cols$value, collapse = ", "))
    }
  })
  
  # 将结果分配给selected_cols变量
  observe({
    if (!is.null(selected_cols$value)) {
      assign("selected_cols", selected_cols$value, envir = .GlobalEnv)
    }
  })
}

# 运行Shiny应用
ht_getcol <- function(){
  shinyApp(ui, server)
}
