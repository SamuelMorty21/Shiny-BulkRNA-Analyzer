
# 安装和加载必要的包
required_packages <- c("shiny", "dplyr","shinyjs","shinyFiles")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = FALSE)
  }
  library(pkg, character.only = TRUE)
}

file_path <<- c()

# Shiny应用程序UI部分
sample_info_ui <- fluidPage(
  navbarPage(
    "ShinyBulk",
    
    tabPanel("信息填写",
             titlePanel("选择文件夹并加载标签信息"),
             
             shinyDirButton("folder_input", "选择文件夹", "选择一个文件夹"),
             
             
             
             actionButton("load_button", "加载标签信息"),
             
             
             verbatimTextOutput("load_status"),
             
             verbatimTextOutput('rawInputValue'),
             verbatimTextOutput('filepaths'),
             sidebarLayout(
               sidebarPanel(
                 textInput("sample_name", "样本名:", ""),
                 textInput("group", "Group:", ""),
                 textInput("batch", "Batch信息:", ""),  # 新增Batch信息输入
                 radioButtons("radio_choice", "是否用于差异分析:",
                            choices = c("Yes", "No"),
                            selected = "No"),
                 actionButton("add_sample", "添加样本"),
                 actionButton("delete_last_row", "删除末行"),
                 actionButton("save_data", "保存数据"),
               ),
               mainPanel(
                 verbatimTextOutput("initial_print"),
                 dataTableOutput("output_table")
               )
             )
    )
    # tabPanel("选择样本",
    #          titlePanel("选择样本"),
    #          
    #          # 多选框
    #          checkboxGroupInput("sample_choice", "选择样本:",
    #                             choices = NULL,
    #                             selected = NULL),
    #          
    #          verbatimTextOutput("selected_samples"),
    #          
    #          actionButton(inputId = "confirm_button", label = "确认")
    # )
  )
)

# Shiny应用程序服务器端部分
sample_info_server <- function(input, output, session) {
  
  sample_data <- reactiveValues(samples = data.frame(Sample = character(0), Group = character(0), Batch = character(0), Select = character(0)))
  
  shared_data <- reactiveValues(files = NULL, file_path = NULL)
  
  values <- reactiveValues(folder_path = NULL)
  
  roots_tmp = c(home = sub("/$", "", normalizePath("/", winslash = "/", mustWork = FALSE)))
  
  
  # if (file.exists(file.path(roots, "sample_info.rds"))){
  #   sample_data$samples <- readRDS(file.path(roots, "sample_info.rds"))
  # }
  # 

  
  
  shinyDirChoose(input, "folder_input", roots = roots_tmp, filetypes=c('', 'txt'), session = session)
  
  observe({
    if (!is.na(input$folder_input['path'])){
      output$rawInputValue <- renderPrint(paste0("Your selected path is ",
                                                 roots_tmp, paste(unlist(input$folder_input['path']), collapse = "/"),
                                                 "/"))
      file_path <<- paste0(roots_tmp, paste0(paste(unlist(input$folder_input['path']), collapse = "/")), "/")
      if (file.exists(file.path(file_path, "sample_info.rds"))){
        sample_data$samples <- readRDS(file.path(file_path, "sample_info.rds"))
        output$output_table <- renderDataTable({
          sample_data$samples
        })
      }
    }
  })
  
  observeEvent(input$load_button, {
    shared_data$file_path <- paste0(roots, paste0(paste(unlist(input$folder_input['path']), collapse = "/")), "/")
    shared_data$files <- list.files(path = paste0(file_path,"raw_material/"))
    
    # pathm = paste0(file_path, "sample_info.rds")
    # if(file.exists(pathm)){
    #   sample_info <<- readRDS(pathm)
    #   good_data <- reactiveValues(samples = sample_info[,1:3])
    # }
    # sample_data$samples <- rbind(sample_data$samples, good_data$samples)
    # output$output_table <- renderDataTable({
    #   sample_data$samples
    # })
    
    updateCheckboxGroupInput(session, "sample_choice", choices = shared_data$files, selected = NULL)
    showModal(
      modalDialog(
        h4("标签信息加载完成！"),
        footer = modalButton("关闭")
      )
    )
  })
  ##### page 2 #####
  
    # 添加Batch列
  
  output$initial_print <- renderText({
    initial_text <- character(0)
    
    if (!is.null(shared_data$files)) {
      for (file in shared_data$files) {
        initial_text <- c(initial_text, paste("Printing file:", file))
        initial_text <- c(initial_text, "\n")  # 添加一个空行作为换行
      }
      
      return(initial_text)
    } else {
      initial_text <- "Please ensure that the file has been loaded on the first page."
    }
  })
  
  observeEvent(input$add_sample, {
    sample_name <- input$sample_name
    group <- input$group
    batch <- input$batch  # 获取输入的Batch信息
    
    if (sample_name != "" && group != "") {
      new_row <- data.frame(Sample = sample_name, Group = group, Batch = batch, Select = input$radio_choice)  # 添加Batch信息
      sample_data$samples <- rbind(sample_data$samples, new_row)
    }
  })
  
  observeEvent(input$delete_last_row, {
    if (nrow(sample_data$samples) > 0) {
      sample_data$samples <- sample_data$samples[-nrow(sample_data$samples), ]
    }
  })
  
  observeEvent(input$save_data, {
    if (ncol(sample_data$samples)==3){
      sample_data_out <- reactiveValues(samples = data.frame(sample_data$samples))
    }else{
      sample_data_out <- reactiveValues(samples = data.frame(sample_data$samples))
    }
     
    if (nrow(sample_data$samples) > 0) {
      assign("sample_info", cbind(sample_data_out$samples), envir = .GlobalEnv)
      saveRDS(sample_data_out$samples, paste0(file_path, "sample_info.rds"))
      assign("group_info", sample_data_out$samples[sample_data_out$samples$Select=="Yes",], envir = .GlobalEnv)
      saveRDS(sample_data_out$samples[sample_data_out$samples$Select=="Yes",], paste0(file_path, "group_info.rds"))
      updateCheckboxGroupInput(session, "sample_choice", choices = sample_data$samples$Sample, selected = NULL)
      showNotification("保存成功！", type = "message")
    }
    else{
      
    }
  })
  
  
  output$output_table <- renderDataTable({
    sample_data$samples
  })
  
  # #显示用户选择的样本
  # output$selected_samples <- renderText({
  #   selected_samples <- input$sample_choice
  #   if (length(selected_samples) > 0) {
  #     paste("你选择了：", paste(selected_samples, collapse = ", "))
  #   } else {
  #     "你还没有选择样本。"
  #   }
  # })
  # #确认后存入样本
  # observeEvent(input$confirm_button, {
  #   selected_samples <- input$sample_choice
  #   if (length(selected_samples) > 0) {
  #     if (ncol(sample_data$samples)==3){
  #       sample_data_out <- reactiveValues(samples = data.frame(sample_data$samples))
  #     }else{
  #       sample_data_out <- reactiveValues(samples = data.frame(sample_data$samples))
  #     }
  #     sample_data_out$samples$Select[sample_data_out$samples$Sample %in% selected_samples] <- "Yes"
  #     assign("sample_info", sample_data_out$samples, envir = .GlobalEnv)
  #     assign("group_info", sample_data_out$samples[sample_data_out$samples$Select=="Yes",], envir = .GlobalEnv)
  #     saveRDS(sample_data_out$samples, paste0(file_path, "sample_info.rds"))
  #     saveRDS(sample_data_out$samples[sample_data_out$samples$Select=="Yes",], paste0(file_path, "group_info.rds"))
  #     }
  # })
  
  
}

load_sample_info <- function(){
  shinyApp(sample_info_ui, sample_info_server)
}
load_sample_info()          

