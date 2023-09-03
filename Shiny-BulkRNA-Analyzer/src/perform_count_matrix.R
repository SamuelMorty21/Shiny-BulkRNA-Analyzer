perform_count_matrix_and_logtpm <- function(){
  # dir：featurecount的存放路径
  # label_names：最终输出的文件的样本标签
  
  ## install and load the dependencies
  
  
  ## load the featurecount and generate the count matrix
  
  # 获取父级目录路径并创建文件夹
  
  folder_names <- c("count_matrix", "results", "plots")
  for (folder_name in folder_names) {
    folder_path <- file.path(file_path, folder_name)  # 构建文件夹路径
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)  # 创建文件夹
    }
  }
  files <<- list.files(path = file.path(file_path,"raw_material/"))
  if (is.null(files)){
    stop("Can not find the material files! Please check it again!")
  }
  expr <- lapply(files, function(x){
    tmp <- read.table(file = file.path(paste0(file_path,"/raw_material/"),x), header = T, check.names = FALSE) [,c(1:7)]
    return(tmp)})
  df <- do.call(cbind,expr) 
  df <- na.omit(df) 
  df <- df[,c(1:6,seq(7,ncol(df),by=7))] 
  df <- column_to_rownames(df,var = "Geneid")
  count <- df[,5:length(colnames(df))]
  
  
  # load_label_info()
  
  colnames(count)[2:length(colnames(count))] <- sample_info$Sample
  count2tpm <- function(df){
    kb <- df$Length /1e3
    rpk <-  df[,2:length(colnames(df))]/kb
    tpm <- t(t(rpk)/colSums(rpk, na.rm = TRUE)*1e6) %>% as.data.frame() %>% na.omit()
    return(tpm)
    
  }
  tpm <- count2tpm(count)
  logtpm <- log2(tpm+1)
  ######## Cor plot ###########
  Corela <- cor(logtpm)
  width = 0.7*ncol(logtpm)+3.5
  height = 0.75*width
  gg_Corela <- as.ggplot(
    pheatmap(Corela,
             display_numbers = TRUE,
             number_format = "%.3f",
             number_color = "black",
             cluster_cols = TRUE,
             angle_col = "45",
             fontsize = 10,                # 设置字体大小
             fontsize_col = 10,             # 设置列标签字体大小
             fontsize_row = 10,             # 设置行标签字体大小
             cluster_rows = TRUE,
             border_color = "black",
             cellwidth = 45,         # 设置单元格宽度
             cellheight = 30
    )
  )
  ggsave(file.path(file_path,"plots/Correlationship.pdf"), plot = gg_Corela,width = width, height=height)
  
  # 写出count matrix和tpm的matrix
  write.csv(count,file.path(file_path,"count_matrix/count_matrix.csv"))
  write.csv(logtpm,file.path(file_path,"count_matrix/logtpm_matrix.csv"))
  write.csv(tpm,file.path(file_path,"count_matrix/tpm_matrix.csv"))
}