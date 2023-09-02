
main <- function(x){
  file_path = "D:/R/test/" # 最后一定要加/，不然会报错
  files <<- list.files(path = paste0(file_path,"raw_material/"))
  source("D:/R/src/rna_seq_pipeline.R")
  source("D:/R/src/label_io.R")
  source("D:/R/src/sample_info_io_server.R")
  source("D:/R/src/sample_selection_server.R")
  load_label_info()
  perform_count_matrix_and_logtpm()
  load_sample_info()
  perform_PCA()
  load_sample_selection()
  load_sample_info()
  perform_DEGs()
  perform_Heatmap()
}
