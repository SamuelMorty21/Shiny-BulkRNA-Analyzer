preload_dependencies <- function(srcpath){
  source_list <- list.files(paste0(dirname(srcpath),"/src"), full.names = TRUE, pattern = "\\.R$")
  
  for (source in source_list) {
    source(source)
  }
  
  required_packages <- c("dplyr","org.Hs.eg.db","limma","ggplot2","circlize","tibble","pheatmap","gridExtra","ggrepel",
                         "DESeq2", 'ComplexHeatmap',"ggplotify","clusterProfiler","cowplot")
  # 逐个加载包
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = F)
    }
    library(pkg, character.only = TRUE)
  }
  
  if (exists("file_path")){
    if (file.exists(paste0(file_path, "genes_info.rds"))){
      genes_info <- readRDS(paste0(file_path, "genes_info.rds"))
    }
    
  }
}

