setwd("D:/Shiny-BulkRNA-Analyzer/Shiny-BulkRNA-Analyzer/test")
roots <- getwd()
source(file.path(dirname(roots), "preload_dependencies.R"))
preload_dependencies(srcpath = roots)

# you can also redefine the file_path, sample info and group info by execute "load_sample_info()"

file_path <- roots
sample_info <- readRDS(file.path(roots, "sample_info.rds"))
group_info <- readRDS(file.path(roots, "group_info.rds"))

# You can load the gene info you are interested by execute gene_info()

genes_info <- readRDS(file.path(roots, "genes_info.rds"))

main <- function(x){
  perform_count_matrix_and_logtpm()
  perform_PCA()
  perform_DEGs(log2fold = 1) # 2fold
  perform_DEGs(log2fold = 2) # 4fold
  perform_DEGs(log2fold = 3)# 8 fold
  Heatmap_result <- perform_Heatmap()
  heatmap_fig <- create_custom_heatmap(
    data = Heatmap_result[[1]][, selected_cols = colnames(Heatmap_result[[1]])], # You can also select you interested cols by ht_col()
    col_data = Heatmap_result[[3]],
    genes_info = genes_info,
    show_column_names = TRUE,
    rect_gp = gpar(col = "white", lwd = 2)
    # rect_gp = gpar(col = NA)
  )
  pdf(file.path(file_path, "plots/Heatmap.pdf"),width = Heatmap_result[[4]], height=Heatmap_result[[5]])
  draw(heatmap_fig)
  dev.off()
  perform_GOs(
    log2fold = 1, mode = 2)
  perform_GOs(
    log2fold = 2, mode = 2)
}

main()