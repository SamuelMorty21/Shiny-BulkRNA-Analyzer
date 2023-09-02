perform_Heatmap <- function(genes=NULL){
  ########### Heatmap of m6a-related genes #########
  genes <- genes_info$GeneName
  ht = read.csv(paste0(file_path,"count_matrix/logtpm_matrix.csv"), header = TRUE, row.names = 1 , check.names = FALSE)
  print(colnames(ht))
  ht <- as.matrix(ht)
  ht <- ht[genes,]
  col_fun = colorRamp2(c(min(ht), (min(ht)+max(ht))/2, max(ht)), c("#38629E","#F3F4B7", "#BF4431"))
  color <- colorRamp2(c(min(ht), (min(ht)+max(ht))/2, max(ht)), c("#38629E","#F3F4B7", "#BF4431"))
  width = 0.75*length(colnames(ht))+1
  height = 0.15*nrow(ht)+1
  return(list(ht,col_fun,color,width,height))
  
  
}