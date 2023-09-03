perform_GOs <- function(log2fold, mode = 2){
  
  if (substr(file_path, nchar(file_path), nchar(file_path)) != "/") {
    # 如果最后一个字符不是斜杠，添加斜杠
    file_path <- paste0(file_path, "/")
  }
  
  group_info <- sample_info[sample_info$Select=="Yes",]
  fold = 2^log2fold
  path = paste0(file_path, "results/DEG/",
                paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1],"_all_DEGs.csv")
  )
  resSig <- read.csv(path,header = T,row.names = 1, check.names = FALSE)
  library(org.Hs.eg.db)
  ###### up DEGs #######
  up <- subset(resSig, log2FoldChange > log2fold)
  up$id <- mapIds(org.Hs.eg.db,keys=row.names(up),column="ENTREZID",
                  keytype="SYMBOL",multiVals="first")
  which(duplicated(up$id))
  id1 <- unique(na.omit(up$id)) 
  go.up <- enrichGO(id1, OrgDb = org.Hs.eg.db,
                    pvalueCutoff = 0.05, 
                    qvalueCutoff = 0.2,keyType = 'ENTREZID',
                    readable = T)
  filename1 = paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1],"_upDEGs_",fold,"fold")
  p1 <- dotplot(go.up, title = filename1)
  write.csv(go.up@result,paste0(paste0(file_path, "results/","GO_",filename1,".csv")))
  
  ################################################################################
  down <- subset(resSig, log2FoldChange < (-log2fold))
  down$id <- mapIds(org.Hs.eg.db,keys=row.names(down),column="ENTREZID",
                    keytype="SYMBOL",multiVals="first")
  which(duplicated(down$id))
  id2 <- unique(na.omit(down$id)) 
  go.down <- enrichGO(id2, OrgDb = org.Hs.eg.db,
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.2,keyType = 'ENTREZID',
                      readable = T)
  filename2 = paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1],"_downDEGs_",fold,"fold")
  p2 <- dotplot(go.down, title = filename2)
  write.csv(go.down@result,paste0(paste0(file_path, "results/","GO_",filename2,".csv")))
  
  if (mode == 2) {
    combined_plot <- plot_grid(p1, p2, ncol = 2)
    pdf(paste0(file_path, "plots/GO_", unique(group_info$Group)[2], "_vs_", unique(group_info$Group)[1], "_", fold, "fold.pdf"), width = 14, height = 7)
    print(combined_plot)
    dev.off()
  }
  
  if (mode == 1) {
    combined_plot <- plot_grid(p1, p2, ncol = 2)
    pdf(paste0(file_path, "plots/GO_", unique(group_info$Group)[2], "_vs_", unique(group_info$Group)[1], "_", fold, "fold.pdf"), width = 14, height = 7)
    print(combined_plot)
    dev.off()
  }
  
  
}

