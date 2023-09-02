perform_DEGs <- function(log2fold = 2, padjvalue = 0.05){
  deg.count <- read.csv(paste0(file_path,"count_matrix/count_matrix.csv"),header = TRUE, row.names = 1, check.names = FALSE)[,-1]
  
  fold = 2^log2fold
  # sample_names <<- colnames(deg.count) # samples_names传入
  # load_sample_selection()
  deg.count <- deg.count[,sample_info$Select=="Yes"]
  deg.count <- deg.count[apply(deg.count,1,function(x) sum(x>1) > 0),]  %>% as.matrix()
  group_info <- sample_info[sample_info$Select=="Yes",]
  info <- data.frame(row.names = colnames(deg.count), group = group_info$Group)
  info$group <- factor(info$group, levels = unique(info$group)) 
  library(DESeq2)
  dds1 <- DESeqDataSetFromMatrix(deg.count,colData = info, design = ~group)
  dds1 <- DESeq(dds1)
  resultsNames(dds1)
  res1 <- results(dds1)
  resOrdered1 <- res1[order(res1$log2FoldChange),] %>% as.data.frame()
  resSig <- subset(resOrdered1, padj < padjvalue)
  dir.create(paste0(file_path,"results/DEG/"))
  #
  # dir.create(paste0(file_path,"results/DEG/",unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1]))
  filename1 = paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1],"_all_DEGs",".csv")
  filename2 = paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1],"_upDEGs_",fold,"fold",".csv")
  filename3 = paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1],"_downDEGs_",fold,"fold",".csv")
  filename4 = paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1],"_DEGs_",fold,"fold",".csv")
  write.csv(resSig,paste0(file_path, "results/DEG/",filename1))
  write.csv(subset(resSig, log2FoldChange > log2fold ),paste0(file_path, "results/DEG/",filename2))
  write.csv(subset(resSig, log2FoldChange < -log2fold ),paste0(file_path, "results/DEG/",filename3))
  write.csv(rbind(subset(resSig, log2FoldChange > log2fold ),subset(resSig, log2FoldChange < -log2fold )),paste0(file_path, "results/DEG/",filename4))
}