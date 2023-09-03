
perform <- function(tpm_matrix, file_DEG, log2folds){
  
  
  library(dplyr)
  DEG_UP <- rownames(subset(file_DEG, log2FoldChange> log2folds))
  DEG_DOWN <- rownames(subset(file_DEG, log2FoldChange< (-log2folds)))
  
  col <- "DEG0.05"
  tpm_matrix[col] <- 'no diff'
  
  tpm_matrix[DEG_UP,col] <- "Up(0.05)"
  tpm_matrix[DEG_DOWN,col] <- "Down(0.05)"
  
  
  col0<-rep("#00000022",nrow(tpm_matrix))
  col0[tpm_matrix$"DEG0.05"=="Down(0.05)"] <- "#FF000090" #红色 降低
  col0[tpm_matrix$"DEG0.05"=="Up(0.05)"]<- "#0000FF90" #蓝色 升高
  pch0<-rep(21,nrow(tpm_matrix))
  
  deg_samples <-sample_info[sample_info$Select=="Yes",]
  for(i in unique(deg_samples$Group)){
    tmp <- (subset(deg_samples, Group==i))$Sample
    names <- paste0("mean_",i)
    tpm_matrix[names] <- rowSums(tpm_matrix[tmp])/2
  }
  table_deg <- as.data.frame(t(as.data.frame(table(tpm_matrix$DEG0.05))))
  
  colnames(table_deg) <- table_deg[1,]
  table_deg <- table_deg[-1,]
  # 创建一个包含所有可能列名的向量
  
  all_columns <- c("Down(0.05)", "no diff", "Up(0.05)")
  
  # 创建一个初始值为0的数据框
  result_data <- data.frame(matrix(0, ncol = length(all_columns), nrow = 1))
  colnames(result_data) <- all_columns
  
  # 匹配数据表的列名并填充数据框
  matching_columns <- intersect(colnames(table_deg), all_columns)
  result_data[, matching_columns] <- table_deg[, matching_columns]
  
  return(list(tpm_matrix,col0,pch0,result_data))
}



perform_dotplot <- function(){
  matrix <- read.csv(file.path(file_path, "count_matrix/logtpm_matrix.csv"), row.names = 1, check.names = F)
  
  filename1 = paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1],"_all_DEGs",".csv")
  
  DEG_result <- read.csv(file.path(file_path, "results/DEG/",filename1), row.names = 1, check.names = F)
  
  matrix1 <- perform(matrix, DEG_result, 1)
  matrix2 <- perform(matrix, DEG_result, 2)
  matrix3 <- perform(matrix, DEG_result, 3)
  
  group_info <- sample_info[sample_info$Select=="Yes",]
  
  pdf(paste0(file_path,"plots/",unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1],"_dotplot.pdf"), width = 13, height = 4)
  
  
  par(mfrow=c(1,3))
  plot(matrix1[[1]][,ncol(matrix1[[1]])-1],
       matrix1[[1]][,ncol(matrix1[[1]])],
       bg=matrix1[[2]],pch=matrix1[[3]],col=NA,cex=1,axes=FALSE,xlab="",ylab="")
  title(xlab=unique(group_info$Group)[2],ylab=unique(group_info$Group)[1],line= 2.3,cex.lab= 1.3)
  title(main=paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1]," (q_val<0.05)"),line= 0.5)
  box();axis(1,at=c(0,5,10));axis(2,at=c(0,5,10))
  legend("topleft", c(paste0("Up (n=",matrix1[[4]]$`Up(0.05)`,")"),
                      paste0("Down (n=",matrix1[[4]]$`Down(0.05)`,")")), pt.bg=c("#0000FF90","#FF000090"), pch=c(21,21), cex=0.8)
  abline(c(.2,1),lty=2,lwd=2)
  abline(c(-.2,1),lty=2,lwd=2)
  
  plot(matrix2[[1]][,ncol(matrix2[[1]])-1],
       matrix2[[1]][,ncol(matrix2[[1]])],
       bg=matrix2[[2]],pch=matrix2[[3]],col=NA,cex=1,axes=FALSE,xlab="",ylab="")
  title(xlab=unique(group_info$Group)[2],ylab=unique(group_info$Group)[1],line= 2.3,cex.lab= 1.3)
  title(main=paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1]," (q_val<0.05)"),line= 0.5)
  box();axis(1,at=c(0,5,10));axis(2,at=c(0,5,10))
  legend("topleft", c(paste0("Up (n=",matrix2[[4]]$`Up(0.05)`,")"),
                      paste0("Down (n=",matrix2[[4]]$`Down(0.05)`,")")), pt.bg=c("#0000FF90","#FF000090"), pch=c(21,21), cex=0.8)
  abline(c(.2,1),lty=2,lwd=2)
  abline(c(-.2,1),lty=2,lwd=2)
  
  plot(matrix3[[1]][,ncol(matrix3[[1]])-1],
       matrix3[[1]][,ncol(matrix3[[1]])],
       bg=matrix3[[2]],pch=matrix3[[3]],col=NA,cex=1,axes=FALSE,xlab="",ylab="")
  title(xlab=unique(group_info$Group)[2],ylab=unique(group_info$Group)[1],line= 2.3,cex.lab= 1.3)
  title(main=paste0(unique(group_info$Group)[2],"_vs_",unique(group_info$Group)[1]," (q_val<0.05)"),line= 0.5)
  box();axis(1,at=c(0,5,10));axis(2,at=c(0,5,10))
  legend("topleft", c(paste0("Up (n=",matrix3[[4]]$`Up(0.05)`,")"),
                      paste0("Down (n=",matrix3[[4]]$`Down(0.05)`,")")), pt.bg=c("#0000FF90","#FF000090"), pch=c(21,21), cex=0.8)
  abline(c(.2,1),lty=2,lwd=2)
  abline(c(-.2,1),lty=2,lwd=2)
  
  dev.off()
}
