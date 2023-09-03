perform_PCA <- function( mode="n"){
  if (!file.exists(file.path(file_path,"count_matrix/count_matrix.csv"))) {
    stop("count matrix不存在。函数无法继续执行。")
  }
  linear <- FALSE
  # load_sample_info()
  for (i in unique(sample_info$Group)){
    tmp <- sample_info[which(sample_info$Group==i),]
    if(length(unique(tmp$Batch))==1){
      linear <- TRUE
      break
    }
  }
  if(any(sample_info$Batch=="")||linear){
    # 不去batch effect的模式
    if(linear){
      print("由于数据线性，不去除批次")
    }
    print("您使用的是不去除批次模式")
    count <- read.csv(file.path(file_path,"count_matrix/count_matrix.csv"),header = TRUE, row.names = 1, check.names = FALSE)[,-1]
    info <- data.frame(row.names = sample_info$Sample, 
                       group = sample_info$Group,
                       sample = sample_info$Sample)
    info$group <- as.factor(info$group)
    count <- count[apply(count,1,function(x) sum(x>1) > 0),]  %>% as.matrix()
    ddsOrig <- DESeqDataSetFromMatrix(countData=count, colData=info, 
                                      design=
                                        ~group
                                      #~sample
                                      #~1
    )
    dds <- DESeq(ddsOrig)
    vstd <- vst(dds, blind=FALSE)#表达量标准化（大样本量用 vst?小样本量用 rlog?）
    #rld <- rlog(dds)
    ##############################
  }else if (!any(sample_info$Batch=="")&&!linear){
    print("您使用的是批次模式")
    count <- read.csv(file.path(file_path,"count_matrix/count_matrix.csv"),header = TRUE, row.names = 1, check.names = FALSE)[,-1]
    info <- data.frame(row.names = sample_info$Sample, 
                       celltype = sample_info$Group,
                       batch = sample_info$Batch, #假设有两个批次，每个批次各有三种 cell type
                       sample = sample_info$Sample
    ) 
    #将向量转换为factor:
    info$batch <- as.factor(info$batch)
    info$celltype <- as.factor(info$celltype)
    info$sample <- as.factor(info$sample)
    
    count <- count[apply(count,1,function(x) sum(x>1) > 0),]  %>% as.matrix()
    library(DESeq2)
    ddsOrig <- DESeqDataSetFromMatrix(countData=count, colData=info, 
                                      design= ~batch + celltype 
                                      #请解释这一步design参数的设置的奥妙
                                      #~sample
                                      #~1
    )
    dds <- DESeq(ddsOrig)
    vstd <- vst(dds, blind=FALSE)#表达量标准化（大样本量用 vst;小样本量用 rlog）
    #rld <- rlog(dds)
    pcaData <- plotPCA(vstd, intgroup=c('group'), returnData=F) #此时得到没有去除批次效应的结果
    
    ## Remove batch effect by limma ##
    mat <- assay(vstd) #取出被 DESeq2 标准化过的表达矩阵
    mm <- model.matrix(~celltype, colData(vstd))
    mat <- limma::removeBatchEffect(mat, batch=vstd$batch, design=mm)
    assay(vstd) <- mat #把去除过 batch effect的表达矩阵放回 vstd bject中，覆盖掉原本的。
  }else{
    stop("You had chosen the unrecognized mode")
  }
  ############################
  pcaData <- plotPCA(vstd, intgroup=c('group'), returnData=T)
  
  ##获取PCA分析中各PC百分比
  pca <- prcomp(t(assay(vstd)))
  standard_deviations <- pca$sdev
  # 计算方差解释比例
  variance_explained <- standard_deviations^2 / sum(standard_deviations^2)
  variance_explained <- round(variance_explained, digits = 4)
  loadings <- as.data.frame(pca$rotation)
  loadings<-loadings[order(loadings$PC1,decreasing = T),]
  ##画图
  xmin <- min(pcaData$PC1)
  xmax <- max(pcaData$PC1)
  xmin <- xmin - 0.3*(xmax - xmin)
  plot_pca <- ggplot(pcaData, aes(PC1, PC2, 
                                  shape=group,
                                  color=group )) + 
    geom_point(aes(color=group),size = 2)+
    geom_text_repel(aes(label=name)) +
    #xlim(xmin, NA) +
    theme_bw()+
    theme(panel.grid = element_blank())+
    labs(x=paste0("PC1:",variance_explained[1]*100,"%"),y=paste0("PC2:",variance_explained[2]*100,"%"))
  ggsave(file.path(file_path,"plots/PCA.pdf"), plot_pca, width = 8, height = 7, dpi = 300)
  write.csv(loadings,file.path(file_path,"results/PCA_data_all.csv"))
}