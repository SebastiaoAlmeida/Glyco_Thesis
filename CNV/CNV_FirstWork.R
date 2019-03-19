source("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Functions/CNV_Functions.R")
library(gplots)


CNVdata <- load_cancer_CNVdata("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/CNV_Data")
old_relgenes_csvname <- "December14 38 genes list.csv"
old_rel_genes <- read.csv(file = paste("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Gene_Lists/",old_relgenes_csvname,sep=""), colClasses = "character")
old_rel_genes <- names(old_rel_genes[1,])

relgenes_csvname <- "Project_Genes.csv"
rel_genes <- read.csv(file = paste("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Gene_Lists/",relgenes_csvname,sep=""), colClasses = "character")
rel_genes <- names(rel_genes)
#Bandaid solution
rel_genes <- rel_genes[rel_genes %in% rownames(CNVdata[[1]])]

colours <- c("steelblue4","steelblue1","white","lightcoral","indianred")
breaks <- c(-2.1,-1.1,-0.1,0.1,1.1,2.1)

setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Plots/CNV/Basic_Heatmaps")

pdf("Basic_CNV_Heatmaps_ProjectGenes.pdf")
for (n in 1:length(CNVdata)){
  cancer_data <- CNVdata[[n]]
  cancer_data <- cancer_data[rel_genes,]
  cancer_symbol <- names(CNVdata)[n]
  heatmap(cancer_data,Rowv = NA,Colv=NA,labCol=F,col=colours,breaks=breaks,scale="none",main=cancer_symbol)
}
dev.off()

pdf("Basic_CNV_Heatmaps_OriginalGenes.pdf")
for (n in 1:length(CNVdata)){
  cancer_data <- CNVdata[[n]]
  cancer_data <- cancer_data[old_rel_genes,]
  cancer_symbol <- names(CNVdata)[n]
  heatmap(cancer_data,Rowv = NA,Colv=NA,labCol=F,col=colours,breaks=breaks,scale="none",main=cancer_symbol)
}
dev.off()
