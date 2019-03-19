source("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Functions/Mutations_First_Work_Functions.R")
source("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Functions/CNV_Functions.R")
library(plyr)

old_relgenes_csvname <- "December14 38 genes list.csv"
old_rel_genes <- read.csv(file = paste("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Gene_Lists/",old_relgenes_csvname,sep=""), colClasses = "character")
old_rel_genes <- names(old_rel_genes[1,])

relgenes_csvname <- "Project_Genes.csv"
rel_genes <- read.csv(file = paste("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Gene_Lists/",relgenes_csvname,sep=""), colClasses = "character")
rel_genes <- names(rel_genes[1,])


TCGA_data <- load_TCGA_data("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Initial_TCGA_Data/gdc_WESdata")
TCGA_filtered_data <- lapply(TCGA_data, function(x) x[x$Hugo_Symbol %in% rel_genes,c("Hugo_Symbol","Tumor_Sample_Barcode","IMPACT") ])
TCGA_impact_data <- remove_low_impact(TCGA_filtered_data)
TCGA_impact_data <- lapply(TCGA_impact_data, function(x) {x$Tumor_Sample_Barcode <- substr(x$Tumor_Sample_Barcode,1,16);x})
TCGA_impact_data <- lapply(TCGA_impact_data, function(x) unique(x[,c("Hugo_Symbol","Tumor_Sample_Barcode")]))
TCGA_impact_data <- lapply(TCGA_impact_data, function(x) {x$Hugo_Symbol <- as.character(x$Hugo_Symbol);x})
TCGA_impact_data <- TCGA_impact_data[order(names(TCGA_impact_data))]

CNVdata <- load_cancer_CNVdata("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/CNV_Data")
CNV_filtered_data <- lapply(CNVdata, function(x) x[rel_genes, ])
#Replace . with - to match the IDs from mutation data
CNV_filtered_data <- lapply(CNV_filtered_data, function(x) {colnames(x) <- gsub("\\.","-",colnames(x));x})
#Keep only the first 16 characters to allow for the matching
CNV_filtered_data <- lapply(CNV_filtered_data, function(x) {colnames(x) <- substr(colnames(x),1,16);x})

#FALAR COM A ANA, ESTOU A ASSUMIR QUE OS PACIENTES DE CNV QUE NÃO TÊM MUTAÇÃO NÃO TÊM MESMO, MAS PODEM SIMPLESMENTE
#NÃO TER SIDO SEQUENCIADOS? VOU TENTAR CHECKAR SEMPRE

length(TCGA_impact_data)
length(CNV_filtered_data)



colours <- c("blue","white","red","black","forestgreen","black")
breaks <- c(-2.1,-0.1,0.1,2.1,4.1,5.1,7.1)

setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Plots/CNV/Basic_Heatmaps")
pdf("Basic_Genomic_Alterations_Heatmaps_Orig_Genes.pdf")

for (m in 1:length(TCGA_impact_data)){
  thisMUTdata <- TCGA_impact_data[[m]]
  cancer_symbol <- names(TCGA_impact_data)[m]
  thisCNVdata <- CNV_filtered_data[[cancer_symbol]]
  for (j in 1:nrow(thisMUTdata)){
    gene <- thisMUTdata[j,"Hugo_Symbol"]
    patient <- thisMUTdata[j,"Tumor_Sample_Barcode"]
    if ((patient %in% colnames(thisCNVdata)) & gene %in% rownames(thisCNVdata)){
      thisCNVdata[gene,patient] <- thisCNVdata[gene,patient]+5
    }
  }
  heatmap(thisCNVdata,Rowv = NA,Colv=NA,labCol=F,col=colours,breaks=breaks,scale="none",main=cancer_symbol)
  CNV_filtered_data[[cancer_symbol]] <- thisCNVdata
}

dev.off()



#Little extra just to get the new list in the format as the previous ones

setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Initial_TCGA_Data")
genes <- read.table(file = "Just_Project_Genes.txt", sep = "\t",fill=TRUE)
genes <- as.character(genes[,1])
genes <- unique(genes[genes != ""])
genes_to_remove <- c("B3GNT10","C1GALT1C1L")
genes <- genes[!genes %in% genes_to_remove]
genes <- sort(mapvalues(genes,from=c("B3GLCT","GALNT17","HPF1","LARGE1","LARGE2","OGA","PUM3","RXYLT1"),
                        to=c("B3GALTL","WBSCR17","C4orf27","LARGE","GYLTL1B","MGEA5","KIAA0020","TMEM5")))
write.table(genes,file="Project_Genes.csv",quote=F,row.names=F,col.names=F,eol=",")

all_genes <- c()

for (g in 1:length(TCGA_data)){
  all_genes <- c(all_genes,as.character(TCGA_data[[g]]$Hugo_Symbol))
}

all_genes <- unique(all_genes)

genes[!genes %in% all_genes]

"B3GALTL" %in% all_genes
"MGEA5" %in% all_genes

genes

#Genes to replace: B3GLCT with B3GALTL. B3GNT10 removed (pseudogene, new, not in data). C1GALT1C1L removed.
#GALNT17 replaced with WBSCR17. HPF1 replaces with C4orf27. LARGE1 replaced with LARGE. LARGE2 replaced with GYLTL1B 
#OGA replaced with MGEA5. PUM3 replaced with KIAA0020. RXYLT1 replaced with TMEM5



####BARPLOTS############################################################################
########################################################################################
#############################################

###Barplot - How many genes altered in at least 50% of patients, per cancer################

CNVandMUT_datatest <- lapply(CNVandMUT_data, function(x) abs(x))

CNVandMUT_datatest <- lapply(CNVandMUT_datatest, function(x) {x[x == 2] <- 1;x})
CNVandMUT_datatest <- lapply(CNVandMUT_datatest, function(x) {x[x == -2] <- -1;x})

tt <- CNVandMUT_datatest[[1]]

length(CNVandMUT_datatest)

results <- data.frame(matrix(ncol=2,nrow=length(CNVandMUT_datatest)))
  
for (i in 1:length(CNVandMUT_datatest)){
  cdata <- CNVandMUT_datatest[[i]]
  results[i,1] <- names(CNVandMUT_datatest[i])
  results[i,2] <- sum((apply(cdata,1, function(x) sum(x))) > (ncol(cdata)/2))
}

rownames(results) <- results$X1

counts <- results$X2
names(counts) <- results$X1

names(counts) %in% col_tissue[,1]

counts <- counts[names(counts) %in% col_tissue[,1]]

col_tissue
colors <- col_tissue[1:31,2]

barplot(counts,las=2,ylab="number of genes mutated in over half npatients",col=colors)

?barplot

load("GlycoColor.RData")


###Barplot - How many patients with at least 1  genomic alteration in glycosilation genes, per cancer##

