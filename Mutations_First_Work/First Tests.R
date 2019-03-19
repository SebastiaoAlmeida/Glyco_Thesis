#September 19th - First tests####

setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Mutations_First_Work")

BRCA_Initial_Data <- read.table(file = "../../Datasets/Initial_TCGA_Data/gdc_WESdata/995c0111-d90b-4140-bee7-3845436c3b42/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz", 
                             comment.char = "", header = TRUE, sep = "\t", quote = "", skip = 5);
  

nrow(BRCA_Initial_Data)

BRCA_Initial_Data[0,]
View(head(BRCA_Initial_Data))

table(BRCA_Initial_Data$Variant_Type)
table(BRCA_Initial_Data$IMPACT)

#Filter out variants without impact
BRCA_Initial_Data <- BRCA_Initial_Data[BRCA_Initial_Data$Variant_Classification != "Silent",]

#How many variants on each gene
nvariants_allgenes <- table(BRCA_Initial_Data$Hugo_Symbol)
nvariants_allgenes <- sort(nvariants_allgenes, decreasing = T); nvariants_allgenes[1:20] #Top 20 most mutated genes

#How many patients with variants for each gene
table(BRCA_Initial_Data$Hugo_Symbol)
#Preciso de saber qual a coluna com os IDs das pessoas
#Tumour_Sample_Barcode

#Defining what are the relevant genes####

#Load the list
#List of ~60
#relevant_genes <- read.csv(file = "../selected genes.csv") #Import the csv file with the gene info
#List of ~20
rel_genes <- read.csv(file = "../../Datasets/Initial_TCGA_Data/October30 33 genes list.csv", colClasses = "character")

#Eventually, clean the list
#rel_genes <- relevant_genes[2:nrow(relevant_genes),] #Remove the first row, it contains no info
#rel_genes <- relevant_genes$official.Symbol #Retrieve a vector with all the gene symbols
#rel_genes <- rel_genes[rel_genes != ""] #remove empty elements
#rel_genes <- as.character(rel_genes) #Turn into a vector
rel_genes <- names(rel_genes[1,])

rel_genes %in% BRCA_Initial_Data$Hugo_Symbol

BRCA_Initial_Data[BRCA_Initial_Data$Entrez_Gene_Id  == "29071",]

for (gene in 1:length(rel_genes)){
  if (rel_genes[gene] %in% BRCA_Initial_Data$Hugo_Symbol == F){
    print(rel_genes[gene])
  }
}


#Now I will filter the TCGA data to include only our relevant genes

BRCA_relevant_data <- BRCA_Initial_Data[BRCA_Initial_Data$Hugo_Symbol %in% rel_genes,]
nrow(BRCA_relevant_data) 

#How many variants on each gene
nvariants_relgenes <- table(BRCA_relevant_data$Hugo_Symbol)
nvariants_relgenes <- sort(nvariants_relgenes, decreasing = T); nvariants_relgenes[1:20] #Top 20 most mutated genes
                                                                                        #among our relevant genes so far


setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Initial_TCGA_Data/gdc_WESdata")
folders <- list.dirs()[2:length(list.dirs())]
full_cancer_data <- vector("list", length(folders))

maf_file1 <- list.files(folders[1])[list.files(folders[1]) != "annotations.txt"]

cancer_symbol <- unlist(strsplit(maf_file1,"[.]"))[2]; cancer_symbol #Get the cancer symbol

assign(cancer_symbol,read.table(file = paste(folders[1],maf_file1,sep="/"), 
                                comment.char = "", header = TRUE, sep = "\t", quote = "", skip = 5)) 

?assign

full_cancer_data[[1]] <- assign(cancer_symbol,read.table(file = paste(folders[1],maf_file1,sep="/"), 
                                                         comment.char = "", header = TRUE, sep = "\t", quote = "", skip = 5)) 
names(full_cancer_data)[1] <- cancer_symbol
#The previous steps created a dataframe whose name is the cancer symbol and whose content are the mutations found
#in that cancer's samples

for (f in 1:length(full_cancer_data)){
  maf_file <- list.files(folders[f])[list.files(folders[f]) != "annotations.txt"]
  cancer_symbol <- unlist(strsplit(maf_file,"[.]"))[2]
  full_cancer_data[[f]] <- read.table(file = paste(folders[f],maf_file,sep="/"), 
                                                           comment.char = "", header = TRUE, sep = "\t", quote = "", skip = 5)
  names(full_cancer_data)[f] <- cancer_symbol
  print(f)
}

#After all these tests I was able to create the function that loads the TCGA data into a list. It is a bit slow though.
#testing it:

setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Mutations_First_Work")
source("Mutations_First_Work_Functions.R")

start.time <- Sys.time()
TCGA_data <- load_TCGA_data("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Initial_TCGA_Data/gdc_WESdata")
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#Filtrar dataframes pelos genes relevantes (lapply)
#Gravar com função save. Carregar é load

#This will create a list with the information for all cancers pertaining exclusively to variants in the relevant
#glycosilation genes
TCGA_filtered_data <- lapply(TCGA_data, function(x) x[x$Hugo_Symbol %in% rel_genes, ])

lapply(TCGA_data, nrow)
lapply(TCGA_filtered_data, nrow)
mapply('/', lapply(TCGA_filtered_data, nrow), lapply(TCGA_data, nrow), SIMPLIFY = FALSE) 

barplot(unlist(mapply('/', lapply(TCGA_filtered_data, nrow), lapply(TCGA_data, nrow), SIMPLIFY = FALSE)),las=2)
?barplot

setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Mutations_First_Work")
save(TCGA_filtered_data, file = "../../RDataObjects/Mutations_First_Work/TCGA_October30filtered.Rdata")
save(TCGA_data, file = "../../RDataObjects/Mutations_First_Work/TCGA_October30unfiltered.Rdata")
save(TCGA_data, TCGA_filtered_data, file = "../../RDataObjects/Mutations_First_Work/TCGA_October30both.Rdata")
save(rel_genes, file = "../../RDataObjects/Mutations_First_Work/33relgenes.Rdata")

?save

#October Task 1 - Barplot impact for each cancer####
setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Mutations_First_Work")
source("Functions.R")
#Load the data
load("TCGA_both.Rdata")

TCGA_impact_data <- remove_low_impact(TCGA_data)


#Prepare the data - Full Data
cancer_labels <- names(TCGA_data)
impact_labels <- names(table(TCGA_data[[1]]$IMPACT))
impact_abs_matrix <- matrix(unlist(lapply(TCGA_data, function(x) table(x$IMPACT))), nrow = 4)
impact_rel_matrix <- sweep(impact_abs_matrix,2,colSums(impact_abs_matrix),'/')
library(RColorBrewer)

#Do the plot - Full Data
barplot(impact_rel_matrix, names.arg = cancer_labels, las=2, col = brewer.pal(nrow(impact_rel_matrix), "Paired"),
        xlim=c(0,ncol(impact_rel_matrix)+10.5), legend.text = impact_labels, ylab = "Relative occurence",
        main = "Distribution of mutation impact for each cancer - Full Data",
        args.legend = list(x = ncol(impact_rel_matrix)+15, y=max(colSums(impact_rel_matrix)), bty = "n"))
?barplot


#Prepare the Data - Filtered Data
impact_labelsf <- names(table(TCGA_filtered_data[[1]]$IMPACT))
impact_abs_matrixf <- matrix(unlist(lapply(TCGA_filtered_data, function(x) table(x$IMPACT))), nrow = 4)
impact_rel_matrixf <- sweep(impact_abs_matrixf,2,colSums(impact_abs_matrixf),'/')
library(RColorBrewer)

#Do the plot - Filtered Data
barplot(impact_rel_matrixf, names.arg = cancer_labels, las=2, col = brewer.pal(nrow(impact_rel_matrixf), "Paired"),
        xlim=c(0,ncol(impact_rel_matrixf)+10.5), legend.text = impact_labelsf, ylab = "Relative occurence",
        main = "Distribution of mutation impact for each cancer - Filtered Data",
        args.legend = list(x = ncol(impact_rel_matrixf)+15, y=max(colSums(impact_rel_matrixf)), bty = "n"))


#Ocober Task 2 - % of patients with mutations in glycosilation genes####

table(TCGA_data[[5]]$Tumor_Sample_Barcode)
table(TCGA_filtered_data[[5]]$Tumor_Sample_Barcode)

median(table(TCGA_filtered_data[[1]]$Tumor_Sample_Barcode)/table(TCGA_data[[1]]$Tumor_Sample_Barcode)*100)
table(TCGA_filtered_data[[5]]$Tumor_Sample_Barcode)/table(TCGA_data[[5]]$Tumor_Sample_Barcode)*100

lapply(TCGA_filtered_data, function(x) table(x$Tumor_Sample_Barcode))
lapply(TCGA_data, function(x) table(x$Tumor_Sample_Barcode))

medians <- unlist(lapply(mapply('/', lapply(TCGA_filtered_data, function(x) table(x$Tumor_Sample_Barcode)), 
       lapply(TCGA_data, function(x) table(x$Tumor_Sample_Barcode)), SIMPLIFY = FALSE), median),use.names = F)

boxplot(mapply('/', lapply(TCGA_filtered_data, function(x) table(x$Tumor_Sample_Barcode)), 
               lapply(TCGA_data, function(x) table(x$Tumor_Sample_Barcode)), SIMPLIFY = FALSE))


#split_n_get_first <- function(name){
#  return(strsplit(name, ".",fixed=T)[[1]][1])
#}
#c_labels <- sapply(names(medians),split_n_get_first, USE.NAMES = F)

barplot(medians, names.arg = cancer_labels, las=2, ylab = " Median percentage",
        main = "Median % of glycogene mutations per patient among all present mutations, per cancer")

mut_ratios <- mapply('/', lapply(TCGA_filtered_data, function(x) table(x$Tumor_Sample_Barcode)), 
       lapply(TCGA_data, function(x) table(x$Tumor_Sample_Barcode)), SIMPLIFY = FALSE)

first <- mut_ratios[[1]]
ffirst <- unname(first, force = T)
ffirst <- sort(ffirst)

xx <- seq(0,1,by=(1/(length(ffirst)-1))); xx

1/(length(ffirst)-1)

plot(xx, ffirst, xlim=c(0,20))

values <- c()
for (i in 1:length(mut_ratios)){
  new_values <- sort(unname(mut_ratios[[i]]))
  values <- c(values,new_values)
}

xx <- seq(0,1,by=(1/(length(values)-1))); xx

plot(xx, values, ylim=c(0,0.25), pch=20)


for (i in 1:length(TCGA_data)){
  symbols <- TCGA_data[[i]]$Hugo_Symbol
  #print(rel_genes %in% symbols)
}




#November Task 1 - Plot, for each gene, the percentage of patients with mutations, before and after filtering out low impact mutations####

setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Mutations_First_Work")
source("Mutations_First_Work_Functions.R")
library(RColorBrewer)
load("../../RDataObjects/Mutations_First_Work/TCGA_October30both.Rdata")
load("../../RDataObjects/Mutations_First_Work/33relgenes.Rdata")
#g <- rel_genes[2]
cancer_labels <- names(TCGA_data)

TCGA_data <- lapply(TCGA_data, function(x) cbind(x, PatientID=substring(x$Tumor_Sample_Barcode, 1,12)))
TCGA_filtered_data <- lapply(TCGA_filtered_data, function(x) cbind(x, PatientID=substring(x$Tumor_Sample_Barcode, 1,12)))
TCGA_impact_data <- remove_low_impact(TCGA_data)

npatients <- unlist(lapply(TCGA_data, function(x) nrow(table(x$PatientID))))

save(TCGA_impact_data, npatients, cancer_labels, file = "C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/RDataObjects/Mutations_First_Work/TCGA_November_impact.Rdata")

#t <- TCGA_data[[2]]


ratios_d <- unlist(lapply(TCGA_data, function(x) sum(table(x[x$Hugo_Symbol == g,]$PatientID))))/npatients*100
ratios_i <- unlist(lapply(TCGA_impact_data, function(x) sum(table(x[x$Hugo_Symbol == g,]$PatientID))))/npatients*100




#gg <- barplot(ratios_d, las=2, ylim=c(0,15),col = colorRampPalette(brewer.pal(length(ratios_d), "Accent"))(length(ratios_d)),
#        main = paste("% of patients with at least 1 mutation on gene",g,"for each cancer type"),
#        ylab="percentage of patients with mutation")

#barplot(ratios_i, las=2, ylim=c(0,15),col = colorRampPalette(brewer.pal(length(ratios_d), "Accent"))(length(ratios_d)))

#text(gg,ratios_d+2,labels=npatients,cex = 0.6)

pdf("../../Plots/Mutations_First_Work/Patients with mutation per cancer - 33 genes - November 6th.pdf")

par(mfrow=c(2,1))

for (i in 1:length(rel_genes)){
  
  ratios_d <- unlist(lapply(TCGA_data, function(x) length(unique(x[x$Hugo_Symbol == rel_genes[i],]$PatientID))))/npatients*100
  ratios_i <- unlist(lapply(TCGA_impact_data, function(x) length(unique(x[x$Hugo_Symbol == rel_genes[i],]$PatientID))))/npatients*100
  
  gg <- barplot(ratios_d, las=2, ylim=c(0,15),col = colorRampPalette(brewer.pal(length(ratios_d), "Accent"))(length(ratios_d)),
          main = paste(rel_genes[i],"(full data)"),
          ylab="patients with mutation (%)")
  
  grid(nx = NA, ny = NULL)
  
  if (i == 1){
    text(gg,ratios_d+2,labels=npatients,cex = 0.6)
  }
  
  ggg <- barplot(ratios_i, las=2, ylim=c(0,15),col = colorRampPalette(brewer.pal(length(ratios_i), "Accent"))(length(ratios_i)),
          main = paste(rel_genes[i],"no low-impact mutations"),
          ylab="patients with mutation (%)")
  
  grid(nx = NA, ny = NULL)
  
  if (i == 1){
    text(ggg,ratios_d+2,labels=npatients,cex = 0.6)
  }
}

par(mfrow=c(1,1))

dev.off()


# November Task 2 - Expand task 1#### 

#In order to present the data to Prof. Paula, and to try to extract some more biological information from the data ,
#I will do some variations on the graphs from Task 1. First, I will create a PDF exactly like the first one, but with
#only the graphs created from the filtered (no low impact) data. Then, instead of looking gene by gene, I will plot
#the same information but for all genes at once, that is, the percentage of patients with mutations on at least
#1 of the glycosilation genes. Then, finally, I will group the genes according to the step of the pathway they act on,
#and plot for each group.

#2.2

rel_genes_noF <- rel_genes[-length(rel_genes)]

ratios_i_2 <- unlist(lapply(TCGA_impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes_noF,]$PatientID))))/npatients*100

new_order <- order(ratios_i_2)
ratios_i_2 <- ratios_i_2[new_order]
npatients_neworder <- npatients[new_order]

pdf("../../Plots/Mutations_First_Work/Patients with mutation per cancer (no low impact) - ALL 33 genes combined - November 6th.pdf")

g3 <- barplot(ratios_i_2, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(ratios_i_2), "Accent"))(length(ratios_i_2)),
        main = "% of patients with mutated glycosilation genes, per cancer",
        ylab="patients with mutation (%)")

grid(nx = NA, ny = NULL)
text(g3,ratios_i_2+2,labels=npatients_neworder,cex = 0.6)

dev.off()

#2.1
pdf("../../Plots/Mutations_First_Work/Patients with mutation per cancer (no low impact) - 33 genes - November 6th.pdf")
par(mfrow=c(2,1))
for (i in 1:length(rel_genes)){
  ratios_i1 <- unlist(lapply(TCGA_impact_data, function(x) length(unique(x[x$Hugo_Symbol == rel_genes[i],]$PatientID))))/npatients*100
  ratios_i1 <- ratios_i1[new_order]
  gg <- barplot(ratios_i1, las=2, ylim=c(0,15),col = colorRampPalette(brewer.pal(length(ratios_i1), "Accent"))(length(ratios_i1)),
                main = rel_genes[i], ylab="patients with mutation (%)")
  grid(nx = NA, ny = NULL)
  if (i == 1){
    text(gg,ratios_i1+2,labels=npatients_neworder,cex = 0.6)
  }
}
par(mfrow=c(1,1))
dev.off()



#test <- lapply(TCGA_impact_data, function(x) table(x[x$Hugo_Symbol %in% rel_genes,]$PatientID))
#test <- unlist(test)
#tapply(test, names(test), sum)

#2.3

#Now I need to assign each gene name to its step in the metabolic pathway. To do that, I will simply name each element
#of the vector of genes with the step they belong to.

names(rel_genes) <- c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2",
                     "2","2","3","3","5","12","14","14","15", "FUT8")

pathway_genes <- split(rel_genes,names(rel_genes)); pathway_genes
pathway_genes <- pathway_genes[c("1","2","3","5","12","14","15","FUT8")]
pathway_steps <- names(pathway_genes)

#tr <- ratios_i3 <- unlist(lapply(TCGA_impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes[[8]],]$PatientID))))/npatients*100
#tr <- tr[new_order]
#tr
#g3 <- barplot(tr, las=2, ylim=c(0,15),col = colorRampPalette(brewer.pal(length(tr), "Accent"))(length(tr)),
#              main = pathway_steps[8], ylab="patients with mutation (%)")

pdf("../../Plots/Mutations_First_Work/Patients with mutation per cancer (no low impact) - 33 genes BY PATHWAY STEP - November 7th.pdf")

for (n in 1:length(pathway_genes)){
  
  ratios_i3 <- unlist(lapply(TCGA_impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% pathway_genes[[n]],]$PatientID))))/npatients*100
  ratios_i3 <- ratios_i3[new_order]
  
  g3 <- barplot(ratios_i3, las=2, ylim=c(0,40),col = colorRampPalette(brewer.pal(length(ratios_i3), "Accent"))(length(ratios_i3)),
                main = paste(pathway_steps[n], "-", length(pathway_genes[[n]]), "gene(s)"), ylab="patients with mutation (%)")
  grid(nx = NA, ny = NULL)
  if (n == 1){
    text(g3,ratios_i3+2,labels=npatients_neworder,cex = 0.6)
  }
  
}

dev.off()


#November Task 3 - Patients vs. Gene - Mutated or not - Heatmap####
setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Mutations_First_Work")
library(plyr)
source("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Mutations_First_Work/Mutations_First_Work_Functions.R")
library(RColorBrewer)
library(gplots)
load("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/RDataObjects/Mutations_First_Work/TCGA_November_impact.Rdata")
load("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/RDataObjects/Mutations_First_Work/33relgenes.Rdata")

#To do this heatmap, I first need to create a binary matrix, where the rows represent each gene and the columns each
#patient. The values will be 1 if the patient has a mutation on that gene, 0 if not.

#I will start by doing it for UCEC

#UCEC <- TCGA_impact_data[["UCEC"]]
patients <- as.character(unique(UCEC$PatientID))

matriz <- matrix(nrow = length(rel_genes), ncol = length(patients))
rownames(matriz) <- rel_genes

#mutated_patients <- unique(as.character(UCEC[UCEC$Hugo_Symbol == rel_genes[1],"PatientID"]))
#matriz[1,] <- sapply(patients, function(x) x %in% mutated_patients)


for (n in 1:length(rel_genes)){
  mutated_patients <- unique(UCEC[UCEC$Hugo_Symbol == rel_genes[n],"PatientID"])
  matriz[n,] <- as.integer(sapply(patients, function(x) x %in% mutated_patients))
}

#matriz <- matriz[, colSums(matriz != 0) > 0]
#matriz <- matriz[,order(colSums(matriz))]

path_steps <- c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2",
                "2","2","3","3","5","12","14","14","15", "FUT8")
path_colors <- mapvalues(path_steps,c("1","2","3","5","12","14","15","FUT8"),brewer.pal(8,"Set3"))

#heatmap(matriz, Rowv=NA, Colv=NA, scale="none", breaks=c(-1,0.1,1), col=c("white", "red"),RowSideColors = path_colors)

#matriz_s <- matriz[,1:100]
#heatmap(matriz_s, Rowv=NA, Colv=NA, scale="none", breaks=c(-1,0.1,1), col=c("white", "red"),RowSideColors = path_colors)

#palette cores
#FAzer para os cancros todos, 2 plots, um só com os 100 mais mutados

#display.brewer.all()

pdf("../../Plots/Mutations_First_Work/Cancer mutations Heatmaps - 12nd November.pdf")

for (n in 1:length(TCGA_impact_data)){
  
  #print(n)
  label <- names(TCGA_impact_data)[n]
  #print(label)
  cancer_data <- TCGA_impact_data[[n]]
  patients <- as.character(unique(cancer_data$PatientID))
  matriz <- matrix(nrow = length(rel_genes), ncol = length(patients))
  rownames(matriz) <- rel_genes
  for (m in 1:length(rel_genes)){
    mutated_patients <- unique(cancer_data[cancer_data$Hugo_Symbol == rel_genes[m],"PatientID"])
    matriz[m,] <- as.integer(sapply(patients, function(x) x %in% mutated_patients))
  }
  
  matriz <- matriz[, colSums(matriz != 0) > 0]
  if (is.null(ncol(matriz))){
    matriz <- cbind(matriz,rep(0,length(rel_genes)))
  }
  matriz <- matriz[,order(colSums(matriz))]
  
  path_steps <- c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2",
                  "2","2","3","3","5","12","14","14","15","FUT8")
  path_colors <- mapvalues(path_steps,c("1","2","3","5","12","14","15","FUT8"),brewer.pal(8,"Set3"))

  heatmap(matriz, Rowv=NA, Colv=NA, scale="none", breaks=c(-1,0.1,1), col=c("white", "red"),
          RowSideColors = path_colors, main = paste(label,"- All mutated patients"),cexRow = 0.8)

  
  
  if (ncol(matriz) > 100) {
    matriz_s <- matriz[,(ncol(matriz)-100):ncol(matriz)]
    heatmap(matriz_s, Rowv=NA, Colv=NA, scale="none", breaks=c(-1,0.1,1), col=c("white", "red"),
            RowSideColors = path_colors, main = paste(label,"- Top 100 mutated patients"),cexRow = 0.8)
  }
  
  
}



dev.off()

names(TCGA_impact_data)[n]

?heatmap.2

#mutação num dos genes no BRCA muito constante



#I will now repeat the heatmaps, but color the mutations by impact
'''

UCEC <- TCGA_impact_data[["SARC"]]
patients <- as.character(unique(UCEC$PatientID))
matriz <- matrix(0,nrow = length(rel_genes), ncol = length(patients))
colnames(matriz) <- patients
rownames(matriz) <- rel_genes
mutated_patientsi <- UCEC[UCEC$Hugo_Symbol == rel_genes[16],c("PatientID","IMPACT")]; mutated_patientsi

#This next bit of code will keep only the most severe mutation classification for each patient
mutated_patientsi$IMPACT <- factor(mutated_patientsi$IMPACT, levels = c("HIGH","MODIFIER","MODERATE"))
mutated_patientsi <- mutated_patientsi[order(mutated_patientsi$IMPACT),]
mutated_patientsi <- mutated_patientsi[!duplicated(mutated_patientsi[,"PatientID"]),]
#Replacing the impact names for numbers -> Moderate - 1; Modifier - 2; High - 3
mutated_patientsi$IMPACT <- mapvalues(mutated_patientsi$IMPACT, from=c("MODERATE", "MODIFIER", "HIGH"), 
                                      to=c("1", "2", "3"))

matriz[16,as.character(mutated_patientsi$PatientID)] <- as.numeric(levels(mutated_patientsi$IMPACT))[mutated_patientsi$IMPACT]


for (m in 1:length(rel_genes)){
  mutated_patientsi <- UCEC[UCEC$Hugo_Symbol == rel_genes[m],c("PatientID","IMPACT")]
  mutated_patientsi$IMPACT <- factor(mutated_patientsi$IMPACT, levels = c("HIGH","MODIFIER","MODERATE"))
  mutated_patientsi <- mutated_patientsi[order(mutated_patientsi$IMPACT),]
  mutated_patientsi <- mutated_patientsi[!duplicated(mutated_patientsi[,"PatientID"]),]
  mutated_patientsi$IMPACT <- mapvalues(mutated_patientsi$IMPACT, from=c("MODERATE", "MODIFIER", "HIGH"), 
                                        to=c("1", "2", "3"))
  matriz[m,as.character(mutated_patientsi$PatientID)] <- as.numeric(levels(mutated_patientsi$IMPACT))[mutated_patientsi$IMPACT]
}

matriz <- matriz[, colSums(matriz != 0) > 0]
matriz
matriz <- matriz[,order(colSums(matriz != 0))]

matriz <- apply(matriz, c(1,2), as.integer)
label <- "UCEC"
path_steps <- c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2",
                "2","2","3","3","5","12","14","14","15","FUT8")
path_colors <- mapvalues(path_steps,c("1","2","3","5","12","14","15","FUT8"),brewer.pal(8,"Set3"))


heatmap(matriz, Rowv=NA, Colv=NA, scale="none", breaks=c(-1,0.1,1.1,2.1,3.1), col=c("white", "yellow","orange","red"),
        RowSideColors = path_colors, main = paste(label,"- All mutated patients"),cexRow = 0.8)
'''

pdf("../../Plots/Mutations_First_Work/Cancer mutations Heatmaps Impact-Colored - 15th November.pdf")

for (n in 1:length(TCGA_impact_data)){
  
  label <- names(TCGA_impact_data)[n]
  cancer_data <- TCGA_impact_data[[n]]
  patients <- as.character(unique(cancer_data$PatientID))
  matriz <- matrix(0,nrow = length(rel_genes), ncol = length(patients))
  colnames(matriz) <- patients
  rownames(matriz) <- rel_genes
  
  for (m in 1:length(rel_genes)){
    mutated_patientsi <- cancer_data[cancer_data$Hugo_Symbol == rel_genes[m],c("PatientID","IMPACT")]
    mutated_patientsi$IMPACT <- factor(mutated_patientsi$IMPACT, levels = c("HIGH","MODIFIER","MODERATE"))
    mutated_patientsi <- mutated_patientsi[order(mutated_patientsi$IMPACT),]
    mutated_patientsi <- mutated_patientsi[!duplicated(mutated_patientsi[,"PatientID"]),]
    mutated_patientsi$IMPACT <- mapvalues(mutated_patientsi$IMPACT, from=c("MODERATE", "MODIFIER", "HIGH"), 
                                          to=c("1", "2", "3"))
    matriz[m,as.character(mutated_patientsi$PatientID)] <- as.numeric(levels(mutated_patientsi$IMPACT))[mutated_patientsi$IMPACT]
  }
  
  matriz <- matriz[, colSums(matriz != 0) > 0]
  if (is.null(ncol(matriz))){
    matriz <- cbind(matriz,rep(0,length(rel_genes)))
  }
  matriz <- matriz[,order(colSums(matriz != 0))]
  
  path_steps <- c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2",
                  "2","2","3","3","5","12","14","14","15","FUT8")
  path_colors <- mapvalues(path_steps,c("1","2","3","5","12","14","15","FUT8"),brewer.pal(8,"Set3"))
  
  heatmap(matriz, Rowv=NA, Colv=NA, scale="none", breaks=c(-1,0.1,1.1,2.1,3.1), col=c("white", "yellow2","orange","red"),
          RowSideColors = path_colors, main = paste(label,"- All mutated patients"),cexRow = 0.8)
  
  
  
  if (ncol(matriz) > 100) {
    matriz_s <- matriz[,(ncol(matriz)-100):ncol(matriz)]
    heatmap(matriz_s, Rowv=NA, Colv=NA, scale="none", breaks=c(-1,0.1,1.1,2.1,3.1), col=c("white", "yellow2","orange","red"),
            RowSideColors = path_colors, main = paste(label,"- Top 100 mutated patients"),cexRow = 0.8)
  }
}


dev.off()


#November Task 4 - Expand Task 2, include extra info####

#Some of the results of Task 1 and 2 were possibly promising. For the purpose of presenting them to Paula, it's 
#interesting to explore multiple ways of presenting them while giving for information about them. Specifically, the 
#plot with the percentage of patients with mutations on the glycosilation genes, for each cancer, revealed some 
#interesting leads. Some of the cancers had really high percentages of patients with mutated glycosilation genes.
#Is that a result worth exploring, or does it simply happen because those cancers have very high mutation rates anyway
#or a very low amount of patients. To try to answer this, we will repeat the plot that shows the percentage of
#patients with mutated glycosilation genes, but we will also include, in the same page, another plor with extra info,
#with one simply being the number of patients per cancer, another the number of total mutations divided by the
#number of patients, a third one with the number of mutations on glycosilation genes divided by the number of total 
#mutations, and maybe the fraction of total mutations that are in the glycosilation genes

#4.1 % patients mutated + n

rel_genes_noF <- rel_genes[-length(rel_genes)]

#Retrieve, for each cancer, the number of patients with at least 1 mutated gene of interest and divide it by the total
#number of patients for that cancer, multiplying by 100 to get a percentage.
ratios_4_1 <- unlist(lapply(TCGA_impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes_noF,]$PatientID))))/npatients*100

#For most graphs from now onwards, I want the columns with the cancer labels to be ordered from the one with the 
#lowest % of mutated patients to the highest. Here, I store the indexes of that new order.
new_order <- order(ratios_4_1)

#Here I re-order the ratios
ratios_4_1 <- ratios_4_1[new_order]

#Here I re-order the vector with the total number of patients per cancer, so that it matches the ratios one
npatients_neworder <- npatients[new_order]


pdf("../../Plots/Mutations_First_Work/Patients with mutation per cancer and Npatients - November 15th.pdf")
par(mfrow=c(2,1))

barplot(npatients_neworder, las=2, ylim=c(0,1000),col = colorRampPalette(brewer.pal(length(ratios_4_1), "Accent"))(length(ratios_4_1)),
        main = "total number of patients",
        ylab="number of patients")
grid(nx = NA, ny = NULL)


barplot(ratios_4_1, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(ratios_4_1), "Accent"))(length(ratios_4_1)),
            main = "% of patients with mutated glycosilation genes, per cancer",
            ylab="patients with mutation (%)")
grid(nx = NA, ny = NULL)

par(mfrow=c(1,1))
dev.off()



#4.2 % patients mutated + nmuttot/npatients

#Here I obtain the ratio between the number of total mutations and the number of patients, for each cancer type
ratios_4_2 <- unlist(lapply(TCGA_impact_data, function(x) nrow(x)))/npatients*100
ratios_4_2 <- ratios_4_2[new_order]

pdf("../../Plots/Mutations_First_Work/Patients with mutation per cancer and Nmuttot div Ndoentes - November 15th.pdf")
par(mfrow=c(2,1))

barplot(ratios_4_2, las=2, ylim=c(0,150000),col = colorRampPalette(brewer.pal(length(ratios_4_2), "Accent"))(length(ratios_4_2)),
        main = "Nmutations / Ndoentes",
        ylab="Nmutations/Npatients")
grid(nx = NA, ny = NULL)


barplot(ratios_4_1, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(ratios_4_1), "Accent"))(length(ratios_4_1)),
        main = "% of patients with mutated glycosilation genes, per cancer",
        ylab="patients with mutation (%)")
grid(nx = NA, ny = NULL)

par(mfrow=c(1,1))
dev.off()



#4.3 %patients mutated + nmutrelgenes/nmuttot

#Here I obtain the ratio between the number of total mutations and the number of patients, for each cancer type
ratios_4_3 <- unlist(lapply(TCGA_impact_data, function(x) nrow(x[x$Hugo_Symbol %in% rel_genes_noF,])))/unlist(lapply(TCGA_impact_data, function(x) nrow(x)))

ratios_4_3 <- ratios_4_3[new_order]

pdf("../../Plots/Mutations_First_Work/Patients with mutation per cancer and Nmut-relgenes div Nmut_total - November 19th.pdf")
par(mfrow=c(2,1))

barplot(ratios_4_3, las=2, ylim=c(0,0.30),col = colorRampPalette(brewer.pal(length(ratios_4_3), "Accent"))(length(ratios_4_3)),
        main = "Nmutations in relevant genes/N total mutations",
        ylab="Nmut_relgenes/Nmut_total")
grid(nx = NA, ny = NULL)

barplot(ratios_4_1, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(ratios_4_1), "Accent"))(length(ratios_4_1)),
        main = "% of patients with mutated glycosilation genes, per cancer",
        ylab="patients with mutation (%)")

grid(nx = NA, ny = NULL)

par(mfrow=c(1,1))
dev.off()



#4.4 %patients mutated + nmutuniquerelgenes/nmutuniquetot

ratios_4_4 <- unlist(lapply(TCGA_impact_data, function(x) nrow(unique(x[x$Hugo_Symbol %in% rel_genes_noF,c("Hugo_Symbol","PatientID")]))))/unlist(lapply(TCGA_impact_data, function(x) nrow(unique(x[,c("Hugo_Symbol","PatientID")]))))*100

ratios_4_4 <- ratios_4_4[new_order]


pdf("../../Plots/Mutations_First_Work/Patients with mut per cancer and Nmut_uniquerelgenes div Nmutunique_total - November22nd.pdf")
par(mfrow=c(2,1))

barplot(ratios_4_4, las=2, ylim=c(0,0.30),col = colorRampPalette(brewer.pal(length(ratios_4_4), "Accent"))(length(ratios_4_4)),
        main = "N unique mut. in rel. genes/N unique total mut.",
        ylab="Nmut_relgenes/Nmut_total")
grid(nx = NA, ny = NULL)

barplot(ratios_4_1, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(ratios_4_1), "Accent"))(length(ratios_4_1)),
        main = "% of patients with mutated glycosilation genes, per cancer",
        ylab="patients with mutation (%)")

grid(nx = NA, ny = NULL)

par(mfrow=c(1,1))
dev.off()



#4.5 %patients mutated + nrelgenesmutated/ntotgenesmutated

ratios_4_5 <- unlist(lapply(TCGA_impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes_noF,"Hugo_Symbol"]))))/unlist(lapply(TCGA_impact_data, function(x) length(unique(x[,"Hugo_Symbol"]))))*100

ratios_4_5 <- ratios_4_5[new_order]


pdf("../../Plots/Mutations_First_Work/Patients with mut per cancer and Nmutrelgenes div Ntotalmutgenes - November22nd.pdf")
par(mfrow=c(2,1))

barplot(ratios_4_5, las=2, ylim=c(0,0.30),col = colorRampPalette(brewer.pal(length(ratios_4_5), "Accent"))(length(ratios_4_5)),
        main = "Nmutrelgenes/Ntotalmutgenes.",
        ylab="Nmut_relgenes/Nmut_total")
grid(nx = NA, ny = NULL)

barplot(ratios_4_1, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(ratios_4_1), "Accent"))(length(ratios_4_1)),
        main = "% of patients with mutated glycosilation genes, per cancer",
        ylab="patients with mutation (%)")

grid(nx = NA, ny = NULL)

par(mfrow=c(1,1))
dev.off()


#FAZER MEDIAN NUMBER MUTATIONS PER PATIENT

symbols <- lapply(TCGA_data, function(x) x$Hugo_Symbol)
symbols <- unlist(symbols)
"ST8SIA6" %in% symbols
SKCM <- TCGA_data[["SKCM"]]
SKCMa <- SKCM[SKCM$Hugo_Symbol == "B4GALNT2",]
SKCMa$IMPACT

SKCM <- TCGA_impact_data[["SKCM"]]
SKCMa <- SKCM[SKCM$Hugo_Symbol == "B4GALNT2",]
SKCMa$IMPACT
