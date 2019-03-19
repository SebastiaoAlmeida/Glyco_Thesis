source("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Functions/Mutations_First_Work_Functions.R")
library(RColorBrewer)
library(gplots)
library(plyr)

#November 22nd - Final Version of the first analysis on the TCGA mutation data####

setwd("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Mutations_First_Work")


#Load the list of relevant genes
relgenes_csvname <- "December14 38 genes list.csv"
rel_genes <- read.csv(file = paste("../../Datasets/Gene_lists/",relgenes_csvname,sep=""), colClasses = "character")
rel_genes <- names(rel_genes[1,])

#Load the TCGA mutation data onto a list of dataframes
TCGA_data <- load_TCGA_data("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Initial_TCGA_Data/gdc_WESdata")
#Remove the rows in each dataframe pertaining to variants in genes not in our list of relevant genes
TCGA_filtered_data <- lapply(TCGA_data, function(x) x[x$Hugo_Symbol %in% rel_genes, ])
#Remove the rows in each dataframe pertaining to variants labeled as "LOW" impact.
TCGA_impact_data <- remove_low_impact(TCGA_filtered_data)
#Create a dataframe with the variants from all genes, but remove low impact ones
TCGA_unfiltered_impact_data <- remove_low_impact(TCGA_data)

plot_directory <- "C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Plots/Mutations_First_Work/"

make_impact_barplot(TCGA_data, "Full Data",plot_directory)
make_impact_barplot(TCGA_filtered_data, "Filtered Data",plot_directory)

#The TCGA data samples have long barcodes that contain much information. For our analysis, its important
#to create a column containing just the subset of the barcode that constitutes the "PatientID".
TCGA_data <- lapply(TCGA_data, function(x) cbind(x, PatientID=substring(x$Tumor_Sample_Barcode, 1,12)))
TCGA_filtered_data <- lapply(TCGA_filtered_data, function(x) cbind(x, PatientID=substring(x$Tumor_Sample_Barcode, 1,12)))
TCGA_impact_data <- lapply(TCGA_impact_data, function(x) cbind(x, PatientID=substring(x$Tumor_Sample_Barcode, 1,12)))
#From this newly created column, we count the number of different patients in the data, for each cancer,
#creating a named vector with the number of different patients per cancer. We use the original dataframe
#to make sure all patients are included
npatients <- unlist(lapply(TCGA_data, function(x) nrow(table(x$PatientID))))

#This here is a particular adjustment for the gene set used in this analysis. One of the genes in the analysis
#was added because it was a gene of interest, but, when it comes to the metabolic pathways, is unrelated to all
#the other genes, which are all from the beginning of the O-glycosilation pathway. So, for some analysis, it is
#necessary to have the same gene set, but without FUT8.
rel_genes_noF <- rel_genes[-length(rel_genes)]
#For future plots, it is interesting to not only always present the different cancer types in the same order, but
#make it so the order reflects, from the least to most, the cancer types with a higher mutation ratio, which I
#defined simply as the percentage of patients for said cancer with at least 1 mutation in the early glycosilation genes.
#Here, I used the gene set with no FUT8.
cancer_mutation_ratios <- unlist(lapply(TCGA_impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes_noF,]$PatientID))))/npatients*100
cancer_order <- order(cancer_mutation_ratios)


plot_comparison_affected_patient_percentage_per_gene(TCGA_data,"Full Data",TCGA_impact_data,"No Low Impact Mutations",plot_directory,"Comparison of the percentage of patients with 1+ genetic variants when removing low-impact data",npatients,rel_genes,cancer_order)
#Analyzing the obtained plot lead to the conclusion that removing the low impact variants, as suspected,
#doesn't change much about the plots obtained, validating the choice to remove them.

plot_affected_patient_percentage_per_gene(TCGA_impact_data,plot_directory,npatients,rel_genes,cancer_order)

warnings()

plot_affected_patient_percentage_all_genes(TCGA_impact_data,plot_directory,npatients,cancer_order,rel_genes_noF)
plot_affected_patient_percentage_bypathwaystep(TCGA_impact_data,plot_directory,npatients,cancer_order,rel_genes)
heatmap_genes_patients_mutations_bycancer(TCGA_impact_data,plot_directory,rel_genes)
heatmap_genes_patients_mutations_bycancer(TCGA_impact_data,plot_directory,rel_genes,TRUE)
plot_percentage_mutated_and_n(TCGA_impact_data,plot_directory,npatients,cancer_order,rel_genes_noF)
plot_percentage_mutated_and_nmutovernpat(TCGA_impact_data,TCGA_unfiltered_impact_data,plot_directory,npatients,cancer_order,rel_genes_noF)
plot_percentage_mutated_and_nmutrelgenesovernmuttot(TCGA_impact_data,TCGA_unfiltered_impact_data,plot_directory,npatients,cancer_order,rel_genes_noF)
plot_percentage_mutated_and_nuniqmutrelgenesovernuniqmuttot(TCGA_impact_data,TCGA_unfiltered_impact_data,plot_directory,npatients,cancer_order,rel_genes_noF)
