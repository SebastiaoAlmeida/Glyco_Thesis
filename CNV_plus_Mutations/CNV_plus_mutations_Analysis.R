source("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Functions/Mutations_First_Work_Functions.R")
source("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Functions/CNV_Functions.R")
source("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Functions/CNV_plus_Mutations_Functions.R")
library(gsubfn)

plot_directory <- "C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Plots/CNV+Mutations/" 
TCGA_data <- load_TCGA_data("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Initial_TCGA_Data/gdc_WESdata")
CNVdata <- load_cancer_CNVdata("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/CNV_Data")
old_rel_genes <- names(read.csv(file = "C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Gene_Lists/December14 38 genes list.csv", colClasses = "character")[1,])
rel_genes <- names(read.csv(file = "C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Gene_Lists/Project_Genes.csv", colClasses = "character")[1,])
prepared_data <- prepare_CNV_and_variant_data_matching(TCGA_data,CNVdata,old_rel_genes)

genomic_aterations_matrix <- return_matrix_and_plot_basic_genomic_alterations_heatmaps(prepared_data[[1]],prepared_data[[2]],plot_directory,"Original Genes")

