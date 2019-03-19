source("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Functions/Mutations_First_Work_Functions.R")
source("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Code/Functions/Survival_Analysis_Functions.R")
library("survival")
library("survminer")
library(plyr)
library(dplyr)
library(data.table)
library(stats)


#Step 1 - Load Clinical Data
clinical_data <- load_clinical_TCGA_data("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Clinical_Data")
variant_data <- load_TCGA_data("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Initial_TCGA_Data/gdc_WESdata")
variant_data <- remove_low_impact(variant_data)
variant_data <- lapply(variant_data, function(x) cbind(x, PatientID=substring(x$Tumor_Sample_Barcode, 1,12)))
relgenes_csvname <- "December14 38 genes list.csv"
rel_genes <- read.csv(file = paste("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Datasets/Initial_TCGA_Data/",relgenes_csvname,sep=""), colClasses = "character")
rel_genes <- names(rel_genes[1,])

#I have to remove from the clinical data the cancers with no variant data
key <- names(clinical_data) %in% names(variant_data)
clinical_data <- clinical_data[key]

#I need to order the data alphabetically
variant_data <- variant_data[order(names(variant_data))]

#pdf("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Plots/Survival_Analysis/Test.pdf")
rel_pvalues <- rep(NA,2000)
r <- 1
p_values <- data.frame("pvalue"=rep(NULL,2000),"descr"=rep(NULL,2000))

pdf("C:/Users/Sebastião Almeida/Desktop/Tese/Main_Work/Plots/Survival_Analysis/TestMaior5.pdf")

for (c in 1:length(variant_data)){
  for (g in 1:length(rel_genes)){
    clinical_df <- clinical_data[[c]];
    variant_df <- variant_data[[c]]
    clin_patient_id <- clinical_df$bcr_patient_barcode
    variant_patient_id <- tolower(unique(as.character(variant_df$PatientID)))
    patient_id <- clin_patient_id[clin_patient_id %in% variant_patient_id]
    curr_gene <- rel_genes[g]
    rownames(clinical_df) <- clinical_df$bcr_patient_barcode
    time <- apply(clinical_df[patient_id, c("days_to_death","days_to_last_followup")], 1, function(x) sum(as.numeric(na.omit(x))))/365
    mutated <- rep(5,length(patient_id))
    for (n in 1:length(patient_id)){
      genes <- variant_df[variant_df$PatientID == toupper(patient_id[n]), "Hugo_Symbol"]
      mutated[n] <- as.numeric(curr_gene %in% as.character(genes))
    }
    if(sum(mutated) < 6){
      next
    }
    data <- data.frame("time"=time, "status"=as.numeric(mapvalues(clinical_df[patient_id ,"vital_status"], c("alive", "dead"),c(0,1))), "group"=mutated, stringsAsFactors = FALSE)
    surv_data <- Surv(data$time, event=data$status)
    surv_fit <- survfit(surv_data~group,data=data)
    tit <- paste("CANCER TYPE - ", names(clinical_data)[c], "; GENE - ",curr_gene,sep = "")
    g <- ggsurvplot(surv_fit,data=data,pval = TRUE,title=tit,palette=c("lightblue3","red"))
    print(g)
    #survdiff(surv_data~group,data=data)
    pvalue <- surv_pvalue(surv_fit, data = data, method = "survdiff")[1,2]
    p_values[r,1] <- pvalue
    p_values[r,2] <- tit
    r <- r+1
  }
}

dev.off()

rel_pvalues <- na.omit(rel_pvalues)
rel_pvalues

adjusted <- p.adjust(p_values$V1,"BH")
sum(na.omit(p_values$V1 < 0.05))
sum(na.omit(adjusted < 0.05))
min(na.omit(adjusted))

####################################################################################################
#Survival analisys comparing patients with early stage (i and ii) and late stage (iii and iv) cancer

clinical_info <- clinical_data[[5]]
variantdfteste <- variant_data[[5]]
curr_gene <- "GALNT8"
rownames(clinical_info) <- clinical_info$bcr_patient_barcode

#Remove the patients with tumor stage x and NA
clinical_info <- clinical_info[!(clinical_info$pathologic_stage %in% c("stage x", NA)),]

clin_patient_id_teste <- clinical_info$bcr_patient_barcode
variant_patient_id_teste <- tolower(unique(as.character(variantdfteste$PatientID)))
patient_id <- clin_patient_id_teste[clin_patient_id_teste %in% variant_patient_id_teste]

#
mutated <- rep(5,length(patient_id))
#Create a vector that, for each patient, tells us whether they are mutated or not for our gene
for (n in 1:length(patient_id)){
  genes <- variantdfteste[variantdfteste$PatientID == toupper(patient_id[n]), "Hugo_Symbol"]
  mutated[n] <- as.numeric(curr_gene %in% as.character(genes))
}

#Divide the samples into early stage and late stage tumors
tumorStage <- clinical_info[patient_id,"pathologic_stage"]
tumorStage <- replace(tumorStage,tumorStage %in% c("stage iia","stage iib","stage ia","stage i",
                                                   "stage ib", "stage ii"),"Early Stage")
tumorStage <- replace(tumorStage,tumorStage %in% c("stage iiic","stage iiia","stage iiib","stage iv"),"Late Stage")

data <- data.frame(mutated, tumorStage) 
contTable <- table(data$mutated, data$tumorStage)
contTable
fisher_result <- fisher.test(t(contTable))



##############Stages I-II vs. Stages III-IV#############################################################

r <- 1
p_values <- data.frame("pvalue"=rep(NULL,2000),"descr"=rep(NULL,2000))

for (c in 1:length(variant_data)){

  for (g in 1:length(rel_genes)){
    clinical_df <- clinical_data[[c]]
    variant_df <- variant_data[[c]]
    curr_gene <- rel_genes[g]
    rownames(clinical_df) <- clinical_df$bcr_patient_barcode
    clinical_df <- clinical_df[!(clinical_df$pathologic_stage %in% c("stage x", NA)),]
    clin_patient_id <- clinical_df$bcr_patient_barcode
    variant_patient_id <- tolower(unique(as.character(variant_df$PatientID)))
    patient_id <- clin_patient_id[clin_patient_id %in% variant_patient_id]
    mutated <- rep(5,length(patient_id))
    for (n in 1:length(patient_id)){
      genes <- variant_df[variant_df$PatientID == toupper(patient_id[n]), "Hugo_Symbol"]
      mutated[n] <- as.numeric(curr_gene %in% as.character(genes))
    }
    if(sum(mutated) < 6){
      next
    }
    tumorStage <- clinical_df[patient_id,"pathologic_stage"]
    tumorStage <- replace(tumorStage,tumorStage %in% c("stage iia","stage iib","stage ia","stage i",
                                                       "stage ib", "stage ii"),"Early Stage")
    tumorStage <- replace(tumorStage,tumorStage %in% c("stage iiic","stage iiia","stage iiib","stage iv"),"Late Stage")
    data <- data.frame(mutated, tumorStage) 

    contTable <- table(data$mutated, data$tumorStage)
    fisher_result <- fisher.test(t(contTable))
    pvalue <- fisher_result[[1]]
    tit <- paste("CANCER TYPE - ", names(clinical_data)[c], "; GENE - ",curr_gene,sep = "")
    p_values[r,1] <- pvalue
    p_values[r,2] <- tit
    r <- r+1
   
  }
}

adjusted <- p.adjust(p_values$V1,"BH")
min(adjusted)


##############Stage IV vs. Stages I-III################################################################

r <- 1
p_values2 <- data.frame("pvalue"=rep(NULL,2000),"descr"=rep(NULL,2000))

for (c in 1:length(variant_data)){
  
  for (g in 1:length(rel_genes)){
    clinical_df <- clinical_data[[c]]
    variant_df <- variant_data[[c]]
    curr_gene <- rel_genes[g]
    rownames(clinical_df) <- clinical_df$bcr_patient_barcode
    clinical_df <- clinical_df[!(clinical_df$pathologic_stage %in% c("stage x", NA)),]
    clin_patient_id <- clinical_df$bcr_patient_barcode
    variant_patient_id <- tolower(unique(as.character(variant_df$PatientID)))
    patient_id <- clin_patient_id[clin_patient_id %in% variant_patient_id]
    mutated <- rep(5,length(patient_id))
    for (n in 1:length(patient_id)){
      genes <- variant_df[variant_df$PatientID == toupper(patient_id[n]), "Hugo_Symbol"]
      mutated[n] <- as.numeric(curr_gene %in% as.character(genes))
    }
    if(sum(mutated) < 6){
      next
    }
    tumorStage <- clinical_df[patient_id,"pathologic_stage"]
    tumorStage <- replace(tumorStage,tumorStage %in% c("stage iia","stage iib","stage ia","stage i",
                                                       "stage ib", "stage ii","stage iiic","stage iiia","stage iiib"),"Early Stage (1-3)")
    data <- data.frame(mutated, tumorStage) 
    
    contTable <- table(data$mutated, data$tumorStage)
    fisher_result <- fisher.test(t(contTable))
    pvalue <- fisher_result[[1]]
    tit <- paste("CANCER TYPE - ", names(clinical_data)[c], "; GENE - ",curr_gene,sep = "")
    p_values2[r,1] <- pvalue
    p_values2[r,2] <- tit
    r <- r+1
    
  }
}

adjusted <- p.adjust(p_values2$V1,"BH")
min(adjusted)



























##############Drafts############################################################################

clinical_info <- clinical_data[[3]]
variantdfteste <- variant_data[[3]]

clin_patient_id_teste <- clinical_info$bcr_patient_barcode
variant_patient_id_teste <- tolower(unique(as.character(variantdfteste$PatientID)))

p <- clin_patient_id_teste[clin_patient_id_teste %in% variant_patient_id_teste]
p %in% variant_patient_id_teste
variant_patient_id_teste %in% clin_patient_id_teste

patient_id <- clin_patient_id_teste

curr_gene <- "GALNT8"

rownames(dfteste) <- dfteste$bcr_patient_barcode

time <- apply(dfteste[patient_id, c("days_to_death","days_to_last_followup")], 1, function(x) sum(as.numeric(na.omit(x))))/365

mutated <- rep(5,length(patient_id))

for (n in 1:length(patient_id)){
  genes <- variantdfteste[variantdfteste$PatientID == toupper(patient_id[n]), "Hugo_Symbol"]
  mutated[n] <- as.numeric(curr_gene %in% as.character(genes))
}

gene

#str(dfteste)

data <- data.frame("time"=time, "status"=as.numeric(mapvalues(dfteste[patient_id ,"vital_status"], c("alive", "dead"),c(0,1))), "group"=mutated, stringsAsFactors = FALSE)

#str(data)

#?data.frame

surv_data <- Surv(data$time, event=data$status)
surv_fit <- survfit(surv_data~group,data=data)
ggsurvplot(surv_fit,data=data,pval = TRUE, title = "hey")
#survdiff(surv_data~group,data=data)
pvalue <- surv_pvalue(surv_fit, data = data, method = "survdiff")[1,4]
pvalue <- as.numeric(substr(pvalue,nchar(pvalue)-3,nchar(pvalue)))

?ggsurvplot

dev.off()





