#####REQUIRES: A directory containing a folder for each type of cancer, each one containing the
###############.maf.gz file with the mutation data. 
#####GOAL: Returns a list of dataframes, each one named with the symbol of the cancer type whose information
###########it contains. Each row of the dataframe represents a variant, and each column contains a piece
###########of information about that variant.
load_TCGA_data <- function(directory){
  setwd(directory)
  folders <- list.dirs()[2:length(list.dirs())]
  full_cancer_data <- vector("list", length(folders))
  for (f in 1:length(full_cancer_data)){
    maf_file <- list.files(folders[f])[list.files(folders[f]) != "annotations.txt"]
    cancer_symbol <- unlist(strsplit(maf_file,"[.]"))[2]
    full_cancer_data[[f]] <- read.table(file = paste(folders[f],maf_file,sep="/"), 
                                        comment.char = "", header = TRUE, sep = "\t", quote = "", skip = 5)
    names(full_cancer_data)[f] <- cancer_symbol
    }
  return(full_cancer_data)
}

#####REQUIRES: A dataframe or a list of dataframes with genomic variant information, each row a variant,
############## and one of the columns named "IMPACT", with information about the predicted impact of said
############## variant, where "LOW" means a variant with predicted low or no impact.
#####GOAL: Return the same dataframe or list of dataframes without the variants predicted to have a low impact
remove_low_impact <- function(data){
  if (class(data) == "list"){
    return(lapply(data, function(x) x[x$IMPACT != "LOW",]))
  }
  return(data[data$IMPACT != "LOW",])
}

#####REQUIRES: A data frame (data) with genomic variant information, each row a variant, and one of the columns named "IMPACT", with information about the predicted impact of said variant. 
############## A label (data_label) that will feature in the title of the pdf, which is meant to identify the data said heatmap was made from
############## A directory where the pdf with the barplot will be stored (plot_directory)
############## The RColorBrewer library
#####GOAL: Print a barplot to a pdf where each bar represents data for one type of cancer, and said bar is
#####segmented into different colored segments, each one representing variant data with a specific impact label.
#####The relative frequence of each impact-level variant is measured in the y-axis
make_impact_barplot <- function(data, data_label, plot_directory){
  pdf(paste(plot_directory,"Variant Impact per Cancer - ",data_label,".pdf",sep=""))
  cancer_labels <- names(data)
  impact_labels <- names(table(data[[1]]$IMPACT))
  impact_abs_matrix <- matrix(unlist(lapply(data, function(x) table(x$IMPACT))), nrow = 4)
  impact_rel_matrix <- sweep(impact_abs_matrix,2,colSums(impact_abs_matrix),'/')
  barplot(impact_rel_matrix, names.arg = cancer_labels, las=2, col = brewer.pal(nrow(impact_rel_matrix), "Paired"),
          xlim=c(0,ncol(impact_rel_matrix)+12.5), legend.text = impact_labels, ylab = "Relative frequence",
          main = paste("Distribution of variant impact for each cancer -",data_label),
          args.legend = list(x = ncol(impact_rel_matrix)+18, y=max(colSums(impact_rel_matrix)), bty = "n"))
  dev.off()
}

#####REQUIRES: The RColorBrewer library
############### Two data frames, data1 and data2, with genomic variant information, each row a variant.
############### A name/label for each of the dataframes, data1label and data2label
############### The name you want to give to the pdf containing the plots (plot_name), and the directory it should be save in (plot_directory)
############### A vector with the number of patients for each cancer, in the same order as the cancer dataframes
############### relevant_genes -> a vector with the names of the genes to consider in the analysis
############### cancer_order -> A numeric vector with the desired position of each cancer in the original order
#####GOAL: Create a PDF with a page for each gene in the relevant_genes vector, where each page contains
#####two barplots, with a bar for each cancer type, that reflects the percentage of patients with that
#####type of cancer that contain at least one mutation in the gene in question.
plot_comparison_affected_patient_percentage_per_gene <- function(data1,data1label,data2,data2label,
                                                        plot_directory,plot_name,number_of_patients,relevant_genes,cancer_order){
  pdf(paste(plot_directory,plot_name,".pdf",sep=""))
  par(mfrow=c(2,1))
  for (i in 1:length(relevant_genes)){
    ratios_1 <- unlist(lapply(data1, function(x) length(unique(x[x$Hugo_Symbol == relevant_genes[i],]$PatientID))))/number_of_patients*100
    ratios_1 <- ratios_1[cancer_order]
    ratios_2 <- unlist(lapply(data2, function(x) length(unique(x[x$Hugo_Symbol == relevant_genes[i],]$PatientID))))/number_of_patients*100
    ratios_2 <- ratios_2[cancer_order]
    gg <- barplot(ratios_1, las=2, ylim=c(0,15),col = colorRampPalette(brewer.pal(length(ratios_1), "Accent"))(length(ratios_1)),
                  main = paste(relevant_genes[i]," (",data1label,")",sep=""), ylab="patients with 1+ variants (%)")
    grid(nx = NA, ny = NULL)
    number_of_patients_new_order <- number_of_patients[cancer_order]
    if (i == 1){text(gg,ratios_1+2,labels=number_of_patients_new_order,cex = 0.6)}
    ggg <- barplot(ratios_2, las=2, ylim=c(0,15),col = colorRampPalette(brewer.pal(length(ratios_2), "Accent"))(length(ratios_2)),
                   main = paste(relevant_genes[i]," (",data2label,")",sep=""),ylab="patients with 1+ variants (%)")
    grid(nx = NA, ny = NULL)
    if (i == 1){text(ggg,ratios_1+2,labels=number_of_patients_new_order,cex = 0.6)}
  }
  par(mfrow=c(1,1))
  dev.off()
}


#####REQUIRES: 
#####GOAL: 
plot_affected_patient_percentage_per_gene <- function(data,plot_directory, number_of_patients,relevant_genes,cancer_order){
  pdf(paste(plot_directory,"Percentage of patients with 1+ genetic variants per cancer, per gene.pdf",sep=""))
  par(mfrow=c(2,1))
  for (i in 1:length(relevant_genes)){
    ratios <- unlist(lapply(data, function(x) length(unique(x[x$Hugo_Symbol == relevant_genes[i],]$PatientID))))/number_of_patients*100
    ratios <- ratios[cancer_order]
    gg <- barplot(ratios, las=2, ylim=c(0,15),col = colorRampPalette(brewer.pal(length(ratios), "Accent"))(length(ratios)),
                  main = relevant_genes[i], ylab="patients with 1+ variants (%)")
    grid(nx = NA, ny = NULL)
    number_of_patients_neworder <- number_of_patients[cancer_order]
    if (i == 1){text(gg,ratios+2,labels=number_of_patients_neworder,cex = 0.6)}
  }
  par(mfrow=c(1,1))
  dev.off()
}


#####REQUIRES: 
#####GOAL: 
plot_affected_patient_percentage_all_genes <- function(data,plot_directory,number_of_patients,cancer_order,relevant_genes){
  ratios <- unlist(lapply(data, function(x) length(unique(x[x$Hugo_Symbol %in% relevant_genes,]$PatientID))))/number_of_patients*100
  ratios <- ratios[cancer_order]
  number_of_patients_neworder <- number_of_patients[cancer_order]
  pdf(paste(plot_directory,"Percentage of patients with 1+ genetic variants per cancer, all genes combined.pdf",sep=""))
  g3 <- barplot(ratios, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(ratios), "Accent"))(length(ratios)),
                main = "% of patients with 1+ variants in glycosilation genes, per cancer", ylab="patients with 1+ variants in glycosylation genes (%)")
  grid(nx = NA, ny = NULL)
  text(g3,ratios+2,labels=number_of_patients_neworder,cex = 0.6)
  dev.off()
}

#####REQUIRES: 
#####GOAL: 
plot_affected_patient_percentage_bypathwaystep <- function(impact_data,plot_directory,npatients,cancer_order,rel_genes){
  names(rel_genes) <- c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2;7","2;7","7","7",
                        "3","3","5","5","9","9","12","14","rm","rm","rm","rm","rm","rm")
  pathway_genes <- split(rel_genes,names(rel_genes))
  genes <- unlist(pathway_genes[[which(unlist(lapply(names(pathway_genes), function(x) grepl(";",x))))]])
  steps <- names(pathway_genes[[which(unlist(lapply(names(pathway_genes), function(x) grepl(";",x))))]])[1]
  steps <- unlist(strsplit(steps,";"))
  for(n in 1:length(steps)){
    for(m in 1:length(genes))
      to_append <- genes[m]
    names(to_append) <- steps[n]
    pathway_genes[[steps[n]]] <- append(pathway_genes[[steps[n]]],to_append)
  }
  pathway_genes[[which(unlist(lapply(names(pathway_genes), function(x) grepl(";",x))))]] <- NULL
  pathway_genes[["rm"]] <- NULL
  pathway_steps <- names(pathway_genes[c("1","12","14","2","3","5","7","9")])
  npatients_neworder <- npatients[cancer_order]
  pdf(paste(plot_directory,"Proportions of patients with at least one relevant impact mutation, by pathway step.pdf",sep=""))
  for (n in 1:length(pathway_genes)){
    ratios_i3 <- unlist(lapply(impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% pathway_genes[[n]],]$PatientID))))/npatients*100
    ratios_i3 <- ratios_i3[cancer_order]
    g4 <- barplot(ratios_i3, las=2, ylim=c(0,40),col = colorRampPalette(brewer.pal(length(ratios_i3), "Accent"))(length(ratios_i3)),
                  main = paste(pathway_steps[n], "-", length(pathway_genes[[n]]), "gene(s)"), ylab="patients with at least one variant in glycosylation genes (%)")
    grid(nx = NA, ny = NULL)
    if (n == 1){text(g4,ratios_i3+2,labels=npatients_neworder,cex = 0.6)}
  }
  print(ratios_i3)
  dev.off()}


#####what it does and requires - requer a library ####
#
heatmap_genes_patients_mutations_bycancer <- function(impact_data,plot_directory,rel_genes,colorbyimpact = FALSE){
  if (colorbyimpact) pdf(paste(plot_directory,"Cancer mutations Heatmaps impact-colored.pdf",sep="")) else pdf(paste(plot_directory,"Cancer mutations Heatmaps.pdf",sep=""))
  for (n in 1:length(TCGA_impact_data)){
    label <- names(impact_data)[n]
    cancer_data <- impact_data[[n]]
    patients <- as.character(unique(cancer_data$PatientID))
    matriz <- matrix(0,nrow = length(rel_genes), ncol = length(patients))
    colnames(matriz) <- patients
    rownames(matriz) <- rel_genes
    if(!colorbyimpact){
      for (m in 1:length(rel_genes)){
        mutated_patients <- unique(cancer_data[cancer_data$Hugo_Symbol == rel_genes[m],"PatientID"])
        matriz[m,] <- as.integer(sapply(patients, function(x) x %in% mutated_patients))}}
    else{for (m in 1:length(rel_genes)){
      mutated_patientsi <- cancer_data[cancer_data$Hugo_Symbol == rel_genes[m],c("PatientID","IMPACT")]
      mutated_patientsi$IMPACT <- factor(mutated_patientsi$IMPACT, levels = c("HIGH","MODIFIER","MODERATE"))
      mutated_patientsi <- mutated_patientsi[order(mutated_patientsi$IMPACT),]
      mutated_patientsi <- mutated_patientsi[!duplicated(mutated_patientsi[,"PatientID"]),]
      mutated_patientsi$IMPACT <- mapvalues(mutated_patientsi$IMPACT, from=c("MODERATE", "MODIFIER", "HIGH"), 
                                            to=c("1", "2", "3"))
      matriz[m,as.character(mutated_patientsi$PatientID)] <- as.numeric(levels(mutated_patientsi$IMPACT))[mutated_patientsi$IMPACT]}}
    
    matriz <- matriz[, colSums(matriz != 0) > 0]
    if (is.null(ncol(matriz))){matriz <- cbind(matriz,rep(0,length(rel_genes)))}
    matriz <- matriz[,order(colSums(matriz != 0))]
    path_steps <- c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2",
                    "2","2","3","3","5","12","14","14","15","FUT8")
    path_colors <- mapvalues(path_steps,c("1","2","3","5","12","14","15","FUT8"),brewer.pal(8,"Set3"))
    if (colorbyimpact) assign("col",c("white", "yellow2","orange","red")) else assign("col",c("white","red"))
    if (colorbyimpact) assign("breaks",c(-1,0.1,1.1,2.1,3.1)) else assign("breaks",c(-1,0.1,1))
    heatmap(matriz, Rowv=NA, Colv=NA, scale="none", breaks=breaks, col=col,
            RowSideColors = path_colors, main = paste(label,"- All mutated patients"),cexRow = 0.8)
    if (ncol(matriz) > 100) {matriz_s <- matriz[,(ncol(matriz)-100):ncol(matriz)]
      heatmap(matriz_s, Rowv=NA, Colv=NA, scale="none", breaks=breaks, col=col,
              RowSideColors = path_colors, main = paste(label,"- Top 100 mutated patients"),cexRow = 0.8)
    }
  }
  dev.off()}

#####what it does and requires - requer a library ####
#
plot_percentage_mutated_and_n <- function(impact_data,plot_directory,npatients,cancer_order,rel_genes){
  mutation_ratios <- unlist(lapply(impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes,]$PatientID))))/npatients*100
  mutation_ratios <- mutation_ratios[cancer_order]
  npatients_neworder <- npatients[cancer_order]
  pdf(paste(plot_directory,"Npatients and percentages of patients with mutation, per cancer.pdf",sep=""))
  par(mfrow=c(2,1))
  barplot(npatients_neworder, las=2, ylim=c(0,1.2*max(npatients_neworder)),col = colorRampPalette(brewer.pal(length(mutation_ratios), "Accent"))(length(mutation_ratios)),
          main = "total number of patients",ylab="number of patients")
  grid(nx = NA, ny = NULL)
  barplot(mutation_ratios, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(mutation_ratios), "Accent"))(length(mutation_ratios)),
          main = "% of patients with mutated glycosilation genes, per cancer",ylab="patients with mutation (%)")
  grid(nx = NA, ny = NULL)
  par(mfrow=c(1,1))
  dev.off()}


#####what it does and requires - requer a library ####
#
plot_percentage_mutated_and_nmutovernpat <- function(impact_data,unfiltered_impact_data,plot_directory,npatients,cancer_order,rel_genes){
  mutation_ratios <- unlist(lapply(impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes,]$PatientID))))/npatients*100
  mutation_ratios <- mutation_ratios[cancer_order]
  ratios_mut_tot <- unlist(lapply(unfiltered_impact_data, function(x) nrow(x)))/npatients
  ratios_mut_tot <- ratios_mut_tot[cancer_order]
  pdf(paste(plot_directory,"Nmuttot over Npatients and percentages of patients with mutation, per cancer.pdf",sep=""))
  par(mfrow=c(2,1))
  plot <- barplot(ratios_mut_tot, ylim=c(0,1.2*max(ratios_mut_tot)), las=2,col = colorRampPalette(brewer.pal(length(mutation_ratios), "Accent"))(length(mutation_ratios)),
          main = "N impact mutations / Npatients", ylab="N impact mutations / Npatients")
  grid(nx = NA, ny = NULL)
  text(plot,ratios_mut_tot+max(ratios_mut_tot)/12,labels=npatients[cancer_order],cex = 0.6)
  barplot(mutation_ratios, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(mutation_ratios), "Accent"))(length(mutation_ratios)),
          main = "% of patients with mutated glycosilation genes, per cancer", ylab="patients with mutation (%)")
  grid(nx = NA, ny = NULL)
  par(mfrow=c(1,1))
  dev.off()}

#####what it does and requires - requer a library ####
#
plot_percentage_mutated_and_nmutrelgenesovernmuttot <- function(impact_data,unfiltered_impact_data,plot_directory,npatients,cancer_order,rel_genes){
  mutation_ratios <- unlist(lapply(impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes,]$PatientID))))/npatients*100
  mutation_ratios <- mutation_ratios[cancer_order]
  ratios_rel_mut <- unlist(lapply(impact_data, function(x) nrow(x[x$Hugo_Symbol %in% rel_genes_noF,])))/unlist(lapply(unfiltered_impact_data, function(x) nrow(x)))
  ratios_rel_mut <- ratios_rel_mut[cancer_order]
  pdf(paste(plot_directory,"Nmutrelgenes over Nmuttotal and percentages of patients with mutation, per cancer.pdf",sep=""))
  par(mfrow=c(2,1))
  plot <- barplot(ratios_rel_mut, las=2,ylim=c(0,1.2*max(ratios_rel_mut)),col = colorRampPalette(brewer.pal(length(mutation_ratios), "Accent"))(length(mutation_ratios)),
          main = "N mutations relgenes / Nmutations total", ylab="N mut relgenes / N mut tot")
  grid(nx = NA, ny = NULL)
  text(plot,ratios_rel_mut+max(ratios_rel_mut)/12,labels=npatients[cancer_order],cex = 0.6)
  barplot(mutation_ratios, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(mutation_ratios), "Accent"))(length(mutation_ratios)),
          main = "% of patients with mutated glycosilation genes, per cancer", ylab="patients with mutation (%)")
  grid(nx = NA, ny = NULL)
  par(mfrow=c(1,1))
  dev.off()}


#####what it does and requires - requer a library ####
#
plot_percentage_mutated_and_nuniqmutrelgenesovernuniqmuttot <- function(impact_data,unfiltered_impact_data,plot_directory,npatients,cancer_order,rel_genes){
  mutation_ratios <- unlist(lapply(impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes,]$PatientID))))/npatients*100
  mutation_ratios <- mutation_ratios[cancer_order]
  ratios_unique_rel_mut <- unlist(lapply(TCGA_impact_data, function(x) nrow(unique(x[x$Hugo_Symbol %in% rel_genes_noF,c("Hugo_Symbol","PatientID")]))))/unlist(lapply(TCGA_unfiltered_impact_data, function(x) nrow(unique(x[,c("Hugo_Symbol","PatientID")]))))
  ratios_unique_rel_mut <- ratios_unique_rel_mut[cancer_order]
  print(ratios_unique_rel_mut)
  pdf(paste(plot_directory,"Nuniquemutrelgenes over Nuniquemuttotal and percentages of patients with mutation, per cancer.pdf",sep=""))
  par(mfrow=c(2,1))
  plot <- barplot(ratios_unique_rel_mut, las=2,ylim=c(0,1.2*max(ratios_unique_rel_mut)),col = colorRampPalette(brewer.pal(length(mutation_ratios), "Accent"))(length(mutation_ratios)),
                  main = "N unique mutations relgenes / N unique mutations total", ylab="N uniqmut relgenes / N uniqmut tot")
  grid(nx = NA, ny = NULL)
  text(plot,ratios_unique_rel_mut+max(ratios_unique_rel_mut)/12,labels=npatients[cancer_order],cex = 0.6)
  barplot(mutation_ratios, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(mutation_ratios), "Accent"))(length(mutation_ratios)),
          main = "% of patients with mutated glycosilation genes, per cancer", ylab="patients with mutation (%)")
  grid(nx = NA, ny = NULL)
  par(mfrow=c(1,1))
  dev.off()}



#####what it does and requires - requer a library ####
#
plot_percentage_mutated_and_nrelgenesmutoverntotgenemut <- function(impact_data,unfiltered_impact_data,plot_directory,npatients,cancer_order,rel_genes){
  mutation_ratios <- unlist(lapply(impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes,]$PatientID))))/npatients*100
  mutation_ratios <- mutation_ratios[cancer_order]
  gene_ratios <- unlist(lapply(TCGA_impact_data, function(x) length(unique(x[x$Hugo_Symbol %in% rel_genes_noF,"Hugo_Symbol"]))))/unlist(lapply(TCGA_impact_data, function(x) length(unique(x[,"Hugo_Symbol"]))))
  gene_ratios <- gene_ratios[cancer_order]
  print(ratios_unique_rel_mut)
  pdf(paste(plot_directory,"Nuniquemutrelgenes over Nuniquemuttotal and percentages of patients with mutation, per cancer.pdf",sep=""))
  par(mfrow=c(2,1))
  plot <- barplot(gene_ratios, las=2,ylim=c(0,1.2*max(gene_ratios)),col = colorRampPalette(brewer.pal(length(mutation_ratios), "Accent"))(length(mutation_ratios)),
                  main = "N unique mutations relgenes / N unique mutations total", ylab="N uniqmut relgenes / N uniqmut tot")
  grid(nx = NA, ny = NULL)
  text(plot,gene_ratios+max(gene_ratios)/12,labels=npatients[cancer_order],cex = 0.6)
  barplot(mutation_ratios, las=2, ylim=c(0,50),col = colorRampPalette(brewer.pal(length(mutation_ratios), "Accent"))(length(mutation_ratios)),
          main = "% of patients with mutated glycosilation genes, per cancer", ylab="patients with mutation (%)")
  grid(nx = NA, ny = NULL)
  par(mfrow=c(1,1))
  dev.off()}

