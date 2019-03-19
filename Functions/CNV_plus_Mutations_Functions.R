
###REQUIRES: - A list of dataframes (variant_data), each one named for a cancer symbol, where each dataframe
############ contains genomic variant information, each row a variant, each column an attribute with info
############ about that variant. A "Hugo_Symbol" and a "Tumor_Sample_Barcode" columns are mandatory.
############ - A list of matrixes (CNV_data), each one where the named rows contain gene symbols and the
############ named columns contain the tumour/patient barcode as it was in the CNVdata. For the values, 0
############ represents no CNV alteration, 1 amplification, 2 a "big" amplification, -1 a deletion, -2  "big" deletion.
###GOAL: Return the variant_data and the CNV_data, as the two elements of a list, ready for the matching.
######## the patient IDs with the same structure and the dataframes with the same order
prepare_CNV_and_variant_data_matching <- function(variant_data,CNV_data,relevant_genes){
  
  #Filter to heep only the relevant genes, and the columns we need
  variant_data <- lapply(variant_data, function(x) x[x$Hugo_Symbol %in% relevant_genes,c("Hugo_Symbol","Tumor_Sample_Barcode","IMPACT") ])
  #Remove the variants classified as "LOW" impact
  variant_data <- remove_low_impact(variant_data)
  #Substring the Tumor_Sample_Barcode to keep only the 16 characters to match the CNV data
  variant_data <- lapply(variant_data, function(x) {x$Tumor_Sample_Barcode <- substr(x$Tumor_Sample_Barcode,1,16);x})
  #Keep only one variant per gene/patient pair. After all, we only care if the gene is mutated or not
  variant_data <- lapply(variant_data, function(x) unique(x[,c("Hugo_Symbol","Tumor_Sample_Barcode")]))
  variant_data <- lapply(variant_data, function(x) {x$Hugo_Symbol <- as.character(x$Hugo_Symbol);x})
  #Order the dataframes so the order is the same as the CNV_data
  variant_data <- variant_data[order(names(variant_data))]
  
  CNV_data <- lapply(CNVdata, function(x) x[relevant_genes, ])
  #Replace . with - to match the IDs from mutation data
  CNV_data <- lapply(CNV_data, function(x) {colnames(x) <- gsub("\\.","-",colnames(x));x})
  #Keep only the first 16 characters to allow for the matching
  CNV_data <- lapply(CNV_data, function(x) {colnames(x) <- substr(colnames(x),1,16);x})
  
  return(list(variant_data,CNV_data))
}




###REQUIRES: - CNV data and mutations data after processing with the prepare_CNV_and_variant_data_matching function
############ - The directory the plot will be saved in (plot_directory)
############ - A label that will figure on the pdf file title (title_label)
###GOAL: Plot a heatmap where each column is a pacient and the rows are labeled each one for a gene. The colour
######## of each "cell" represents the state of each gene in each patient. White for no genomic alterations,
######## red for CNV amplification, blue for CNV deletion, green for a mutation and black for a combination of
######## a CNV alteration and a genomic mutation
return_matrix_and_plot_basic_genomic_alterations_heatmaps <- function(prepared_variant_data,prepared_CNV_data,plot_directory,title_label){
  
  colours <- c("blue","white","red","black","forestgreen","black")
  breaks <- c(-2.1,-0.1,0.1,2.1,4.1,5.1,7.1)
  
  merged_data <- prepared_CNV_data
  
  pdf(paste(plot_directory,"Basic Genomic Alterations Heatmaps - ",title_label,".pdf",sep=""))
  
  for (m in 1:length(prepared_variant_data)){
    thisMUTdata <- prepared_variant_data[[m]]
    cancer_symbol <- names(prepared_variant_data)[m]
    thisCNVdata <- prepared_CNV_data[[cancer_symbol]]
    for (j in 1:nrow(thisMUTdata)){
      gene <- thisMUTdata[j,"Hugo_Symbol"]
      patient <- thisMUTdata[j,"Tumor_Sample_Barcode"]
      if ((patient %in% colnames(thisCNVdata)) & gene %in% rownames(thisCNVdata)){
        thisCNVdata[gene,patient] <- thisCNVdata[gene,patient]+5
      }
    }
    #Maybe I should add a legend to the heatmap
    heatmap(thisCNVdata,Rowv = NA,Colv=NA,labCol=F,col=colours,breaks=breaks,scale="none",main=cancer_symbol)
    merged_data[[cancer_symbol]] <- thisCNVdata
  }
  
  dev.off()
  return(merged_data)
}



name <- function()
