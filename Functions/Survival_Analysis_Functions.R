#####what it does and requires####
####
load_clinical_TCGA_data <- function(directory){
  setwd(directory)
  folders <- list.dirs()[2:length(list.dirs())]
  clinical_data <- vector("list", length(folders))
  for (f in 1:length(clinical_data)){
    txtfile <- list.files(folders[f])[list.files(folders[f]) == "All_CDEs.txt"]
    cancer_symbol <- unlist(strsplit(unlist(strsplit(folders[f],"[.]"))[4],"_"))[[2]]
    
    data <- read.table(file = paste(folders[f],txtfile,sep="/"), 
                        comment.char = "", row.names = 1, header = F, sep = "\t", quote = "")
    fdata <- transpose(data)
    colnames(fdata) <- rownames(data)
    clinical_data[[f]] <- fdata
    names(clinical_data)[f] <- cancer_symbol
  }
  return(clinical_data)
}

