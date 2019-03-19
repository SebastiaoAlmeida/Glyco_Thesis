
#Documentation#
load_cancer_CNVdata <- function(directory){
  setwd(directory)
  folders <- list.dirs()[2:length(list.dirs())]
  full_cancer_CNVdata <- vector("list", length(folders))
  for (f in 1:length(full_cancer_CNVdata)){
    cancer_symbol <- unlist(strsplit(unlist(strsplit(folders[f],"_"))[2],"-"))[1]
    df <- read.table(file = paste(folders[f],"all_thresholded.by_genes.txt",sep="/"), 
                     header = TRUE, sep = "\t", quote = "")
    matriz <- as.matrix(df[,-(1:3)])
    rownames(matriz) <- df$Gene.Symbol
    full_cancer_CNVdata[[f]] <- matriz
    names(full_cancer_CNVdata)[f] <- cancer_symbol
  }
  return(full_cancer_CNVdata)
}






