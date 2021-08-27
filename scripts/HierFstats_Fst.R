#!/usr/bin/env Rscript

#   This is a script that calculate the Fst with HierFstats
#   Written by Li Lei
#   August 25, 2017

#   Usage: ./R_hierfst_calculation.R <yourfile.hierfstat>
###calculate Fst using the
library(hierfstat)
#   Take command line arguments
#   Stores arguments into a vector
args <- commandArgs(trailingOnly = TRUE)

#   Function to read in data
readFile <- function(filename) {
  data.file <- read.delim(
    file = filename, # passed as an argument
    header = TRUE, # First line is a header
    sep="\t", #seperate with tab
    fill = TRUE, # Fill empty fields with NAs
    na.strings = "NA"
  )
  return(data.file)
}

#   Function to write out file
writeOutFile <- function(data, outName) {
  write.table(x = data,
              file = outName,
              quote = FALSE,
              sep = "\t",
              eol = "\n",
              col.names = TRUE,
              row.names = FALSE)
}

#   Driver function
main <- function() {
  fst.file <- args[1] # file name end with .hierfstat
  outName_base <- args[2] # Outpstrsplit("a.b.c", ".", fixed = TRUE)ut filename
  data.fst <- readFile(filename = fst.file)
  loci <- data.fst[,-c(1,2)] #get all of the loci information
  Fst<-NULL
  for (i in 1:dim(loci)[[2]]) {
           Fst[i]<-varcomp(data.frame(data.fst[,2],as.numeric(loci[,i])),diploid=FALSE)$F[1,1]
  }
  fst.df <- as.data.frame(Fst)
  rownames(fst.df)<-colnames(loci)
  fst.df$PhyPos <- rownames(fst.df)
  rownames(fst.df) <- NULL
  library(splitstackshape)
  fst.df<-cSplit(fst.df,"PhyPos",".")
  fst.df$chr <- fst.df$PhyPos_1
  fst.df$pos <- fst.df$PhyPos_2
  fst.df$PhyPos_1 <- NULL
  fst.df$PhyPos_2 <- NULL
  fst.df$PhyPos_3 <- NULL
  row.names(fst.df) <- NULL
  fst.df <- data.frame(cbind(as.character(fst.df$chr),fst.df$pos,fst.df$Fst))
  colnames(fst.df) <- c("Chr","Pos","FST")
  row.names(fst.df) <- NULL
  writeOutFile(data = fst.df, outName = paste0(outName_base,".Fst"))
  #writeOutFile(data = overall.fst, outName = paste0(outName_base,".overall.Fst"))
 # writeOutFile(data = F.fst, outName = paste0(outName_base,".F.fst"))

}

main() # Run the program