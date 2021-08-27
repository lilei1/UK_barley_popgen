#!/usr/bin/env Rscript

#   This is a script that can plot the Fst in different chromosomes.
#   Written by Li Lei
#   August 25, 2017

#   Usage: ./Fst_plot.R <yourfile.Fst> <chrX>

#   Take command line arguments
#   Stores arguments into a vector
args <- commandArgs(trailingOnly = TRUE)

#   Function to read in data
readFile <- function(filename) {
  data.file <- read.delim(
    file = filename, # passed as an argument
    header = FALSE, # First line is a header
    sep="\t", #seperate with tab
    fill = TRUE, # Fill empty fields with NAs
    na.strings = "NA"
  )
  return(data.file)
}


#Plot the files:
plot.Fst <- function(data,chrs,xlimits,outDir,outName) {
  pdf(file = paste0(outDir, "/", outName, ".pdf"))
  plot(x = data$V2/1000000,
       y = data$V3,
       ylab="Fst",
       xlab="Physical position(Mb)",
       main=chrs, 
       cex=0.5,
       xlim = c(0, xlimits), 
       ylim = c(0, 0.8)
       )
  dev.off()
  }


#   Driver function
main <- function() {
  fst.file <- args[1] # file name end with .hierfstat
  chrome <- args[2] # chromesome
  x.label.limits <- as.numeric(args[3])
  outDir <- args[4]
  outName <- args[5]
  Fst <- readFile(filename = fst.file)
  plot_Fst <- plot.Fst(data = Fst, chrs = chrome, xlimits = x.label.limits, outDir=outDir, outName = outName)
}

main() # Run the program