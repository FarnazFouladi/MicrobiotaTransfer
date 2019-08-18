#Modded from https://benjjneb.github.io/dada2/tutorial.html
#Farnaz Fouladi
#08-18-2019

#Generating sequence variant table using Dada2 pipeline

rm(list=ls())
library(dada2)

folders<-c("AN34","AN40","AN703","AN81")
for (folder in folders){
  
  path <- paste0("/users/ffouladi/anorexiaMouseTransfer/demultiplexedR1/primerremoved/",folder)
  
  fnFs <- sort(list.files(path, pattern="_R1.fastq.filtered.fastq"))
  
  sample.names <- gsub("_R1.fastq.filtered.fastq", "", fnFs, fixed=TRUE)
  
  fnFs <- file.path(path, fnFs)
  
  filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
  
  filtFs <- file.path(filt_path, paste0(sample.names, "_R1_filt.fastq"))
  
  out <- filterAndTrim(fnFs, filtFs, ,truncLen=200,
                       maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=FALSE) 
  
  filtpath <- paste0("/users/ffouladi/anorexiaMouseTransfer/demultiplexedR1/primerremoved/",folder,"/filtered")
  
  filts <- list.files(filtpath, pattern="fastq", full.names=TRUE) # CHANGE if different file extensions
  
  names(filts) <- sample.names
  
  dds <- vector("list", length(sample.names))
  names(dds) <- sample.names
  
  index <-1 
  setwd("/users/ffouladi/anorexiaMouseTransfer/dada2result")
  for (f in filtFs){
    errF <- learnErrors(f, multithread = TRUE)
    derepFs <- derepFastq(f, verbose=TRUE)
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
    dds[[index]]<-dadaFs
    index<-index+1
  }
  seqtab <- makeSequenceTable(dds)
  
  saveRDS(seqtab, paste0("/users/ffouladi/anorexiaMouseTransfer/dada2result/",folder,"ForwardsReads.rds"))
  write.table(seqtab,file=paste0(folder,"ForwardReads.txt"),sep="\t")
}

