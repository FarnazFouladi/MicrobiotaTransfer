#Modded from https://benjjneb.github.io/dada2/tutorial.html
#Farnaz Fouladi
#08-18-2019

#Merging sequence variants from different runs, removing chimeras, and assigning taxanomy

library(dada2)
# Merge multiple runs 
st1 <- readRDS("/users/ffouladi/anorexiaMouseTransfer/dada2result/AN34ForwardsReads.rds")
st2 <- readRDS("/users/ffouladi/anorexiaMouseTransfer/dada2result/AN40ForwardsReads.rds")
st3 <- readRDS("/users/ffouladi/anorexiaMouseTransfer/dada2result/AN703ForwardsReads.rds")
st4 <- readRDS("/users/ffouladi/anorexiaMouseTransfer/dada2result/AN81ForwardsReads.rds")
seqtab <- mergeSequenceTables(st1, st2, st3,st4)

#Removing chimeras
seqtab <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
setwd ("/users/ffouladi/anorexiaMouseTransfer/dada2result")
saveRDS(seqtab, "/users/ffouladi/anorexiaMouseTransfer/dada2result/ForwardsReads.rds")
write.table(seqtab,file="ForwardReads.txt",sep="\t")

#Assigning taxanomy 
tax <- assignTaxonomy(seqtab, "/users/ffouladi/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
tax<-addSpecies(tax, "/users/ffouladi/silva_species_assignment_v128.fa.gz")

saveRDS(tax, "/users/ffouladi/anorexiaMouseTransfer/dada2result/taxForwardReads.rds")
write.table(tax,file="taxForwardReads.txt",sep="\t")


