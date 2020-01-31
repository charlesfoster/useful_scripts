#!/usr/bin/env Rscript

### Usage: rscript prune_fasta.R
### Requirement: a directory of fasta files, and a text file with the taxa that you want 
### to keep in your alignments ("taxa_to_keep.txt")

library(ape)

keep <- readLines("taxa_to_keep.txt")
files <- list.files(patt="*fasta")
data = lapply(files, function(f) {
wb = read.FASTA(f, type="DNA")})
names(data) <- files

pruned <- lapply(data,function(f){
tmp <- f[names(f) %in% keep]
return(tmp)
})

outnames <- names(pruned)
lapply(outnames,function(f){
fname <- gsub(".fasta","_pruned.fasta", f)
cat("\nPruned: ",fname,"\n")
ofile <- pruned[[grep(f, names(pruned))]]
write.FASTA(file=fname, x=ofile)
})
