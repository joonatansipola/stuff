#!/usr/bin/env Rscript
d<-scan("stdin", quiet=TRUE)
cat(paste0("min ", min(d)), 
    paste0("max ", max(d)), 
    paste0("median ", median(d)), 
    paste0("mean ", mean(d)), 
    sep="\n")
