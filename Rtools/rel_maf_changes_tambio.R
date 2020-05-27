match.rows <- function(x, y, all=F) 
{
  if (all) 
  {
    if(ncol(x) != ncol(y)) return(rep(F, nrow(x)))
    ret <- list()
    for(x.row in 1:nrow(x)) {
      x.row.matches <- vector()
      for(y.row in 1:nrow(y)) {
        if(all(x[x.row,] == y[y.row,])) {
          x.row.matches <- append(x.row.matches, y.row)
        }
      }
      ret[[x.row]] <- x.row.matches
    }
    return(ret)
  }
  
  else {
    x <- apply(x, 1, paste, collapse="cOlLaPsE")
    y <- apply(y, 1, paste, collapse="cOlLaPsE")
    match(x,y)
  }
}

double.grep.idx <- function(pattern1, pattern2, search.from) 
{
  which(grepl(search.from, pattern=pattern1) & grepl(search.from, pattern=pattern2))
}

library(parallel)
setwd("/home/annalam/datasets/ozm-054")

# DATA
patients <- read.delim2("clinical/patient_arms.tsv", header=F, stringsAsFactors=FALSE)[,1]
patients <- gsub(patients, pattern="-", replacement=".")
data <- read.delim("panel/mutations/somatic.vcf", stringsAsFactors=FALSE)
blacklist <- read.delim("panel/mutations/blacklist.tsv", header=FALSE)
ctdna_fractions <- read.delim("panel/ctdna_fractions.tsv", header=FALSE, stringsAsFactors=FALSE)
ctdna_fractions[,1] <- gsub(ctdna_fractions[,1], pattern="-", replacement=".")
TIMEPOINTS <- c("Baseline", "Progression1", "Progression2")

# MIN CTDNA FRACTION FOR A SAMPLE TO BE ANALYZED
CTDNA_THRESHOLD = 5
# NUMBER OF SIMULATIONS. CORRESPONDS TO P-VALUE AS: 
# P-VAL = MIN(NUMBER OF SIMULATED RELATIVE MAFS SMALLER/LARGER THAT OBSERVED RELATIVE MAF) / N_SIMULATIONS
N_SIMULATIONS <- 1000000
# NUMBER OF CORES TO USE
N_CORES = 8

# PRE-CALCULATED INDICES
too.low.ctfr.idx <- match(ctdna_fractions[ctdna_fractions[,2] < CTDNA_THRESHOLD, 1], colnames(data))
not.blacklisted.idx <- seq(nrow(data))[-na.omit(match.rows(blacklist[,1:4], data[,1:4]))]
called.muts.idx <- apply(data, 2, function(sample) {
  called.mutations <- which(sapply(sample, grepl, pattern="\\*"))
  Filter(function(x) x %in% not.blacklisted.idx, called.mutations)
})

results <- c("PATIENT", "TIMEPOINT 1", "TIMEPOINT 2", "MUTATION A", "MUTATION B", "COUNTS 1A", "COUNTS 2A", "COUNTS 1B", "COUNTS 2B", "PVAL")
for (patient in patients) {
  print(patient)
  p.idx = sapply(TIMEPOINTS, function(tp) double.grep.idx(patient, tp, colnames(data)))
  
  for (tp in 1:length(p.idx)-1) {
    tp1 = unlist(p.idx[tp])
    tp2 = unlist(p.idx[tp+1])
    
    if(length(tp1)==0 | length(tp2)==0) next
    if(tp1 %in% too.low.ctfr.idx | tp2 %in% too.low.ctfr.idx) next
    
    p.called.muts.idx = unique(unlist(c(called.muts.idx[tp1], called.muts.idx[tp2])))
    if (length(p.called.muts.idx) <= 1) next
    
    combinations = combn(p.called.muts.idx, 2)
    c <- mclapply(1:ncol(combinations), mc.cores=N_CORES, FUN=function(combn) {
      muta <- combinations[1, combn]
      mutb <- combinations[2, combn]
      
      tp1a <- unlist(strsplit(data[muta, tp1], ":"))
      tp2a <- unlist(strsplit(data[muta, tp2], ":"))
      tp1b <- unlist(strsplit(data[mutb, tp1], ":"))
      tp2b <- unlist(strsplit(data[mutb, tp2], ":"))
      
      tot1a <- as.numeric(tp1a[2]); tot1b <- as.numeric(tp1b[2])
      tot2a <- as.numeric(tp2a[2]); tot2b <- as.numeric(tp2b[2])
      mut1a <- as.numeric(tp1a[1]); mut1b <- as.numeric(tp1b[1])
      mut2a <- as.numeric(tp2a[1]); mut2b <- as.numeric(tp2b[1])
      
      if (tot1a ==  0 | tot2a ==  0 | tot1b ==  0 | tot2b ==  0) NA
      else {
        maf1a <- mut1a / tot1a
        maf2a <- mut2a / tot2a
        maf1b <- mut1b / tot1b
        maf2b <- mut2b / tot2b
        
        if (maf1a < 0.03 & maf2a < 0.03) NA
        else if (maf1b < 0.03 & maf2b < 0.03) NA
        else if (maf1a < 0.03 & maf1b < 0.03) NA
        else if (maf2a < 0.03 & maf2b < 0.03) NA
        else {
          
          nonzero_maf1a <- max(maf1a, 1/tot1a)
          nonzero_maf1b <- max(maf1b, 1/tot1b)
          
          sim_maf1a <- rbinom(N_SIMULATIONS, tot1a, nonzero_maf1a) / tot1a
          sim_maf1b <- rbinom(N_SIMULATIONS, tot1b, nonzero_maf1b) / tot1b
          
          # Null hypothesis: relative MAF does not change.
          # For simulation, we adjust later timepoint MAFs so the
          # null hypothesis is true, while conserving average MAF.
          average = (maf2a + maf2b) / 2
          sim_rel_maf1 = sim_maf1b / sim_maf1a
          sim_rel_maf1 <- sim_rel_maf1[!is.na(sim_rel_maf1)]
          
          expected_maf2a <- sapply(sim_rel_maf1, function(x){ min(1, 2*average/(1+x)) })
          expected_maf2b <- sapply(sim_rel_maf1, function(x){ min(1, 2*average/(1+1/x)) })
          
          sim_maf2a <- sapply(expected_maf2a, function(x) rbinom(1, tot2a, x)) / tot2a
          sim_maf2b <- sapply(expected_maf2b, function(x) rbinom(1, tot2b, x)) / tot2b
          
          sim_rel_maf2 <- sim_maf2a / sim_maf2b
          sim_rel_maf2 <- sim_rel_maf2[!is.na(sim_rel_maf2)]
          obs_rel_maf <- maf2a / maf2b
          
          pval = min(sum(sim_rel_maf2>obs_rel_maf), 
                     sum(sim_rel_maf2<obs_rel_maf)) / length(sim_rel_maf2)
          
          if ( pval < 0.001 ) {
            A <- paste(data[muta,"GENE"], muta, sep=";")
            B <- paste(data[mutb,"GENE"], mutb, sep=";")
            
            A1 <- paste0(mut1a, "/", tot1a, "=", round(maf1a,3)); A2 <- paste0(mut2a, "/", tot2a, "=", round(maf2a,3))
            B1 <- paste0(mut1b, "/", tot1b, "=", round(maf1b,3)); B2 <- paste0(mut2b, "/", tot2b, "=", round(maf2b,3))
            
            c(patient, TIMEPOINTS[tp], TIMEPOINTS[tp+1], A, B, A1, A2, B1, B2, pval)
          } else NA
        }
      }
    })
    c <- c[lengths(c) == 10]
    results <- rbind(results, do.call(rbind, c))
  }
}

write.table(results, "~/tmp/sig_mutation_pairs.tsv",
            quote=F, col.names=F, row.names=F, sep="\t")

