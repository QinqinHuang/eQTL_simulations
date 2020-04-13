############################################
# Which Replicates:
whichreps <- 11:100
howmanyreps <- length(whichreps)
############################################

library(foreach)
library(doMC); registerDoMC(cores = 48)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")


#----------------------------------------------------------------------
# Run FastQTL - Beta approximate scheme
#----------------------------------------------------------------------
# All scenarios
nothing <- foreach(SS = rep(c(5000,2000,1000,500,200,100), each = 6*howmanyreps), IdxRep = rep(rep(whichreps, each = 6), 6), MAF = rep(c(50,25,10,5,1,0.5), 6*howmanyreps) ) %dopar% {
  
  # Directory of current scenario
  setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
  
  # Check whether FastQTL has finished the job for this scenario
  flag <- F
  if(!file.exists("Permutation_1k.log")){flag <- T} else {
    lastline <- unlist(strsplit(readLines(pipe("tail -n 1 Permutation_1k.log")), split = " "))
    if(!"Running" %in% lastline) {flag <- T}
  }
  
  # If FastQTL hasn't finished
  if(flag) {
    try( system("fastQTL.static --vcf ../Genotype/filtered.vcf.gz --bed gene_expression.bed.gz --permute 1000 --window 1000000 --out Permutation_1k.txt.gz --log Permutation_1k.log --chunk 1 1") )
  }
  
  return(NULL)
}



#----------------------------------------------------------------------
# Run FastQTL - Adaptive scheme for the first 10 replicates
#----------------------------------------------------------------------
# Which Replicates:
whichreps <- 1:10
howmanyreps <- length(whichreps)

# All scenarios
nothing <- foreach(SS = rep(c(5000,2000,1000,500,200,100), each = 6*howmanyreps), IdxRep = rep(rep(whichreps, each = 6), 6), MAF = rep(c(50,25,10,5,1,0.5), 6*howmanyreps) ) %dopar% {
  
  # Directory of current scenario
  setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
  
  # Check whether FastQTL has finished the job for this scenario
  flag <- F
  if(!file.exists("Permutation_100-10k.log")){flag <- T} else {
    lastline <- unlist(strsplit(readLines(pipe("tail -n 1 Permutation_100-10k.log")), split = " "))
    if(!"Running" %in% lastline) {flag <- T}
  }
  
  # If FastQTL hasn't finished
  if(flag) {
    try( system("fastQTL.static --vcf ../Genotype/filtered.vcf.gz --bed gene_expression.bed.gz --permute 100 10000 --window 1000000 --out Permutation_100-10k.txt.gz --log Permutation_100-10k.log --chunk 1 1") )
  }
  
  return(NULL)
}





#----------------------------------------------------------------------
# Run FastQTL after excluding some genes with sample size is 100
#----------------------------------------------------------------------
# Some genes failed when sample size is 100.
SS <- 100

########## Exclude some genes ##########
nothing <- foreach(IdxRep = rep(11:100, each = 6), MAF = rep(c(50,25,10,5,1,0.5), 90)) %do% {
  # Directory of current scenario
  setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
  lastline <- unlist(strsplit(readLines(pipe("tail -n 1 Permutation_1k.log")), split = " "))
  if(!"Running" %in% lastline) {
    # Read log
    genetoexc <- readLines(pipe("cat Permutation_1k.log|grep gene|tail -n 1"))
    cat(IdxRep, SS, MAF, genetoexc, "\t")
    genename <- unlist(strsplit(genetoexc, split = "[[]"))[2]
    genename <- unlist(strsplit(genename, split = "[]]"))
    cat(genename, "\n")
    #write(genename, "Permutation_FastQTL_exclude_genes.txt", append = T)
  }
  return(NULL)
}

# All simulations
nothing <- foreach(IdxRep = rep(11:100, each = 6), MAF = rep(c(50,25,10,5,1,0.5), 90)) %dopar% {
  # Directory of current scenario
  setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
  flag <- F
  if(!file.exists("Permutation_1k.log")){flag <- T} else {
    lastline <- unlist(strsplit(readLines(pipe("tail -n 1 Permutation_1k.log")), split = " "))
    if(!"Running" %in% lastline) {flag <- T}
  }
  
  # If FastQTL hasn't finished
  if(flag) {
    if(file.exists("Permutation_FastQTL_exclude_genes.txt")) {
      cat(IdxRep, SS, MAF, readLines("Permutation_FastQTL_exclude_genes.txt"), "\n")
      system("fastQTL.static --vcf ../Genotype/filtered.vcf.gz --bed gene_expression.bed.gz --permute 1000 --window 1000000 --out Permutation_1k.txt.gz --log Permutation_1k.log --chunk 1 1 --exclude-phenotypes Permutation_FastQTL_exclude_genes.txt")
      
      # Or adaptive scheme
      #system("fastQTL.static --vcf ../Genotype/filtered.vcf.gz --bed gene_expression.bed.gz --permute 100 10000 --window 1000000 --out Permutation_100-10k.txt.gz --log Permutation_100-10k.log --chunk 1 1 --exclude-phenotypes FastQTL_exclude_gene_adaptive.txt")
    }
  }
  return(NULL)
}








