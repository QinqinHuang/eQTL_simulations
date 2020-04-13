############################################
# Which Replicates:
whichreps <- 11:100
howmanyreps <- length(whichreps)
############################################

library(foreach)
library(doMC); registerDoMC(cores = 25)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")


#----------------------------------------------------------------------
# Run eigenMT
#----------------------------------------------------------------------
# Six sample sizes
nothing <- foreach(SS = rep(c(5000,2000,1000,500,200,100), each = howmanyreps), IdxRep = rep(whichreps, 6), .inorder = F) %dopar% {
  # Genotype directory
  setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/Genotype"))
  
  # Check whether we have eigenMT output
  if(file.exists("eigenMT_output")) {return(NULL)}
  
  # Run eigenMT
  system("time python ~/bin/eigenMT/eigenMT.py --CHROM 22 --QTL ../SS*maf50/MatrixEQTL_all_cis_tests_pval.txt --GEN SNP.txt --GENPOS snpsloc.txt --PHEPOS ~/Simulation_eQTL_replicates/hg19_gene_annotation/geneloc_for_ME.txt --OUT eigenMT_output")
  system("rm ../SS*maf50/MatrixEQTL_all_cis_tests_pval.txt")
  
  return(NULL) 
}



