############################################
# Which Reps:
whichreps <- 1:10
howmanyreps <- length(whichreps)
typeeffectsize <- "Unfixed"
# r2 threshold:
threshold_r2 <- 0.8
############################################
# 2018-03-20 [FINISHED]
# Aim: to determine TP eGenes, and to calcualte FDR/TPR of different methods;
#      try differnt LD r2 threshold when determining TP eGenes.

library(foreach)
library(doMC); registerDoMC(cores = 20)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")

# Create a folder
mainfolder <- "~/Simulation_eQTL_replicates/scenario_files/MultipleTestingCorrection_comparison/FDR_TPR_UnfixedEffectSize_Rep10/"


############################################

# An eGene is a TP only if at least one of its eSNPs were in high LD with the causal eSNP.

# Caculate FDR and TPR based on all discoveries across 100 simulations
# Return number of TPs or significant eGenes for each simulation scenario
results_10sim <- foreach(IdxRep = rep(whichreps, each = 36), SS = rep(rep(c(5000,2000,1000,500,200,100), each = 6), howmanyreps), MAF = rep(c(50,25,10,5,1,0.5), 6*howmanyreps), .combine = rbind) %dopar% {
  
  # Scenario directory
  setwd(paste0("~/Simulation_eQTL_replicates/", whichfolder(IdxRep), "/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
  
  # Read signficant associations
  eAsso <- read.table("cis_eAssociations_modelLinear_CP05.txt.gz", header = T)
  
  # Return results of each method
  ddreturn <- foreach(mm = 4:(ncol(eAsso)-1), .combine = rbind) %do% {
    sigeAsso <- eAsso[complete.cases(eAsso[,mm]),]
    thismet <- data.frame(met = colnames(eAsso)[mm],
                          nasso = nrow(sigeAsso),
                          TP = length(which(sigeAsso$R2_causaleSNP >= threshold_r2)),
                          ngene = length(unique(sigeAsso$gene)),
                          TPgene = length(unique(sigeAsso$gene[which(sigeAsso$R2_causaleSNP >= threshold_r2)]))  )
  }
  
  # Scenario parameters
  ddreturn$Rep <- IdxRep
  ddreturn$SS <- SS
  ddreturn$MAF <- MAF
  
  return(ddreturn)
}

# Replace "." in method names by "_"
results_10sim$met <- sub("[.]", "_", results_10sim$met)
# Save the results
setwd(mydir)
results_10sim <- results_10sim[,c("Rep","SS","MAF","met","nasso","TP","ngene","TPgene")]
write.table(results_10sim, paste0("Number_sig_TP_each_simulation_method_r2threshold_", threshold_r2, "_" ,typeeffectsize, "EffectSize_Rep10.txt"), quote = F, sep = "\t", row.names = F)

# Calcualte sum across simulations
# Unique scenarios and methods:
uniqscen <- unique(results_10sim[,c("SS","MAF","met")])
uniqscen_data <- foreach(ii = 1:nrow(uniqscen), .combine = rbind) %do% {
  # A specific scenario and method
  dd <- results_10sim[which(results_10sim$SS == uniqscen$SS[ii] & results_10sim$MAF == uniqscen$MAF[ii] & results_10sim$met == uniqscen$met[ii]),]
  # Return sum of the last four columns
  returndd <- uniqscen[ii,]
  dd <- dd[,-1:-4]
  returndd <- cbind(returndd, t(apply(dd, 2, sum)))
  # Number of simulations
  returndd$n_sims <- nrow(dd)
  return(returndd)
}

# Calculate FDR
uniqscen_data$FDR <- (uniqscen_data$ngene - uniqscen_data$TPgene) / uniqscen_data$ngene
# TPR
uniqscen_data$TPR <- uniqscen_data$TPgene / (uniqscen_data$n_sims * 200)

# Save
write.table(uniqscen_data, paste0("FDR_TPR_all_methods_r2threshold_",threshold_r2,"_",typeeffectsize, "EffectSize_Rep10.txt"), quote = F, sep = "\t", row.names = F)




