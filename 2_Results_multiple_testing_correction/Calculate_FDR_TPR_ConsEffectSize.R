############################################
############# Basic parameters #############
# Which Reps:
whichreps <- 1:100
howmanyreps <- length(whichreps)
typeeffectsize <- "Cons"
# r2 threshold:
threshold_r2 <- 0.8
############################################
# 2017-10-24 ()
# Aim: to determine TP eGenes, and to calcualte FDR/TPR of different methods;
#      try differnt LD r2 threshold when determining TP eGenes.

library(foreach)
library(doMC); registerDoMC(cores = 20)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")

# Create a folder
mainfolder <- "~/Simulation_eQTL_replicates/scenario_files/MultipleTestingCorrection_comparison/"
mydir <- paste0(mainfolder, "FDR_TPR_",typeeffectsize,"EffectSize/")
if(!file.exists(mydir)) {dir.create(mydir)}

############################################

# An eGene is a TP only if at least one of its eSNPs were in high LD with the causal eSNP.

# Caculate FDR and TPR based on all discoveries across 100 simulations
# Return number of TPs or significant eGenes for each simulation scenario
results_100sim <- foreach( IdxRep = rep(whichreps, each = 36*4), SS = rep(rep(c(5000,2000,1000,500,200,100), each = 6*4), howmanyreps), MAF = rep(rep(c(50,25,10,5,1,0.5), each = 4), 6*howmanyreps), effectsize = rep(c(0.25, 0.5, 1, 1.5), 36*howmanyreps), .combine = rbind ) %dopar% {
  
  # Scenario directory
  setwd(paste0("~/Simulation_eQTL_replicates/", whichfolder(IdxRep), "/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF, "/Constant_effect_size/effectsize", effectsize))
  
  # Read signficant associations
  eAsso <- read.table("cis_eAssociations_modelLinear_CP05.txt.gz", header = T)
  
  # Return results of each method (didn't run BPerm1k in all 100 simulations)
  ddreturn <- foreach(mm = paste0(c("eigenMT","Bonferroni"),".BH"), .combine = rbind) %do% {
    sigeAsso <- eAsso[complete.cases(eAsso[,mm]),]
    thismet <- data.frame(Method = mm,
                          nasso = nrow(sigeAsso),
                          TP = length(which(sigeAsso$R2_causaleSNP >= threshold_r2)),
                          ngene = length(unique(sigeAsso$gene)),
                          TPgene = length(unique(sigeAsso$gene[which(sigeAsso$R2_causaleSNP >= threshold_r2)]))  )
  }
  
  # Scenario parameters
  ddreturn$Rep <- IdxRep
  ddreturn$SS <- SS
  ddreturn$MAF <- MAF
  ddreturn$effectsize <- effectsize
  
  return(ddreturn)
}

# Replace "." in method names by "_"
results_100sim$Method <- sub("[.]", "_", results_100sim$Method)
# Save the results
setwd(mydir)
results_100sim <- results_100sim[,c("Rep","SS","MAF","effectsize","Method","nasso","TP","ngene","TPgene")]
write.table(results_100sim, paste0("Number_sig_TP_each_simulation_method_r2threshold_", threshold_r2, "_" ,typeeffectsize, "EffectSize.txt"), quote = F, sep = "\t", row.names = F)

# Calcualte sum across simulations
# Unique scenarios and methods:
uniqscen <- unique(results_100sim[,c("SS","MAF","effectsize","Method")])
uniqscen_data <- foreach(ii = 1:nrow(uniqscen), .combine = rbind) %dopar% {
  # A specific scenario and method
  dd <- results_100sim[which(results_100sim$SS == uniqscen$SS[ii] & results_100sim$MAF == uniqscen$MAF[ii] & results_100sim$effectsize == uniqscen$effectsize[ii] & results_100sim$Method == uniqscen$Method[ii]),]
  # Return sum of the last four columns
  returndd <- uniqscen[ii,]
  dd <- dd[,-1:-5]
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
write.table(uniqscen_data, paste0("FDR_TPR_all_methods_r2threshold_",threshold_r2,"_",typeeffectsize, "EffectSize.txt"), quote = F, sep = "\t", row.names = F)





