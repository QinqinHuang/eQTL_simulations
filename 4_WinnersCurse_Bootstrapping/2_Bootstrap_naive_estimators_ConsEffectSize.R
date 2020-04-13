############################################
############# Basic parameters #############
# Which Reps:
whichreps <- 1:10
howmanyreps <- length(whichreps)
typeeffectsize <- "Cons"
localcormet <- "eigenMT"
globalcormet <- "BH"
############################################
# 2017-10-19 [FINISHED]
# Aim: to compare Bootstrap estimators as well as naive estimator

# Version of BootstrapQTL
versionbootsQ <- "0.7.3"
n_bootstraps <- 200

library(foreach)
library(doMC); registerDoMC(cores = 20)
library(reshape)

# Data under:
mydir <- "~/Simulation_eQTL_replicates/scenario_files/Bootstrap_reestimation_effect_size/BootstrapQTL_v0.7.3_eigenMT-BH_nomthres_200_ConsEffectSize/"
setwd(mydir)

# Log file, save version info
logfilename <- paste0(mydir,"LogFile_Bootstrap_estimator_rep",paste0(range(whichreps),collapse = "-"),".txt")
write(paste0("***** BootstrapQTL Version: ",versionbootsQ, " ***** "), file = logfilename, append = T)
write(paste0("      Date: ",Sys.time()), file = logfilename, append = T)
write(paste0("      Hierarchical correction method: ",localcormet,"-",globalcormet," \n"), file = logfilename, append = T)

# Function to write logs
writelog <- function(newinfo, IdxRep, SS, MAF, effectsize) {
  newline <- paste0("* Rep", IdxRep, ", Sample size=", SS, ", MAF=", MAF, "%, effect size=", effectsize, ": ", newinfo)
  write(newline, file = logfilename, append = T)
}

############################################


# Keep all 200 true eGenes of which true effect size is known.
# For TPs, return all three bootstrap estimatros;
# for false negatives, return only naive estimator.

# Bootstrap results
bootsest <- foreach( IdxRep = rep(whichreps, each = 36*4), SS = rep(rep(c(5000,2000,1000,500,200,100), each = 6*4), howmanyreps), MAF = rep(rep(c(50,25,10,5,1,0.5), each = 4), 6*howmanyreps), effectsize = rep(c(0.25, 0.5, 1, 1.5), 36*howmanyreps), .combine = rbind ) %dopar% {
  # Top SNP for all genes:
  topSNPfilename <- paste0("~/Simulation_eQTL_replicates/scenario_files/Top_SNP_causal_SNP/Top_SNP_all_genes_causalSNP_info_",localcormet,"_",globalcormet,"_",typeeffectsize,"EffectSize/output_files/Rep",IdxRep,"SS",SS,"MAF",MAF,"effectsize",effectsize,"_",typeeffectsize,"_Top_SNP_all_genes_causalSNP_info_",localcormet,globalcormet,".txt")
  topSNP <- read.table(topSNPfilename, header = T)
  
  # True eGenes including TPs and false negatives
  topSNP200true <- topSNP[which(!is.na(topSNP$causal_snp)),]
  
  # 1.Naive estimator for all true eGenes
  returndd <- topSNP200true[,c(paste0("TP_",localcormet,globalcormet),"nominal_beta")]
  # Not identified or TP
  returndd$Estimator <- "not identified"
  returndd$Estimator[which(returndd[,paste0("TP_",localcormet,globalcormet)])] <- "Naive"
  returndd <- returndd[,-1]
  colnames(returndd)[1] <- "estimates"
  # scenario parameters
  returndd$Rep <- IdxRep; returndd$SS <- SS; returndd$MAF <- MAF; returndd$effectsize <- effectsize
  
  # If no TP, only return naive estimator for not significant true eGene.
  n_TP <- sum(topSNP200true[,paste0("TP_",localcormet,globalcormet)])
  if(n_TP == 0) { return(returndd) }
  
  # Read BootstrapQTL output
  setwd(mydir); setwd("output_files")
  bootsfilename <- paste0("Rep",IdxRep,"SS",SS,"MAF",MAF,"effectsize",effectsize,"_BootstrapQTL",n_bootstraps,"_",typeeffectsize,"EffectSize_output.txt")
  if(!file.exists(bootsfilename)) {
    writelog(newinfo = "Error - no BootstrapQTL output.", IdxRep, SS, MAF, effectsize)
    return(NULL)
  }
  bootsout <- read.table(bootsfilename, header = T)
  
  # Check results of BootstrapQTL and my data
  bootsout <- bootsout[order(abs(bootsout$statistic), decreasing = T),]
  bootsout_top <- bootsout[!duplicated(bootsout$eGene),]
  # (i). number of significant eGenes
  if(sum(topSNP[,paste0("sig_",localcormet,globalcormet)]) != nrow(bootsout_top)) {  writelog(newinfo = "Error - BootstrapQTL not the same number of significant eGenes.", IdxRep, SS, MAF, effectsize) }
  # (ii). top SNPs and nominal betas
  if(!identical( formatC(topSNP$nominal_beta[match(bootsout_top$eGene, topSNP$gene)], digits = 5), formatC(bootsout_top$nominal_beta, digits = 5) )) { writelog(newinfo = "Error - nominal betas are not identical.", IdxRep, SS, MAF, effectsize) }
  # (iii). number of TP eGenes
  bootsout_TP <- bootsout_top[which(bootsout_top$eGene %in% topSNP$gene[which(topSNP[,paste0("TP_",localcormet,globalcormet)])]),]
  if(n_TP != nrow(bootsout_TP)) { writelog(newinfo = "Error - BootstrapQTL not the same number of TPs.", IdxRep, SS, MAF, effectsize) }
  
  # 2. Bootstrap estimator for True Positive eGenes
  bootest <- bootsout_TP[,c("shrinkage","outof","weighted")]
  bootest <- melt(bootest, id.vars = NULL)
  colnames(bootest) <- c("Estimator","estimates")
  bootest$Rep <- IdxRep; bootest$SS <- SS; bootest$MAF <- MAF; bootest$effectsize <- effectsize
  
  # Return shrinkage estimators
  returndd <- rbind(bootest, returndd)
  
  return(returndd)
}

setwd(mydir)
write.table(bootsest, paste0("BootstrapQTL_estimators_Naive_not_identified_",localcormet,globalcormet,"_",typeeffectsize,"EffectSize_sim",howmanyreps,".txt"), quote = F, sep = "\t", row.names = F)




