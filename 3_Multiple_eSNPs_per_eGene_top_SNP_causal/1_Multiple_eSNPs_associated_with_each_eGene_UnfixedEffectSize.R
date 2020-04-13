############################################
############# Basic parameters #############
# Which Reps:
whichreps <- 1:100
howmanyreps <- length(whichreps)
typeeffectsize <- "Unfixed"
localcormet <- "eigenMT"
globalcormet <- "BH"
# r2 threshold:
threshold_r2 <- 0.8
############################################
# 2017-10-24 [FINISHED]
# Aim: to calculate the average number of significant eSNPs per TP eGene,
# as well as per significant eGenes (including FPs).

library(foreach)
library(doMC); registerDoMC(cores = 20)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")

# Create a new folder - outputs, log files, Rscript
mainfolder <- "~/Simulation_eQTL_replicates/scenario_files/Multiple_eSNPs_per_TPeGene/"
mydir <- paste0(mainfolder, "Multiple_eSNP_",localcormet,"_",globalcormet,"_",typeeffectsize,"EffectSize/")
if(!file.exists(mydir)) {dir.create(mydir)}
setwd(mydir)

############################################

# Return (1) the number of TP eGenes and their eSNPs
#        (2) the number of significant eGenes and their eSNPs
# so that we can calculate the sum across all simulations and calculate the average

num_eAsso_allscenarios <- foreach(IdxRep = rep(whichreps, each = 36), SS = rep(rep(c(5000,2000,1000,500,200,100), each = 6), howmanyreps), MAF = rep(c(50,25,10,5,1,0.5), 6*howmanyreps), .combine = rbind) %dopar% {
  # Directory of this scenario
  setwd(paste0("~/Simulation_eQTL_replicates/", whichfolder(IdxRep), "/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
  # Read significant associations
  eAsso <- read.table("cis_eAssociations_modelLinear_CP05.txt.gz", header = T, stringsAsFactors = F)
  
  # Significant associations by this method
  eAsso <- eAsso[,c("snps","gene","R2_causaleSNP",paste0(localcormet,".",globalcormet))]
  eAsso <- eAsso[which(!is.na(eAsso[,ncol(eAsso)])),]
  
  # TP eGene list (requires eSNP in LD with the causal eSNP)
  TPgenes <- unique(eAsso[which(eAsso$R2_causaleSNP >= threshold_r2),"gene"])
  
  # Return table, scenario parameters
  ddreturn <- data.frame(Rep = IdxRep, SS = SS, MAF = MAF)
  # Number of TP eGenes
  ddreturn$TPeGene <- length(TPgenes)
  # Number of their eSNPs
  ddreturn$eSNP_TPeGene <- sum(eAsso$gene %in% TPgenes)
  # Number of eGenes
  ddreturn$eGene <- length(unique(eAsso$gene))
  # Number of their eSNPs
  ddreturn$eSNP_eGene <- nrow(eAsso)
  
  return(ddreturn)
}

# Calculate the average across 100 simulations
# sum across 100 simulations
sum_data <- calculate_sum_ignoreNAs(num_eAsso_allscenarios, c("SS","MAF"))
# average number of eSNPs
sum_data$n_SNP_TPeGene <- sum_data$eSNP_TPeGene / sum_data$TPeGene
sum_data$n_SNP_sigeGene <- sum_data$eSNP_eGene / sum_data$eGene
ave_nSNP <- sum_data[,c("SS","MAF","n_SNP_TPeGene","n_SNP_sigeGene")]

# Save results
setwd(mydir)
# Add the multiple testing correction method name
colnames(num_eAsso_allscenarios)[4:7] <- paste0(colnames(num_eAsso_allscenarios)[5:8],"_",localcormet,"_",globalcormet)
colnames(ave_nSNP)[3:4] <- paste0(c("n_SNP_TPeGene","n_SNP_sigeGene"),"_",localcormet,"_",globalcormet)
write.table(num_eAsso_allscenarios, paste0("Number_eSNPs_TPgene_siggene_each_scenario_100reps_",localcormet,"_",globalcormet,"_",typeeffectsize,"EffectSize.txt"), row.names = F, quote = F, sep = "\t")
write.table(ave_nSNP, paste0("Average_number_eSNPs_per_TPgeneORsiggene_100reps_",localcormet,"_",globalcormet,"_",typeeffectsize,"EffectSize.txt"), row.names = F, quote = F, sep = "\t")







