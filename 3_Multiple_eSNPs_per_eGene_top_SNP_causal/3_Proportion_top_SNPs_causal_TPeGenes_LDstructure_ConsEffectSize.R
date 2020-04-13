############################################
############# Basic parameters #############
# Which Reps:
whichreps <- 1:100
howmanyreps <- length(whichreps)
typeeffectsize <- "Cons"
localcormet <- "eigenMT"
globalcormet <- "BH"
############################################
# 2017-10-24 [FINISHED]
# Aim: to calcluate the proportion of top SNPs that are causal;
#      to confirm how LD structrue around the causal eSNP affects the identification.
# Focus on only TP eGenes because we know the causal eSNPs.

library(foreach)
library(doMC); registerDoMC(cores = 20)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")

# Find the folder
mainfolder <- "~/Simulation_eQTL_replicates/scenario_files/Top_SNP_causal_SNP/"
mydir <- paste0(mainfolder, "Top_SNP_all_genes_causalSNP_info_",localcormet,"_",globalcormet,"_",typeeffectsize,"EffectSize/")

############################################


# Return: (1) calculate the proportion of top SNPs that are causal
#         (2) LD info around causal eSNPs of TP eGenes

setwd(mydir)

num_top_causal <- foreach( IdxRep = rep(whichreps, each = 36*4), SS = rep(rep(c(5000,2000,1000,500,200,100), each = 6*4), howmanyreps), MAF = rep(rep(c(50,25,10,5,1,0.5), each = 4), 6*howmanyreps), effectsize = rep(c(0.25, 0.5, 1, 1.5), 36*howmanyreps), .combine = rbind ) %dopar% {
  # Read the info about top SNPs
  output_filename <- paste0("./output_files/Rep",IdxRep,"SS",SS,"MAF",MAF,"effectsize",effectsize,"_",typeeffectsize,"_Top_SNP_all_genes_causalSNP_info_",localcormet,globalcormet,".txt")
  topSNP <- read.table(output_filename, header = T)
  # TP eGenes
  topSNP_TP <- topSNP[which(topSNP[,paste0("TP_",localcormet,globalcormet)]),]
  
  # Read LD structure around causal eSNPs
  LDinfofilename <- paste0("~/Simulation_eQTL_replicates/", whichfolder(IdxRep), "/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF, "/True_eQTLs_with_LD_info_around_causaleSNP.txt")
  # Same with unfixed effect size scenarios but TPs are different.
  causaleSNP <- read.table(LDinfofilename, header = T)
  # Focus on TP eGenes
  topSNP_TP <- merge(topSNP_TP, causaleSNP[,c("eGene","num_SNP05_2Mb")], by.x = "gene", by.y = "eGene", all.x = T)
  
  # Return table
  returndd <- data.frame(Rep = IdxRep, SS = SS, MAF = MAF, effectsize = effectsize)
  # (1) number of significant eGenes
  returndd$nsig <- sum(topSNP[,paste0("sig_",localcormet,globalcormet)])
  # (2) number of TPs
  returndd$nTP <- nrow(topSNP_TP)
  # (3) number of top SNPs that are causal, e.g. R2 with causal eSNP is 1
  r2 <- topSNP_TP$R2_causaleSNP    # NA if any NA exists
  returndd$ncausal <- sum(r2[complete.cases(r2)] == 1)
  # (4) sum number of SNPs in 2Mb window in LD with the causal eSNP
  returndd$sum_num_SNP05_2Mb <- sum(topSNP_TP$num_SNP05_2Mb)
  
  return(returndd)
}

# Sum across simulations and calculate proportions
sumrep_num_top_causal <- calculate_sum_ignoreNAs(num_top_causal, scenario_parameters = c("SS","MAF","effectsize"))
# Calculate the proportion of top SNPs that are in perfect LD with the causal SNP
sumrep_num_top_causal$proportion <- sumrep_num_top_causal$ncausal / sumrep_num_top_causal$nTP
# Calcualte the average LD structure
sumrep_num_top_causal$ave_num_SNP05_2Mb <- sumrep_num_top_causal$sum_num_SNP05_2Mb / sumrep_num_top_causal$nTP
# Save
write.table(sumrep_num_top_causal, paste0("Proportion_top_SNP_causal_TPeGenes_LDstructure_",localcormet,"_",globalcormet,"_",typeeffectsize,"EffectSize.txt"), quote = F, sep = "\t", row.names = F)


# Return: LD r2 of top SNPs of TP eGenes
LDr2_allscen <- foreach( IdxRep = rep(whichreps, each = 36*4), SS = rep(rep(c(5000,2000,1000,500,200,100), each = 6*4), howmanyreps), MAF = rep(rep(c(50,25,10,5,1,0.5), each = 4), 6*howmanyreps), effectsize = rep(c(0.25, 0.5, 1, 1.5), 36*howmanyreps), .combine = rbind ) %dopar% {
  # Read the info about top SNPs
  output_filename <- paste0("./output_files/Rep",IdxRep,"SS",SS,"MAF",MAF,"effectsize",effectsize,"_",typeeffectsize,"_Top_SNP_all_genes_causalSNP_info_",localcormet,globalcormet,".txt")
  topSNP <- read.table(output_filename, header = T)
  # TP eGenes
  topSNP_TP <- topSNP[which(topSNP[,paste0("TP_",localcormet,globalcormet)]),]
  # Return table
  if(nrow(topSNP_TP) > 0) {
    returndd <- data.frame(Rep = IdxRep, SS = SS, MAF = MAF, effectsize = effectsize,
                           LDr2 = topSNP_TP$R2_causaleSNP)
    return(returndd)
  } else {   # if no TP eGene, don't return LD r2
    return(NULL)
  }
}

# Save LD r2 of top SNPs of TP eGenes
write.table(LDr2_allscen, paste0("LDr2_list_top_SNP_causal_TPeGenes_",localcormet,"_",globalcormet,"_",typeeffectsize,"EffectSize.txt"), quote = F, sep = "\t", row.names = F)


