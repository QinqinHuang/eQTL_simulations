############################################
############# Basic parameters #############
# Which Reps:
whichreps <- 1:100
howmanyreps <- length(whichreps)
typeeffectsize <- "Cons"
localcormet <- "eigenMT"
globalcormet <- "BH"
# r2 threshold:
threshold_r2 <- 0.8
############################################
# 2017-10-24 [FINISHED]
# Aim: to get top SNP for each gene, LD r2 with the causal eSNP.
# Sort by statistics to make sure the most significant SNPs have minimum pvalues and largest statistics.
# Did not keep all SNPs that are in perfect LD.

library(foreach)
library(doMC); registerDoMC(cores = 20)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")

# Create a new folder - outputs, log files, Rscript
mainfolder <- "~/Simulation_eQTL_replicates/scenario_files/Top_SNP_causal_SNP/"
mydir <- paste0(mainfolder, "Top_SNP_all_genes_causalSNP_info_",localcormet,"_",globalcormet,"_",typeeffectsize,"EffectSize/")
if(!file.exists(mydir)) {dir.create(mydir)}
setwd(mydir)
if(!file.exists("output_files")) {dir.create("output_files")}

# Logfile 
namelogfile <- paste0("LogFile_top_SNP_rep", min(whichreps), "_", max(whichreps), "_", typeeffectsize,"EffectSize.txt")
write(paste0(" ** Date: ",Sys.time()), file = namelogfile, append = T)
write(paste0(" ** Hierarchical correction method: ",localcormet,"-",globalcormet," \n"), file = namelogfile, append = T)

############################################


# Generate a file for each scenario, top SNP of all genes.
nothing <- foreach( IdxRep = rep(whichreps, each = 36*4), SS = rep(rep(c(5000,2000,1000,500,200,100), each = 6*4), howmanyreps), MAF = rep(rep(c(50,25,10,5,1,0.5), each = 4), 6*howmanyreps), effectsize = rep(c(0.25, 0.5, 1, 1.5), 36*howmanyreps), .combine = rbind ) %dopar% {
  
  # Read MatrixEQTL output
  setwd(paste0("~/Simulation_eQTL_replicates/", whichfolder(IdxRep), "/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF, "/Constant_effect_size/effectsize", effectsize))
  me <- readRDS("MatrixEQTL.output.rds")
  # all tests (keep t-statistics)
  cis <- me$cis$eqtls[,c("snps","gene","statistic","pvalue","beta")]
  cis$gene <- as.character(cis$gene)
  colnames(cis)[ncol(cis)] <- "nominal_beta"
  
  # Sort by absolute statistics
  cis <- cis[order(abs(cis$statistic), decreasing = T),]
  # Top SNP with the min pvalue for each gene
  topSNP <- cis[!duplicated(cis$gene),]
  rownames(topSNP) <- 1:nrow(topSNP)
  
  # LD r2 between top SNP and the causal eSNP for significant eGenes
  eAsso <- read.table("cis_eAssociations_modelLinear_CP05.txt.gz", header = T)[,c("snps","gene","R2_causaleSNP",paste0(localcormet,".",globalcormet))]
  eAsso <- eAsso[!is.na(eAsso[,4]),]
  topSNP <- merge(topSNP, eAsso, by = c("gene","snps"), all.x = T)
  
  # Significant (including False Positives)
  topSNP[,ncol(topSNP)] <- !is.na(topSNP[,ncol(topSNP)])
  colnames(topSNP)[ncol(topSNP)] <- paste0("sig_",localcormet,globalcormet)
  
  # True Positives
  TPlist <- unique(eAsso[which(eAsso$R2_causaleSNP>=threshold_r2),"gene"])
  topSNP$TP <- FALSE
  topSNP$TP[which(topSNP$gene %in% TPlist)] <- TRUE
  colnames(topSNP)[ncol(topSNP)] <- paste0("TP_",localcormet,globalcormet)
  
  # Causal eSNP
  eQTLsnp <- read.table("../../True_eQTLs.txt", header = T)
  topSNP <- merge(topSNP, eQTLsnp[,c("eGene","SNP")], by.x = "gene", by.y = "eGene", all.x = T)
  colnames(topSNP)[ncol(topSNP)] <- "causal_snp"
  
  setwd(mydir); setwd("output_files")
  output_filename <- paste0("Rep",IdxRep,"SS",SS,"MAF",MAF,"effectsize",effectsize,"_",typeeffectsize,"_Top_SNP_all_genes_causalSNP_info_",localcormet,globalcormet,".txt")
  write.table(topSNP, output_filename, row.names = F, quote = F, sep = "\t")
  
  return(NULL)
}


#----- 2017.10.17 -----
# R script is under: "~/Simulation_eQTL_replicates/scenario_files/Top_SNP_causal_SNP/Top_SNP_all_genes_causalSNP_info_eigenMT_BH_ConsEffectSize/Get_topSNP_sim100.R"
# Outputs are under ./output_files/

# Running time 20 cores, 100 reps, 4.46 hr
# A sceond time: 6.74 hr using 20 cores, 100 simulations.
#real	404m29.354s
#user	760m28.864s
#sys	258m15.916s



