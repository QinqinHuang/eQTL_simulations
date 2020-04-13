############################################
############# Basic parameters #############
# Which Reps:
whichreps <- 1:10
howmanyreps <- length(whichreps)
typeeffectsize <- "Unfixed"
localcormet <- "eigenMT"
globalcormet <- "BH"
############################################
# 2017-10-17 [FINISHED]
# Aim: to test if using eigenMT-BH to correct for multiple testing,
# instead of using the default Bonf-BH

# Install the updated BootstrapQTL
#library(devtools)
#install_github("InouyeLab/BootstrapQTL")

library(BootstrapQTL)
library(foreach)

# Version: (should be v 0.7.3, 2017-10-17 installed)
# Default of this version: same SNP shrinkage estimator, nominal thresholds(bonfbh) in detection group
#                          return all significant eSNPs
versionbootsQ <- sessionInfo()$otherPkgs$BootstrapQTL$Version
cat(" ** Version of BootstrapQTL:", versionbootsQ, "\n")

# Number of bootstrap
n_bootstraps <- 200
# Number of cores
n_cores <- 55

# Create a new folder for this version
# Outputs, log files, Rscript, and calcualted estimators are under it.
mainfolder <- "~/Simulation_eQTL_replicates/scenario_files/Bootstrap_reestimation_effect_size/"
mydir <- paste0(mainfolder, "BootstrapQTL_v",versionbootsQ,"_",localcormet,"-",globalcormet,"_nomthres_",n_bootstraps,"_",typeeffectsize,"EffectSize/")
if(!file.exists(mydir)) {dir.create(mydir)}
setwd(mydir)
if(!file.exists("output_files")) {dir.create("output_files")}

# Log file, save version info
logfilename <- paste0(mydir,"LogFile_BootstrapQTL_rep",paste0(range(whichreps),collapse = "-"),".txt")
write(paste0("***** BootstrapQTL Version: ",versionbootsQ, " ***** "), file = logfilename, append = T)
write(paste0("      Date: ",Sys.time()), file = logfilename, append = T)
write(paste0("      Hierarchical correction method: ",localcormet,"-",globalcormet," \n"), file = logfilename, append = T)

# Function to write logs
writelog <- function(newinfo) {
  newline <- paste0("* Rep", IdxRep, ", Sample size=", SS, ", MAF=", MAF, "%: ", newinfo)
  write(newline, file = logfilename, append = T)
}

############################################


# Gene location
genepos <- read.table("~/Simulation_eQTL_replicates/hg19_gene_annotation/geneloc_for_ME.txt", header = T, stringsAsFactors = F)

# Six sample sizes
nothing <- foreach(IdxRep = rep(whichreps, each = 6), SS = rep(c(100,200,500,1000,2000,5000), howmanyreps)) %do% {
  # Load genotype data and SNP location
  setwd(paste0("~/Simulation_eQTL_replicates/simulation_results/rep", IdxRep, "/SS", SS, "/Genotype"))
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"
  snps$fileOmitCharacters = "NA"
  snps$fileSkipRows = 1
  snps$fileSkipColumns = 1
  snps$fileSliceSize = 2000
  snps$LoadFile("SNP.txt")
  snpspos <- read.table("snpsloc.txt", header = T, stringsAsFactors = F)

  # Read the output of eigenMT
  output_eigenMT <- read.table("eigenMT_output", header = T)
  output_eigenMT <- output_eigenMT[,c("gene","TESTS")]
  
  # Six MAFs
  nothing2 <- foreach(MAF = c(50,25,10,5,1,0.5) ) %do% {
    
    ##### Check whether bootstrap finished in this scenario #####
    setwd(paste0(mydir,"output_files"))
    if(file.exists(paste0("Rep",IdxRep,"SS",SS,"MAF",MAF,"_BootstrapQTL",n_bootstraps,"_",typeeffectsize,"EffectSize_output.txt"))) { return(NULL) }
    
    ##### Check - significant eGenes #####
    # If there's no significant eGene, next loop
    topSNPfilename <- paste0("~/Simulation_eQTL_replicates/scenario_files/Top_SNP_causal_SNP/Top_SNP_all_genes_causalSNP_info_",localcormet,"_",globalcormet,"_",typeeffectsize,"EffectSize/output_files/Rep",IdxRep,"SS",SS,"MAF",MAF,"_",typeeffectsize,"_Top_SNP_all_genes_causalSNP_info_",localcormet,globalcormet,".txt")
    topSNP <- read.table(topSNPfilename, header = T)
    # Signiciant eGenes
    if(sum(topSNP[,paste0("sig_",localcormet,globalcormet)]) == 0) {
      writelog("no significant eGenes.")
      return(NULL)
    }
    
    # Load expression data
    # Directory of the scenario
    setwd(paste0("~/Simulation_eQTL_replicates/simulation_results/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
    writelog("running bootstraps.")
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 2000
    gene$LoadFile("gene_expression.txt")
    
    ##### Run BootstrapQTL ####
    bootsout_shr <- BootstrapQTL(snps = snps, gene = gene, snpspos = snpspos, genepos = genepos, n_bootstraps = n_bootstraps, n_cores = n_cores, correction_type = "shrinkage", local_correction = "eigenMT", eigenMT_tests_per_gene = output_eigenMT)
    bootsout_oos <- BootstrapQTL(snps = snps, gene = gene, snpspos = snpspos, genepos = genepos, n_bootstraps = n_bootstraps, n_cores = n_cores, correction_type = "out_of_sample", local_correction = "eigenMT", eigenMT_tests_per_gene = output_eigenMT)
    bootsout_wei <- BootstrapQTL(snps = snps, gene = gene, snpspos = snpspos, genepos = genepos, n_bootstraps = n_bootstraps, n_cores = n_cores, correction_type = "weighted", local_correction = "eigenMT", eigenMT_tests_per_gene = output_eigenMT)
    # Merge three estimators into one table
    bootsout <- as.data.frame(bootsout_shr)
    bootsout_oos <- as.data.frame(bootsout_oos)
    bootsout_wei <- as.data.frame(bootsout_wei)
    bootsout$winners_curse <- NULL
    colnames(bootsout)[match(c("corrected_beta","correction_boots"),colnames(bootsout))] <- c("shrinkage","nsigboots_shr")
    
    bootsout <- merge(bootsout, bootsout_oos[,c("eGene","eSNPs","corrected_beta","correction_boots")], by = c("eGene","eSNPs"), all = T)
    colnames(bootsout)[(ncol(bootsout)-1):ncol(bootsout)] <- c("outof","nsigboots_oos")
    
    bootsout <- merge(bootsout, bootsout_wei[,c("eGene","eSNPs","corrected_beta","correction_boots")], by = c("eGene","eSNPs"), all = T)
    colnames(bootsout)[(ncol(bootsout)-1):ncol(bootsout)] <- c("weighted","nsigboots_wei")
    
    bootsout <- bootsout[order(bootsout$eGene_pval, bootsout$eSNP_pval),]
    if(!identical(bootsout$eGene, bootsout_shr$eGene)) { writelog("error: BootstrapQTL sig eAssos are not identical.") }
    
    # Data are under:
    setwd(mydir)
    setwd("output_files")
    write.table(bootsout, paste0("Rep",IdxRep,"SS",SS,"MAF",MAF,"_BootstrapQTL",n_bootstraps,"_",typeeffectsize,"EffectSize_output.txt"), quote = F, sep = "\t", row.names = F)
    
  }# End of 6 MAFs
}# End of 6 sizes


#----- 2017.10.17 -----
# R script is under: "~/Simulation_eQTL_replicates/scenario_files/Bootstrap_reestimation_effect_size/BootstrapQTL_v0.7.3_eigenMT-BH_nomthres_200_UnfixedEffectSize/BootstrapQTL_200_eigenMTBH_nomthre_three_estimators_UnfixedEffectSize_sim1-10.R"
# Outputs are under ./output_files/
# Logs are under the same folder.

# Running time 36.2hr using 55 cores, simulation 1-10, unfixed effect sizes.
#real	2173m8.829s
#user	30436m48.532s
#sys	81121m23.852s

# 2017-10-19 Bug:
# Missing some scenarios where there is no BonfBH sig eGenes but eigenMTBH sig eGenes.


