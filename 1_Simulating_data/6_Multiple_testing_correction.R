############################################
# Which Reps:
whichreps <- 1:100
howmanyreps <- length(whichreps)
############################################

library(foreach)
library(doMC); registerDoMC(cores = 35)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")
options(stringsAsFactors = F)


#----------------------------------------------------------------------------------------------------------------
# Get the list of significant SNP-gene associations using different methods to correct for multiple testing
#----------------------------------------------------------------------------------------------------------------
# (1) pooled method
# (2) hierarchical testing procedure

# Handle 36 scenarios for each replicate 
nothing <- foreach(IdxRep = rep(whichreps, each = 36), SS = rep(rep(c(5000,2000,1000,500,200,100), each = 6), howmanyreps), MAF = rep(c(50,25,10,5,1,0.5), 6*howmanyreps) ) %dopar% {
  
  # Directory of current scenario
  setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
  
  # Check whether MTC has finished for this scenario
  if(file.exists("cis_eAssociations_modelLinear_CP05.txt")) {return(NULL)}
  
  # Keep a log
  write(paste("Working on: Rep", IdxRep, "SS", SS, "MAF", MAF), file = paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/LogFIles/MTC_whichscenario.txt"), append = T)
  
  # Read the output of eigenMT
  eMT200_output <- read.table("../Genotype/eigenMT_output", header = T)[,c("gene","TESTS")]
  
  # Read output of FastQTL
  per1k <- read.table("Permutation_1k.txt.gz", header = F)
  
  # Read true causal eSNPs and eGenes
  eQTLsnp <- read.table("True_eQTLs.txt", header = T)
  # Read LD structure around causal eSNPs
  LDdata001 <- read.table("LD_eSNP_r2_001_window_2Mb.ld.gz", header = T)
  # eGene of the SNP_A (causal eSNP)
  LDdata001 <- merge(LDdata001, eQTLsnp[,c("SNP","eGene")], by.x = "SNP_A", by.y = "SNP", all.x = T)
  
  # Read the output of MatrixEQTL
  me <- readRDS("MatrixEQTL.output.rds")
  # P-values for all tests
  cis <- me$cis$eqtls[, c("snps","gene","pvalue")]
  # Number of genes tested
  ntestedgene <- length(unique(cis$gene))
  
  #--------------------------------------------------------------------
  #----- Part 1: Apply a method to the entire set of the p-values -----
  ##### Pooled ST FDR procedure #####
  STfdr <- cis
  STfdr$pooledST <- qvalue(STfdr$pvalue)$qvalues
  STfdr <- STfdr[which(STfdr$pooledST <= 0.05),]
  eAsso <- STfdr
  
  ##### Pooled BH FDR procedure #####
  fdr <- me$cis$eqtls[,c("snps","gene","FDR")]
  fdr <- fdr[which(fdr$FDR <= 0.05),]
  eAsso <- merge(eAsso, fdr, by = c("snps","gene"), all = T)
  names(eAsso)[ncol(eAsso)] <- "pooledBH"
  
  ##### Pooled BY FDR procedure #####
  BYfdr <- cis
  BYfdr$pooledBY <- p.adjust(BYfdr$pvalue, method = "BY")
  BYfdr <- BYfdr[which(BYfdr$pooledBY <= 0.05), -3]
  eAsso <- merge(eAsso, BYfdr, by = c("snps","gene"), all = T)
  
  ##### Pooled Bonferroni correction #####
  bonf <- cis
  bonf$pooledBonf <- p.adjust(bonf$pvalue, method = "bonferroni")
  bonf <- bonf[which(bonf$pooledBonf <= 0.05), -3]
  eAsso <- merge(eAsso, bonf, by = c("snps","gene"), all = T)
  
  #--------------------------------------------------------------------
  #----- Part 2: Hierarchical testing procedure -----
  ##### FDR controlling procedures used in both steps #####
  # ST-ST
  STST <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "ST", step2_met = "ST")
  eAsso <- merge(eAsso, STST, by = c("snps","gene"), all = T)
  # ST-BH
  STBH <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "ST", step2_met = "BH")
  eAsso <- merge(eAsso, STBH, by = c("snps","gene"), all = T)
  # ST-BY
  STBY <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "ST", step2_met = "BY")
  eAsso <- merge(eAsso, STBY, by = c("snps","gene"), all = T)
  # ST-Bonferroni
  STBonf <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "ST", step2_met = "Bonferroni")
  eAsso <- merge(eAsso, STBonf, by = c("snps","gene"), all = T)
  
  # BH-ST
  BHST <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "BH", step2_met = "ST")
  eAsso <- merge(eAsso, BHST, by = c("snps","gene"), all = T)
  # BH-BH
  BHBH <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "BH", step2_met = "BH")
  eAsso <- merge(eAsso, BHBH, by = c("snps","gene"), all = T)
  # BH-BY
  BHBY <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "BH", step2_met = "BY")
  eAsso <- merge(eAsso, BHBY, by = c("snps","gene"), all = T)
  # BH-Bonferroni
  BHBonf <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "BH", step2_met = "Bonferroni")
  eAsso <- merge(eAsso, BHBonf, by = c("snps","gene"), all = T)
  
  # BY-ST
  BYST <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "BY", step2_met = "ST")
  eAsso <- merge(eAsso, BYST, by = c("snps","gene"), all = T)
  # BY-BH
  BYBH <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "BY", step2_met = "BH")
  eAsso <- merge(eAsso, BYBH, by = c("snps","gene"), all = T)
  # BY-BY
  BYBY <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "BY", step2_met = "BY")
  eAsso <- merge(eAsso, BYBY, by = c("snps","gene"), all = T)
  # BY-Bonferroni
  BYBonf <- mul_corr_FDR_in_hierarchical(alltests = cis, step1_met = "BY", step2_met = "Bonferroni")
  eAsso <- merge(eAsso, BYBonf, by = c("snps","gene"), all = T)
  
  ##### Bonferroni in Step 1; FDR/Bonferroni in Step2 #####
  # Bonferroni-ST
  BonfST <- mul_corr_Bonf_in_hierarchical(alltests = cis, step2_met = "ST")
  eAsso <- merge(eAsso, BonfST, by = c("snps","gene"), all = T)
  # Bonferroni-BH
  BonfBH <- mul_corr_Bonf_in_hierarchical(alltests = cis, step2_met = "BH")
  eAsso <- merge(eAsso, BonfBH, by = c("snps","gene"), all = T)
  # Bonferroni-BY
  BonfBY <- mul_corr_Bonf_in_hierarchical(alltests = cis, step2_met = "BY")
  eAsso <- merge(eAsso, BonfBY, by = c("snps","gene"), all = T)
  # Bonferroni-Bonferroni
  BonfBonf <- mul_corr_Bonf_in_hierarchical(alltests = cis, step2_met = "Bonferroni")
  eAsso <- merge(eAsso, BonfBonf, by = c("snps","gene"), all = T)
  
  ##### eigenMT in Step 1; FDR in Step2 #####
  # eigenMT-ST
  eMTST_returnlist <- mul_corr_eigenMT(alltests = cis, eMT_output = eMT200_output, step2_met = "ST")
  eAsso <- merge(eAsso, eMTST_returnlist[[1]], by = c("snps","gene"), all = T)
  NominalPvalueThreshold <- eMTST_returnlist[[2]]
  # eigenMT-BH
  eMTBH_returnlist <- mul_corr_eigenMT(alltests = cis, eMT_output = eMT200_output, step2_met = "BH")
  eAsso <- merge(eAsso, eMTBH_returnlist[[1]], by = c("snps","gene"), all = T)
  NominalPvalueThreshold <- merge(NominalPvalueThreshold, eMTBH_returnlist[[2]], by = "gene", all = T)
  # eigenMT-BY
  eMTBY_returnlist <- mul_corr_eigenMT(alltests = cis, eMT_output = eMT200_output, step2_met = "BY")
  eAsso <- merge(eAsso, eMTBY_returnlist[[1]], by = c("snps","gene"), all = T)
  NominalPvalueThreshold <- merge(NominalPvalueThreshold, eMTBY_returnlist[[2]], by = "gene", all = T)
  # eigenMT-Bonferroni
  eMTBonf_returnlist <- mul_corr_eigenMT(alltests = cis, eMT_output = eMT200_output, step2_met = "Bonferroni")
  eAsso <- merge(eAsso, eMTBonf_returnlist[[1]], by = c("snps","gene"), all = T)
  NominalPvalueThreshold <- merge(NominalPvalueThreshold, eMTBonf_returnlist[[2]], by = "gene", all = T)
  
  ##### Beta approximation permutations (1k) in Step 1; FDR in Step2 #####
  # BPerm1k-ST
  bper1kST_return_list <- mul_corr_BPerm(alltests = cis, fastqtl_output = per1k, scheme = "BPerm1k", step2_met = "ST", samplesize = SS)
  eAsso <- merge(eAsso, bper1kST_return_list[[1]], by = c("snps","gene"), all = T)
  NominalPvalueThreshold <- merge(NominalPvalueThreshold, bper1kST_return_list[[2]], by = "gene", all = T)
  # BPerm1k-BH
  bper1kBH_return_list <- mul_corr_BPerm(alltests = cis, fastqtl_output = per1k, scheme = "BPerm1k", step2_met = "BH", samplesize = SS)
  eAsso <- merge(eAsso, bper1kBH_return_list[[1]], by = c("snps","gene"), all = T)
  NominalPvalueThreshold <- merge(NominalPvalueThreshold, bper1kBH_return_list[[2]], by = "gene", all = T)
  # BPerm1k-BY
  bper1kBY_return_list <- mul_corr_BPerm(alltests = cis, fastqtl_output = per1k, scheme = "BPerm1k", step2_met = "BY", samplesize = SS)
  eAsso <- merge(eAsso, bper1kBY_return_list[[1]], by = c("snps","gene"), all = T)
  NominalPvalueThreshold <- merge(NominalPvalueThreshold, bper1kBY_return_list[[2]], by = "gene", all = T)
  # BPerm1k-Bonferroni
  bper1kBonf_return_list <- mul_corr_BPerm(alltests = cis, fastqtl_output = per1k, scheme = "BPerm1k", step2_met = "Bonferroni", samplesize = SS)
  eAsso <- merge(eAsso, bper1kBonf_return_list[[1]], by = c("snps","gene"), all = T)
  NominalPvalueThreshold <- merge(NominalPvalueThreshold, bper1kBonf_return_list[[2]], by = "gene", all = T)
  
  
  #----- We also apply method: adaptive permutations and manual permutations in Rep 1-10 -----
  if(IdxRep <= 10) {
    
    ##### Beta approximation permutations (100-10k, adaptive scheme) in Step 1; FDR in Step2 #####
    # Output of FastQTL
    aper10k <- read.table("Permutation_100-10k.txt.gz", header = F)
    
    # APerm10k-ST
    aper10kST_return_list <- mul_corr_BPerm(alltests = cis, fastqtl_output = aper10k, scheme = "APerm10k", step2_met = "ST", samplesize = SS)
    eAsso <- merge(eAsso, aper10kST_return_list[[1]], by = c("snps","gene"), all = T)
    NominalPvalueThreshold <- merge(NominalPvalueThreshold, aper10kST_return_list[[2]], by = "gene", all = T)
    # APerm10k-BH
    aper10kBH_return_list <- mul_corr_BPerm(alltests = cis, fastqtl_output = aper10k, scheme = "APerm10k", step2_met = "BH", samplesize = SS)
    eAsso <- merge(eAsso, aper10kBH_return_list[[1]], by = c("snps","gene"), all = T)
    NominalPvalueThreshold <- merge(NominalPvalueThreshold, aper10kBH_return_list[[2]], by = "gene", all = T)
    # APerm10k-BY
    aper10kBY_return_list <- mul_corr_BPerm(alltests = cis, fastqtl_output = aper10k, scheme = "APerm10k", step2_met = "BY", samplesize = SS)
    eAsso <- merge(eAsso, aper10kBY_return_list[[1]], by = c("snps","gene"), all = T)
    NominalPvalueThreshold <- merge(NominalPvalueThreshold, aper10kBY_return_list[[2]], by = "gene", all = T)
    # APerm10k-Bonf
    aper10kBonf_return_list <- mul_corr_BPerm(alltests = cis, fastqtl_output = aper10k, scheme = "APerm10k", step2_met = "Bonferroni", samplesize = SS)
    eAsso <- merge(eAsso, aper10kBonf_return_list[[1]], by = c("snps","gene"), all = T)
    NominalPvalueThreshold <- merge(NominalPvalueThreshold, aper10kBonf_return_list[[2]], by = "gene", all = T)
    
    ##### Permutations (1k) using MatrixEQTL in Step 1; FDR in Step2 #####
    # Read permutated pvalues
    permuted_pvals <- read.table("Permuted_pvalues_1k.txt", header = F)
    # Read the permutation adjusted pvalues for each gene
    ppval <- read.table("Permutation_adjusted_pvalues_1k_MatrixeQTL.txt", header = T)
    
    # Perm1k-ST
    perm1kST_return_list <- mul_corr_ManualPerm(alltests = cis, allpermutedpval = permuted_pvals, ppval_gene = ppval, step2_met = "ST")
    eAsso <- merge(eAsso, perm1kST_return_list[[1]], by = c("snps","gene"), all = T)
    NominalPvalueThreshold <- merge(NominalPvalueThreshold, perm1kST_return_list[[2]], by = "gene", all = T)
    # Perm1k-BH
    perm1kBH_return_list <- mul_corr_ManualPerm(alltests = cis, allpermutedpval = permuted_pvals, ppval_gene = ppval, step2_met = "BH")
    eAsso <- merge(eAsso, perm1kBH_return_list[[1]], by = c("snps","gene"), all = T)
    NominalPvalueThreshold <- merge(NominalPvalueThreshold, perm1kBH_return_list[[2]], by = "gene", all = T)
    # Perm1k-BY
    perm1kBY_return_list <- mul_corr_ManualPerm(alltests = cis, allpermutedpval = permuted_pvals, ppval_gene = ppval, step2_met = "BY")
    eAsso <- merge(eAsso, perm1kBY_return_list[[1]], by = c("snps","gene"), all = T)
    NominalPvalueThreshold <- merge(NominalPvalueThreshold, perm1kBY_return_list[[2]], by = "gene", all = T)
    # No Perm1k-Bonferroni because the possible lowest permutation adjusted p-value is 1/1001 and nothing will be significant after Bonferroni penalty.
    
  }
  
  
  ##### Add LD r2 between each significant eSNP and the causal eSNP of the eGene #####
  # Merge LD measurements into eAsso
  eAsso <- merge(eAsso, LDdata001[,c("SNP_B","eGene","R2")], by.x = c("gene","snps"), by.y = c("eGene","SNP_B"), all.x = T)
  # R2 is NA could be due to (1) r2 <0.01 or (2) the gene is not among the 200 true eGenes, so no causal eSNP for it.
  names(eAsso)[ncol(eAsso)] <- "R2_causaleSNP"
  
  
  ##### Save results ##### 
  # Save nominal pvalue thresholds
  write.table(NominalPvalueThreshold, "Nominal_thresholds_gene_permutation_eigenMT.txt", quote = F, sep = "\t", row.names = F)
  
  # Save the list of significant cis eQTLs
  write.table(eAsso, "cis_eAssociations_modelLinear_CP05.txt", quote = F, sep = "\t", row.names = F)
  
  return(NULL)
  
}





