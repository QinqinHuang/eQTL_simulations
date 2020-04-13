#--------------------------------------------------------------------------------
# Functions to calculate Bootstrap estimatiors using the same SNP
#--------------------------------------------------------------------------------
# Bootstrap (sameSNP) output is a list, 
# each element contains all significant eGenes and their top SNPs identified in original dataset. 
# Significant eGenes include false positives, we also re-estimate effect sizes for them.
# In each bootstrap, we should decide significant pairs.

library(foreach)

# weight
w <- 0.632

writelog <- function(newinfo) {
  newline <- paste0("* Rep", IdxRep, ", Sample size=", SS, ", MAF=", MAF, "%", ", effect size=", effectsize, ": ", newinfo)
  #write(newline, file = "~/Simulation_eQTL_replicates/results_data_100reps/LogFIles/Bootstrap500_sameSNP_nomthre_BonfBH.txt", append = T)
  cat(newline, "\n")
}


#---------- 1. shrinkage estimator ----------
# "eGenes": Naive estimator for TP eGenes, contains at least "gene", "Naive", and "nom_threshold"
# "boots": a list, each element contains a table ("sameSNP")
# "force_sign": if it is true, we don't just use absolute values, sign will be consistent.
#               e.g. if Naive estimator is negative, we assume betas from detection and estimation group are both negative.
# Return the whole table

shrinkage_sameSNP <- function(eGenes, boots, force_sign = FALSE) {
  # Check gene order
  if(!identical(eGenes$gene, boots[[1]]$gene)) {
    writelog("Error (calculation) - significant eGenes are not the same with that in bootstrap!")
    return(eGenes)
  }
  # Calcualte shrinkage factor/winner's curse effect for each bootstrap, NA if not signficant
  shrinkage.factor <- foreach(ii = 1:length(boots), .combine = cbind) %do% {
    tempdd <- boots[[ii]]
    if(force_sign) {
      tempdd$beta.diff <- abs(tempdd$betaD)*sign(eGenes$Naive) - abs(tempdd$betaE)*sign(eGenes$Naive)
    } else {
      tempdd$beta.diff <- abs(tempdd$betaD) - abs(tempdd$betaE)    # Absolute
    }
    # Keep those are still significant
    tempdd$beta.diff[which(tempdd$pvalue > eGenes$nom_threshold)] <- NA
    return(tempdd$beta.diff)
  }
  # Each row is an eGene and columns are for bootstraps
  # Number of significant bootstraps
  if(!force_sign) { eGenes$nsig_boots <- apply(shrinkage.factor, 1, function(d) {sum(complete.cases(d))} ) }
  # Average of shrinkage factor/winner's curse effect
  if(!force_sign) {
    eGenes$shrinkage_factor <- apply(shrinkage.factor, 1, function(d) {mean(d[complete.cases(d)])})
    eGenes$shrinkage_same <- abs(eGenes$Naive) - eGenes$shrinkage_factor    # Absolute
  } else {
    eGenes$shrinkage_factor_forcesign <- apply(shrinkage.factor, 1, function(d) {mean(d[complete.cases(d)])})
    eGenes$shrinkage_same_forcesign <- eGenes$Naive - eGenes$shrinkage_factor_forcesign
  }
  return(eGenes)
}


#---------- 2. out-of-sample estimator ----------
# "eGenes": Naive estimator, 200 true eGenes, contains at least "gene", "Naive", and "nom_threshold"
# "boots": a list, each element contains a table ("sameSNP")
# "force_sign": if it is true, we don't just use absolute values, sign will be consistent.
#               e.g. if Naive estimator is negative, we assume betas from detection and estimation group are both negative.
# Return the whole table

outof_sameSNP <- function(eGenes, boots, force_sign = FALSE) {
  # Check gene order
  if(!identical(eGenes$gene, boots[[1]]$gene)) {
    writelog("Error (calculation) - significant eGenes are not the same with that in bootstrap!")
    return(eGenes)
  }
  # Calcualte shrinkage factor/winner's curse effect for each bootstrap, NA if not signficant
  betas_estimation <- foreach(ii = 1:length(boots), .combine = cbind) %do% {
    tempdd <- boots[[ii]]
    # Keep those are still significant
    tempdd$betaE[which(tempdd$pvalue > eGenes$nom_threshold)] <- NA
    if(force_sign) {
      return(abs(tempdd$betaE) * sign(eGenes$Naive))
    } else {
      return(abs(tempdd$betaE))   # Absolute
    }
  }
  # Each row is an eGene and columns are for bootstraps
  # Average of out-of-sample estimates
  if(!force_sign) {
    eGenes$outof_same <- apply(betas_estimation, 1, function(d) {mean(d[complete.cases(d)])})
  } else {
    eGenes$outof_same_forcesign <- apply(betas_estimation, 1, function(d) {mean(d[complete.cases(d)])})
  }
  return(eGenes)
}



#---- 2017.09.13 -----
# These two functions are part of ~/Simulation_eQTL_replicates/Scripts/002_MyFunctions_Bootstrap_reestimation.R




