############################################
# Which Replicates:
whichreps <- 11:100
howmanyreps <- length(whichreps)
############################################

library(foreach)
library(doMC); registerDoMC(cores = 15)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")
options(stringsAsFactors = F)


#----------------------------------------------------------------------
# Run MatrixeQTL
#----------------------------------------------------------------------
# Read gene location file
geneloc <- read.table("~/Simulation_eQTL_replicates/hg19_gene_annotation/geneloc_for_ME.txt", header = T)
dd <- data.frame()

# Six sample sizes and replicates
nothing <- foreach(SS = rep(c(5000,2000,1000,500,200,100), each = howmanyreps), IdxRep = rep(whichreps, 6)) %dopar% {
  ## Check whether we have results from MatrixEQTL ##
  setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS))
  if( all(file.exists(paste0("SS",SS,"maf50/MatrixEQTL.output.rds")),
         file.exists(paste0("SS",SS,"maf25/MatrixEQTL.output.rds")),
         file.exists(paste0("SS",SS,"maf10/MatrixEQTL.output.rds")),
         file.exists(paste0("SS",SS,"maf5/MatrixEQTL.output.rds")),
         file.exists(paste0("SS",SS,"maf1/MatrixEQTL.output.rds")),
         file.exists(paste0("SS",SS,"maf0.5/MatrixEQTL.output.rds")) ) ) {
    return(NULL)
  }

  cat("Rep", IdxRep, "SS", SS, "\n")
  
  # Load genotype data for this sample size
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"
  snps$fileOmitCharacters = "NA"
  snps$fileSkipRows = 1
  snps$fileSkipColumns = 1
  snps$fileSliceSize = 2000   # read file in pieces of 2000 rows
  snps$LoadFile("Genotype/SNP.txt")
  
  snpspos <- read.table("Genotype/snpsloc.txt", header = T)
  
  # Six MAFs
  for(MAF in c(50,25,10,5,1,0.5)) {
    # Working directory
    setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
    
    # Whether MatrixeQTL has finished for this scenario.
    if(file.exists("MatrixEQTL.output.rds")) { next }
    
    # Run MatrixeQTL
    if(!file.exists("MatrixEQTL.output.rds")) {
      # Run MatrixEQTL
      me <- RunMatrixeQTL_withloaded_snps(geno = snps, snpspos = snpspos, expr_path = "gene_expression.txt", genepos = geneloc)
      saveRDS(me, file = "MatrixEQTL.output.rds")
      
      # Save results for eigenMT as input for each sample size.
      if(MAF == 50) {
        cis <- me$cis$eqtls[, c("snps","gene","pvalue")]
        write.table(cis, file = "MatrixEQTL_all_cis_tests_pval.txt", quote = F, sep = "\t", row.names = F)
      }
    }
    
  } # End of 6 MAFs
  
  return(NULL)
}


