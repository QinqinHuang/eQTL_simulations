############################################
# Which Replicates:
whichreps <- 1:10
howmanyreps <- length(whichreps)
############################################

library(foreach)
library(doMC); registerDoMC(cores = 15)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")
options(stringsAsFactors = F)


#-----------------------------------------------
# Perform permutations using MatrixeQTL
#-----------------------------------------------
# Read gene location data
geneloc <- read.table("~/Simulation_eQTL_replicates/hg19_gene_annotation/geneloc_for_ME.txt", header = T)

# Six sample sizes and replicates
nothing <- foreach(SS = rep(c(5000,2000,1000,500,200,100), each = howmanyreps), IdxRep = rep(whichreps, 6)) %do% {
  # Genotype are under
  setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/Genotype/"))
  
  # Load genotype data
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"
  snps$fileOmitCharacters = "NA"
  snps$fileSkipRows = 1
  snps$fileSkipColumns = 1
  snps$fileSliceSize = 2000   # read file in pieces of 2000 rows
  snps$LoadFile("SNP.txt")
  
  # SNP location
  snpspos <- read.table("snpsloc.txt", header = T)
  
  # Six MAFs
  for(MAF in c(50,25,10,5,1,0.5)) {
    # Working directory
    setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
    
    # If the job has finished
    if(file.exists("Permutation_adjusted_pvalues_1k_MatrixeQTL.txt") & file.exists("Permuted_pvalues_1k.txt")) {next}
    
    # Load gene expression data
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 2000   # read file in pieces of 2000 rows
    gene$LoadFile("gene_expression.txt")
    
    # Run 1000 permutations
    permpval <- foreach(ii = 1:1000, .combine = rbind) %dopar% {
      # Shuffle sample lables for the expression data
      new.order <- sample(1:SS, SS, replace = F)
      gene.run <- gene$Clone()
      gene.run$ColumnSubsample(new.order)
      
      # Run the Matrix eQTL
      me <- Matrix_eQTL_main(
        snps = snps,
        gene = gene.run,
        pvOutputThreshold = 0,
        output_file_name.cis = NULL,
        pvOutputThreshold.cis = 1,
        snpspos = snpspos,
        genepos = geneloc,
        cisDist = 1000000,
        output_file_name = NULL,
        useModel = modelLINEAR,
        errorCovariance = numeric(),
        verbose = TRUE,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory = FALSE)
      
      # Keep the minium pvalues for each gene, and make it one row and 618 columns
      perm.curr <- as.data.frame(t(me$cis$min.pv.gene))
      return(perm.curr)
    } # End of permutations
    
    #-- Save permuted pvalues --
    write.table(permpval, "Permuted_pvalues_1k.txt", quote = F, sep = "\t", row.names = F)
    
    # Calculate permutation adjusted p-values for the best association of each gene
    # Read all permuted pvalues; header = F, otherwise '-' in gene names would be changed to '.'
    permuted_pvals <- read.table("Permuted_pvalues_1k.txt", header = F)
    permuted_pvals <- as.data.frame(t(permuted_pvals))
    # Check the number of permutations
    if(ncol(permuted_pvals) != 1001 | nrow(permuted_pvals) != 618) {
      write("Check the permuted pvalue file!", "Error.txt", append = T)
    }
    colnames(permuted_pvals)[-1] <- paste0("p", 1:1000)
    colnames(permuted_pvals)[1] <- "gene"
    
    # The observed minimal pvalue per gene
    me.ori <- readRDS("MatrixEQTL.output.rds")
    cis <- me.ori$cis$eqtls[, c("snps", "gene", "pvalue")]
    cis <- cis[which(!duplicated(cis$gene)),]
    cis <- merge(cis, permuted_pvals, by = "gene")
    
    # Calculate the permutation adjusted pvalues
    ppval <- cis[,1:3]
    cis <- as.matrix(cis[,-1:-3])
    ppval$per.pval <- sapply(1:nrow(cis), function(ii) {
      (sum(as.numeric(cis[ii,]) < as.numeric(ppval[ii,3])) + 1) / 1001
    })
    
    #-- Save permutation adjusted pvalues --
    write.table(ppval, "Permutation_adjusted_pvalues_1k_MatrixeQTL.txt", quote = F, sep = "\t", row.names = F)
    
  } # End of the loop for 6 MAFs
  
  return(NULL)
}
 
 


