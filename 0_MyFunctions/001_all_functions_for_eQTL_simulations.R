#----------------------------------------------------------------------------------------------
# Functions in eQTL simulation studies
#----------------------------------------------------------------------------------------------
library(MatrixEQTL)
library(qvalue)
library(foreach)

options(stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------------------------
# Descriptions of functions
#
# ---------- 1. Simulate Data ----------
# //whichfolder(II)
# rep1-10 and rep11-100 are under different folders.
#
# //simulate_exp_data(geno_data, gene_pos, true_eQTLs_eGenes, true_effectsize = "frombetadist", needfastqtl = FALSE)
# Aim: to generate expression data with some of them regulated by eQTLs.
# Retrun expression matrix
# 
# //RunMatrixeQTL(geno_path, expr_path, snppos_path, genepos, whichmodel = modelLINEAR)
# Aim: Load data and run MatrixEQTL.
# Return MatrixEQTL results
# 
# //RunMatrixeQTL_withloaded_snps(geno, snpspos, expr_path, genepos, whichmodel = modelLINEAR)
# Aim: Run MatrixEQTL when genotype data are loaded.
# Return MatrixEQTL results
#
# ---------- 2. Multiplicity Correction ----------
# //mul_corr_eigenMT(alltests, eMT_output, fdr_desired = 0.05, step2_met = "BH")
# Method: eigenMT-ST, eigenMT-BH, eigenMT-BY, or eigenMT-Bonferroni
# Return eAssociations and nominal thresholds
#
# //mul_corr_BPerm(alltests, fastqtl_output, scheme = "BPerm1k", fdr_desired = 0.05, step2_met = "BH", samplesize)
# Method: Beta approximation based permutations (BPerm1k-BH, BPerm1k-ST, BPerm1k-BY, or BPerm1k-Bonferroni; APerm10k-BH, APerm10k-ST, APerm10k-BY, or APerm10k-Bonferroni)
# Return eAssociations and nominal thresholds
#
# //mul_corr_ManualPerm(alltests, allpermutedpval, ppval_gene, fdr_desired = 0.05, step2_met = "BH")
# Method: Permutations using direct scheme (Perm1k-BH, Perm1k-ST, or Perm1k-BY)
# Return eAssociations and nominal thresholds
#
# //mul_corr_FDR_in_hierarchical(alltests, step1_met, fdr_desired = 0.05, step2_met = "BH", mulcores = F)
# Method: ST-ST, ST-BH, ST-BY, ST-Bonferroni, BH-ST, BH-BH, BH-BY, BH-Bonferroni, BY-ST, BY-BH, BY-BY, BY-Bonferroni (can be done in parallel)
# Return eAssociations
#
# //mul_corr_Bonf_in_hierarchical(alltests, fdr_desired = 0.05, step2_met = "BH")
# Method: Bonferroni-ST, Bonferroni-BH, Bonferroni-BY, or Bonferroni-Bonferroni
# Return eAssociations
# //mul_corr_Bonf_in_hierarchical_nomthr(alltests, fdr_desired = 0.05, step2_met = "BH")
# Method: Bonferroni-ST, Bonferroni-BH, Bonferroni-BY, or Bonferroni-Bonferroni
# Return eAssociations and nominal thresholds
#
# ----- Functions used in previous mul_corr functions -----
# //getThreshold(d, opt_fdr = 0.05)
# Return the empirical pvalue threshold based on adjusted pvalue after correction for multiple genes.
#
# //getAssociation(d1, d2, met2)
# After calculating step1-adjusted pvalues, get the significant associations given a method used in step2 (ST/BH/BY/Bonf).
# Return significant associations and nominal pvalues.
#
# //getAssociation.withThresholds <- function(d1, d2)
# Get significant associations based on nominal pvalue thresholds.
# Return significant associations and nominal pvalues.
#
# //getNominalThresholds.bperm(perm, opt_fdr = 0.05, met2, SS)
# Calculate nominal pvalue thresholds based on beta approximation of permutation null distribution
#
# ---------- 3. Calculate Sensitivity and FDR ----------
# ## calculate FDR and TPR for each replicate separately and then take the average ##
# //calculate_TPR(resultsdata, methodlist = "all")
# Return TPR, row order is the same with the input
#
# //calculate_FDR(resultsdata, methodlist = "all")
# Return FDR, row order is the same with the input
#
# //calculate_average_values(data_reps, n_reps)
# Aim: to calculate the average values (FDR/TPR) for each scenario
#
# ## calculate FDR/TPR based on all discoveries from all replicates ##
# //calculateTPRFDR_allreps(resultsdata, methodlist = "all", FDRorTPR)
# Return TPR or FDR
#
# ---------- 4. Calculate average/sum for each scenario ----------
# //calculate_ave_ignoreNAs(data_reps, scenario_parameters = c("SS","MAF"))
# Calculate the average of each column for each scenario separately. 
#
# //calculate_sum_ignoreNAs(data_reps, scenario_parameters = c("SS","MAF"))
# Calculate the sum of each column for each scenario separately.
#
# ---------------------------------------------------------------------------------------------


# -------------------- 1. Simulate Data --------------------
whichfolder <- function(II) {
  if(II <= 10) {return("simulation_results")} else if(II >=11) {return("simulation_90_more")}
}


# Generate expression data based on genotypes and simulated eQTLs
# "geno_data": first column must be SNP id named "SNP", the rest columns are samples; each line is the genotypes of one SNP.
# "gene_pos": contain a column named "geneid"
# "true_eQTLs_eGenes": contain columns named "eGene", "SNP", "POS", "sim_beta"
# "true_effectsize": if it's not provided, use the values from true_eQTLs_eGenes$sim_beta
# "needfastqtl": T/F, if it's true, generate bed.gz file under the current folder. Default is FALSE.

simulate_exp_data <- function(geno_data, gene_pos, true_eQTLs_eGenes, true_effectsize = "frombetadist", needfastqtl = FALSE) {
  # Check the first column is SNP ID
  if(names(geno_data)[1] != "SNP") {
    return(NA)
  }
  # Number of samples
  num_indiv <- ncol(geno_data) - 1 
  
  # Number of genes
  num_genes <- nrow(gene_pos)
  
  # Number of true eGenes
  num_eGenes <- nrow(true_eQTLs_eGenes)
  
  # The error term follows the standard normal distribution
  sim_expr <- as.data.frame(matrix(rnorm(num_genes * num_indiv, mean = 0, sd = 1), nrow = num_genes, ncol = num_indiv))
  
  # The first few genes are eGenes and the rest genes are null genes
  rownames(sim_expr)[1:num_eGenes] <- true_eQTLs_eGenes$eGene
  rownames(sim_expr)[(num_eGenes + 1):num_genes] <- setdiff(gene_pos$geneid, true_eQTLs_eGenes$eGene)
  # Column names are individual IDs
  colnames(sim_expr) <- colnames(geno_data)[2:(num_indiv+1)]
  
  # The genetically regulated component of the expression data
  regulated_expr <- merge(true_eQTLs_eGenes[,c("SNP","POS")], geno_data, by = "SNP")
  # Sort according to the SNP positions (the same order with that in true_eQTLs_eGenes)
  regulated_expr <- regulated_expr[order(regulated_expr$POS),]
  true_eQTLs_eGenes <- true_eQTLs_eGenes[order(true_eQTLs_eGenes$POS),]
  # If the true effect size follows beta distribution; or a fixed value
  if(true_effectsize == "frombetadist") {
    regulated_expr <- regulated_expr[,3:ncol(regulated_expr)] * true_eQTLs_eGenes$sim_beta
  } else if(length(true_effectsize) == 1 & true_effectsize > 0) {
    regulated_expr <- regulated_expr[,3:ncol(regulated_expr)] * true_effectsize
  }
  sim_expr[1:num_eGenes,] <- regulated_expr + sim_expr[1:num_eGenes,]
  
  # If we need to run fastQTL later. Gene stop is 5 bp distance from gene start.
  if(needfastqtl) {
    # Generate the input phenotype bed file of FastQTL
    sim_expr$geneid <- rownames(sim_expr)
    expr.fastQTL <- merge(gene_pos, sim_expr, by = "geneid")
    expr.fastQTL <- expr.fastQTL[,c(2:4,1,5:ncol(expr.fastQTL))]
    # (FastQTL defines cisdistance as genotype_pos[g] - phenotype_start[p], and compares it with ciswindow)
    expr.fastQTL$right <- expr.fastQTL$left + 5
    colnames(expr.fastQTL)[1:4] <- c("#Chr", "start", "end", "ID")
    expr.fastQTL <- expr.fastQTL[order(expr.fastQTL$start),]
    write.table(expr.fastQTL, file = "gene_expression.bed", sep = "\t", row.names = F, quote = F)
    # FastQTL input file
    system("bgzip gene_expression.bed && tabix -p bed gene_expression.bed.gz")
    sim_expr$geneid <- NULL
  }
  
  return(sim_expr)
}


# Load data and run MatrixeQTL
# "geno_path", "expr_path", "snppos_path": the path of the data
# "genepos": gene location data
# "whichmodel": default is linear model, and also use ANOVA

RunMatrixeQTL <- function(geno_path, expr_path, snppos_path, genepos, whichmodel = modelLINEAR) {
  # Load the data
  snps = SlicedData$new()
  snps$fileDelimiter = "\t"
  snps$fileOmitCharacters = "NA"    # missing values
  snps$fileSkipRows = 1
  snps$fileSkipColumns = 1
  snps$fileSliceSize = 2000   # read file in pieces of 2000 rows
  snps$LoadFile(geno_path)
  
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"
  gene$fileOmitCharacters = "NA"
  gene$fileSkipRows = 1
  gene$fileSkipColumns = 1
  gene$fileSliceSize = 2000   # read file in pieces of 2000 rows
  gene$LoadFile(expr_path)
  
  snpspos <- read.table(snppos_path, header = T, stringsAsFactors = F)

  # linear model
  me <- Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    pvOutputThreshold = 0,
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = 1,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = 1000000,
    output_file_name = NULL,
    useModel = whichmodel,
    errorCovariance = numeric(),
    verbose = TRUE,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE)
  
  return(me)
}


# Load only expression data and run MatrixeQTL
# "geno", "snpspos": loaded genotype data and snp location data
# "expr_path": the path of the data
# "genepos": gene location data
# "whichmodel": default is linear model, and also use ANOVA

RunMatrixeQTL_withloaded_snps <- function(geno, snpspos, expr_path, genepos, whichmodel = modelLINEAR) {
  # Load expression data
  gene = SlicedData$new()
  gene$fileDelimiter = "\t"
  gene$fileOmitCharacters = "NA"
  gene$fileSkipRows = 1
  gene$fileSkipColumns = 1
  gene$fileSliceSize = 2000   # read file in pieces of 2000 rows
  gene$LoadFile(expr_path)
  
  # linear model
  me <- Matrix_eQTL_main(
    snps = geno,
    gene = gene,
    pvOutputThreshold = 0,
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = 1,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = 1000000,
    output_file_name = NULL,
    useModel = whichmodel,
    errorCovariance = numeric(),
    verbose = TRUE,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE)
  
  return(me)
}


# -------------------- 2. Multiplicity Correction --------------------
# Multiplicity correction using eigenMT in Step 1; BH or ST in Step 2
# "alltests": output of MatrixEQTL containing as least "[snps]", "[gene]", and "[pvalue]"
# "eMT_output": contains two column named "[gene]" and "[TESTS]"
# "fdr_desired": the desired level of gene-level FDR, default value 5%
# "step2_met": methods to control for multiple genes - "BH" (default), "ST", "BY", or "Bonferroni"

mul_corr_eigenMT <- function(alltests, eMT_output, fdr_desired = 0.05, step2_met = "BH") {
  cat("  ** Step 1: eigenMT with window size 200, Step 2: ", step2_met, "-FDR ** \n", sep = "")
  # Keep the best association for each gene
  eMT_output <- merge(eMT_output, alltests[which(!duplicated(alltests$gene)), c("gene","pvalue")], by = "gene")
  # Calculate eMT adjusted pvalues for each top association
  eMT_output$eMT <- eMT_output$pvalue * eMT_output$TESTS
  eMT_output$eMT[which(eMT_output$eMT > 1)] <- 1
  
  # Correct for multiple genes
  if(step2_met == "BH") {
    eMT_output$gene.fdr <- p.adjust(eMT_output$eMT, "fdr")
  } else if(step2_met == "ST") {
    eMT_output$gene.fdr <- qvalue(eMT_output$eMT)$qvalues
    cat("    calculating qvalues in Step 2, pi0 =", qvalue(eMT_output$eMT)$pi0, "\n")
  } else if(step2_met == "BY") {
    eMT_output$gene.fdr <- p.adjust(eMT_output$eMT, "BY")
  } else if(step2_met == "Bonferroni") {
    eMT_output$gene.fdr <- p.adjust(eMT_output$eMT, "bonferroni")
  }
  
  # Get the empirical p-value threshold
  threshold <- getThreshold(eMT_output[,c("eMT","gene.fdr")], opt_fdr = fdr_desired)
  cat("    The empirical p-value threshold =", threshold, "\n")
  # Calculate the nominal pvalue threshold for each gene
  nthre.eMT <- data.frame(gene = eMT_output$gene, eMT = threshold / eMT_output$TESTS)
  names(nthre.eMT)[2] <- paste0("eigenMT.", step2_met)
  # Get the significant pairs based on the nominal thresholds
  dd <- getAssociation.withThresholds(d1 = alltests[, c("snps","gene","pvalue")], d2 = nthre.eMT)
  names(dd)[3] <- paste0("eigenMT.", step2_met)
  
  return(list(dd,nthre.eMT))
}

#####  # Nonimal p-value thresholds
#####  if(identical(eAsso_eMT_BH_list[[2]], eAsso_eMT_ST_list[[2]])) {
#####    nom_thre <- eAsso_eMT_BH_list[[2]]
#####    names(nom_thre)[2] <- "eMT_BHorST"
#####  } else {
#####    nom_thre <- merge(eAsso_eMT_BH_list[[2]], eAsso_eMT_ST_list[[2]], by = "gene")
#####    names(nom_thre)[-1] <- c("eMT_BH", "eMT_ST")
#####  }


# Multiplicity correction using permutations based on beta approximation, BPerm1k or APerm10k
# "alltests": output of MatrixEQTL containing as least "[snps]", "[gene]", and "[pvalue]"
# "fastqtl_output": 1k permutations or 100-10k permutations
# "scheme": fixed number of permutations ("BPerm1k", default), or 100-10k permutations ("APerm10k"). just for names
# "fdr_desired": the desired level of gene-level FDR, default value 5%
# "step2_met": methods to control for multiple genes - "BH" (default), "ST", "BY", or "Bonferroni"

mul_corr_BPerm <- function(alltests, fastqtl_output, scheme = "BPerm1k", fdr_desired = 0.05, step2_met = "BH", samplesize) {
  cat("  ** Step 1: Permutation (beta approximation), Step 2: ", step2_met, " ** \n", sep = "")
  nthresholds <- getNominalThresholds.bperm(fastqtl_output, opt_fdr = fdr_desired, met2 = step2_met, nindiv = samplesize)
  names(nthresholds)[2] <- paste0(scheme, ".", step2_met)
  dd <- getAssociation.withThresholds(d1 = alltests[, c("snps","gene","pvalue")], d2 = nthresholds, opt_fdr = fdr_desired)
  names(dd)[3] <- paste0(scheme, ".", step2_met)
  return(list(dd, nthresholds))
}


# Multiplicity correction using permutations done by MatrixEQTL
# "alltests": output of MatrixEQTL containing as least "[snps]", "[gene]", and "[pvalue]"
# "allpermutedpval": output of manual permutations, first row is gene names and each of the resting rows is one permutation
# "ppval_gene": permutation adjusted pvalues for each gene, columns containing "[gene]", "[per.pval]"
# "fdr_desired": the desired FDR level, default is 0.05
# "step2_met": methods to control for multiple genes - "BH" (default), "ST", or "BY"

mul_corr_ManualPerm <- function(alltests, allpermutedpval, ppval_gene, fdr_desired = 0.05, step2_met = "BH") {
  cat("  ** Step 1: Permutation 1k (manually), Step 2:", step2_met, "** \n")
  # Make it - each row is one gene and each column is one permutation
  allpermutedpval <- as.data.frame(t(allpermutedpval))
  # Put gene IDs as column names
  rownames(allpermutedpval) <- allpermutedpval$V1
  allpermutedpval <- allpermutedpval[,-1]
  colnames(allpermutedpval) <- paste0("p", 1:1000)
  
  # Step 2 to correct for multiple genes
  if(step2_met == "BH") {
    ppval_gene$gene.fdr <- p.adjust(ppval_gene$per.pval, "fdr")
  } else if(step2_met == "ST") {
    ppval_gene$gene.fdr <- qvalue(ppval_gene$per.pval)$qvalues
    cat("    calculating qvalues in Step 2, pi0 =", qvalue(ppval_gene$per.pval)$pi0, "\n")
  } else if(step2_met == "BY") {
    ppval_gene$gene.fdr <- p.adjust(ppval_gene$per.pval, "BY")
  } #else if(step2_met == "Bonferroni") {
    #ppval_gene$gene.fdr <- p.adjust(ppval_gene$per.pval, "bonferroni")
  #}     # Cannot use Bonferroni - the minimum ppval is 1/1000, 1/1000*618 > 0.05, so none would pass the threshold
  
  # Permutation empirical threshold
  pthreshold <- getThreshold(ppval_gene[,c("per.pval", "gene.fdr")], opt_fdr = fdr_desired)
  cat("    The empirical p-value threshold =", pthreshold, "\n")
  
  # Calculate the nominal threshold for each gene
  # Take the integer part of a number: the ranking of the largest permuted pvalue that are smaller than the original pvalue. The next one is the nominal threshold.
  whichone <- floor(pthreshold * 1001 - 1) + 1
  nthre <- data.frame(gene = rownames(allpermutedpval), thre = 0)
  allpermutedpval <- as.matrix(allpermutedpval)
  nthre$thre <- foreach(ii = 1:nrow(nthre), .combine = "c") %do% sort(as.numeric(allpermutedpval[ii,]))[whichone]
  # Significantly associated SNP-gene pairs
  dd <-  getAssociation.withThresholds(d1 = alltests, d2 = nthre, opt_fdr = fdr_desired)
  # Check the number of identified eGenes
  if(sum(ppval_gene$gene.fdr <= 0.05) != length(unique(dd$gene)) | !all(unique(dd$gene) %in% ppval_gene[which(ppval_gene$gene.fdr <= 0.05), "gene"])) {cat("    Check the eGenes (ST)! \n") ; return(NULL)}
  
  names(dd)[3] <- paste0("Perm1k.", step2_met)
  names(nthre)[2] <- paste0("Perm1k.", step2_met)
  
  return(list(dd, nthre))
}


# Multiplicity correction using FDR procedures in Step 1; BH or ST in Step 2
# "alltests": output of MatrixEQTL containing as least "snps", "gene", and "pvalue"
# "step1_met": "ST", "BH", or "BY" to control multiple cis SNPs for each gene tested
# "fdr_desired": the desired level of gene-level FDR, default value 5%
# "step2_met": "BH" (default) or "ST" to control for multiple genes
# "mulcores": whether we use multiple cores; default is FALSE, or give it an integer as the number of cores we will use.

mul_corr_FDR_in_hierarchical <- function(alltests, step1_met, fdr_desired = 0.05, step2_met = "BH", mulcores = F) {
  cat("  ** Step 1: ", step1_met, ", Step 2: ", step2_met, "-FDR ** \n", sep = "")
  alltests <- alltests[, c("snps","gene","pvalue")]
  # List of genes
  genelist <- as.character(unique(alltests$gene))
  
  # Step 1 - correct for multiple SNPs for each gene in a loop
  # If mulcore is false (only use one core, which is the default):
  if(!mulcores) {
    alltests.fdr <- foreach(ii = 1:length(genelist), .combine = 'rbind') %do% {
      # Tests relating to this gene
      d <- alltests[which(alltests$gene == genelist[ii]),]
      if(step1_met == "ST") { d$adjpval <- qvalue(d$pvalue)$qvalues }
      if(step1_met == "BH") { d$adjpval <- p.adjust(d$pvalue, method = "BH") }
      if(step1_met == "BY") { d$adjpval <- p.adjust(d$pvalue, method = "BY") }
      #cat(ii, "\n")
      return(d)
    }
  }
  # If we choose to use multiple cores:
  if(mulcores > 0) {
    library(doMC); registerDoMC(cores = mulcores)
    alltests.fdr <- foreach(ii = 1:length(genelist), .combine = 'rbind') %dopar% {
      # Tests relating to this gene
      d <- alltests[which(alltests$gene == genelist[ii]),]
      if(step1_met == "ST") { d$adjpval <- qvalue(d$pvalue)$qvalues }
      if(step1_met == "BH") { d$adjpval <- p.adjust(d$pvalue, method = "BH") }
      if(step1_met == "BY") { d$adjpval <- p.adjust(d$pvalue, method = "BY") }
      #cat(ii, "\n")
      return(d)
    }
  }
  # The minimum p-value for each gene
  alltests.fdr <- alltests.fdr[order(alltests.fdr$pvalue),]
  min.fdr.gene <- alltests.fdr[which(!duplicated(alltests.fdr$gene)),]
  
  # Step 2 - correct for multiple genes
  dd <- getAssociation(d1 = min.fdr.gene[,c("gene","adjpval")], d2 = alltests.fdr[,c("snps","gene","pvalue","adjpval")], met2 = step2_met, opt_fdr = fdr_desired)
  names(dd)[3] <- paste0(step1_met, ".", step2_met)
  
  return(dd)
}


# Multiplicity correction using Bonferroni in Step 1; BH or ST in Step 2
# "alltests": output of MatrixEQTL containing as least "snps", "gene", and "pvalue"
# "fdr_desired": the desired level of gene-level FDR, default value 5%
# "step2_met": "BH" (default), "ST", or "Bonferroni" to control for multiple genes

mul_corr_Bonf_in_hierarchical <- function(alltests, fdr_desired = 0.05, step2_met = "BH") {
  cat("  ** Step 1: Bonferroni", ", Step 2: ", step2_met, " ** \n", sep = "")
  # The minimum p-value of each gene
  min.bon.gene <- alltests[which(!duplicated(alltests$gene)),]
  # Number of local SNPs for each gene
  gene.nSNPs <- as.data.frame(table(alltests$gene))
  names(gene.nSNPs) <- c("gene", "nSNPs")
  
  # Step 1 - Bonferroni adjust p-values
  # Best pvalues of genes
  min.bon.gene <- merge(min.bon.gene, gene.nSNPs, by = "gene")
  min.bon.gene$bpval <- min.bon.gene$pvalue * min.bon.gene$nSNPs
  min.bon.gene$bpval[which(min.bon.gene$bpval > 1)] <- 1
  # All pvalues for tests (<0.05)
  tests.bonf <- alltests[which(alltests$pvalue < 0.05),]
  tests.bonf <- merge(tests.bonf, gene.nSNPs, by = "gene")
  tests.bonf$bpval <- tests.bonf$pvalue * tests.bonf$nSNPs
  tests.bonf$bpval[which(tests.bonf$bpval > 1)] <- 1
  
  # Step 2 - correct for multiple genes using ST, BH, or Bonferroni
  dd <- getAssociation(d1 = min.bon.gene[,c("gene","bpval")], d2 = tests.bonf[,c("snps","gene","pvalue","bpval")], met2 = step2_met, opt_fdr = fdr_desired)
  names(dd)[3] <- paste0("Bonferroni.", step2_met)
  return(dd)
}

# Multiplicity correction using Bonferroni in Step 1; BH or ST in Step 2
# "alltests": output of MatrixEQTL containing as least "snps", "gene", and "pvalue"
# "fdr_desired": the desired level of gene-level FDR, default value 5%
# "step2_met": "BH" (default), "ST", or "Bonferroni" to control for multiple genes
# return both significant associations and nominal thresholds

mul_corr_Bonf_in_hierarchical_nomthr <- function(alltests, fdr_desired = 0.05, step2_met = "BH") {
  cat("  ** Step 1: Bonferroni", ", Step 2: ", step2_met, " ** \n", sep = "")
  # The minimum p-value of each gene
  min.bon.gene <- alltests[which(!duplicated(alltests$gene)),]
  # Number of local SNPs for each gene
  gene.nSNPs <- as.data.frame(table(alltests$gene))
  names(gene.nSNPs) <- c("gene", "nSNPs")
  
  # Step 1 - Bonferroni adjust p-values
  # Best pvalues of genes
  min.bon.gene <- merge(min.bon.gene, gene.nSNPs, by = "gene")
  min.bon.gene$bpval <- min.bon.gene$pvalue * min.bon.gene$nSNPs
  min.bon.gene$bpval[which(min.bon.gene$bpval > 1)] <- 1
  
  # Step 2 - correct for multiple genes
  if(step2_met == "BH") {
    min.bon.gene$gene.fdr <- p.adjust(min.bon.gene$bpval, "fdr")
  } else if(step2_met == "ST") {
    min.bon.gene$gene.fdr <- qvalue(min.bon.gene$bpval)$qvalues
    cat("    calculating qvalues in Step 2, pi0 =", qvalue(min.bon.gene$bpval)$pi0, "\n")
  } else if(step2_met == "BY") {
    min.bon.gene$gene.fdr <- p.adjust(min.bon.gene$bpval, "BY")
  } else if(step2_met == "Bonferroni") {
    min.bon.gene$gene.fdr <- p.adjust(min.bon.gene$bpval, "bonferroni")
  }
  # Empirical p-value threshold
  threshold <- getThreshold(min.bon.gene[,c("bpval","gene.fdr")], opt_fdr = fdr_desired)
  cat("    The empirical p-value threshold =", threshold, "\n")
  # Nominal pvalue threshold for each gene
  nthre.Bonf <- data.frame(gene = min.bon.gene$gene, Bonf = threshold / min.bon.gene$nSNPs)
  names(nthre.Bonf)[2] <- paste0("Bonferroni.", step2_met)
  # Get the significant pairs based on the nominal thresholds
  dd <- getAssociation.withThresholds(d1 = alltests[, c("snps","gene","pvalue")], d2 = nthre.Bonf, opt_fdr = fdr_desired)
  names(dd)[3] <- paste0("Bonferroni.", step2_met)
  
  return(list(dd,nthre.Bonf))
}



# ----- Functions used for multiplicity correction -----
# Calculate the empirical p-value threshold at the 0.05 significant level
# "d": a table with 2 columns - "empirical p-values" and "after correction for multiple genes"; order sensitive
# "opt_fdr": the desired level of FDR of eGenes, default value is 5%
getThreshold <- function(d, opt_fdr = 0.05) {
  set0 <- d[which(d[,2] <= opt_fdr),]
  set1 <- d[which(d[,2] > opt_fdr),]
  threshold <- (sort(set1[,1])[1] - sort(-1.0 * set0[,1])[1]) / 2
  return(threshold)
}


# Get significant SNP-gene pairs - Step 2 and 3
# "d1": 2 columns - "[gene]" and "min empirical p-values"; order sensitive
# "d2": 4 columns - "[snps]", "[gene]", "[pvalue]" (original pvalues), and "empirical p-values"; first 3 columns should have those colnames, the 4th column should be empirical p-values.
# "met2": method used in Step 2 - "ST", "BH", "BY", or "Bonferroni"
getAssociation <- function(d1, d2, met2, opt_fdr = 0.05) {
  if(met2 == "ST") {
    # Calculate q-values
    qval <- qvalue(d1[,2])
    cat("    calculating qvalues in Step 2, pi0 =", qval$pi0, "\n")
    d1$gene.qval <- qval$qvalues
  } else if(met2 == "BH") {
    d1$gene.fdr <- p.adjust(p = d1[,2], method = "BH")
  } else if(met2 == "Bonferroni") {
    d1$gene.bpavl <- p.adjust(d1[,2], method = "bonferroni")
  } else if (met2 == "BY") {
    d1$gene.fdr <- p.adjust(p = d1[,2], method = "BY")
  }
  # Get the empirical p-value threshold
  threshold <- getThreshold(d1[,2:3], opt_fdr = opt_fdr)
  cat("    The empirical p-value threshold =", threshold, "\n")
  # Get the significant SNP-gene associations
  d <- d2[which(d2[,4] <= threshold), c("snps", "gene", "pvalue")]
  d <- d[order(d$pvalue),]
  # Check the number of eGenes
  if(length(unique(d$gene)) != sum(d1[,3] <= 0.05) | !all(unique(d$gene) %in% d1[which(d1[,3] <= 0.05), "gene"])) {
    cat("    Check the eGenes! \n")  }
  return(d)
}


# Get significant SNP-gene pairs based on nominal p-value thresholds
# "d1": 3 columns - "snps", "[gene]" and pvalue"
# "d2": nominal thresholds, 2 columns - "[gene]" and "nominal p-value thresholds"
getAssociation.withThresholds <- function(d1, d2, opt_fdr = 0.05) {
  d1 <- d1[which(d1$pvalue < opt_fdr),]
  d1 <- merge(d1, d2, by = "gene")
  d1 <- d1[which(d1$pvalue <= d1[,4]),]
  d1 <- d1[order(d1$pvalue), -4]
  return(d1)
}


# Calculate the nominal p-value threshold for each gene based on beta approximation permutations. (from Olivier Delaneau, FastQTL)
# "perm": the output of FastQTL
# "opt_fdr": the desired FDR level, default value is 0.05
# "met2": the method used in Step 2 - "ST", "BH", "BY", or "Bonferroni"
# "nindiv": sample size
# return a table, first column is gene, second is nominal threshold
getNominalThresholds.bperm <- function(perm, opt_fdr = 0.05, met2, nindiv) {
  cat("    Correlation between Beta approx. and Empirical p-values =", round(cor(perm[,10], perm[,11]), 4), "\n")
  # Correct for multiple genes (Step2)
  if(met2 == "ST") {
    perm$qval <- qvalue(perm[,11])$qvalues
  } else if(met2 == "BH") {
    perm$bhfdr <- p.adjust(perm[,11], method = "fdr")
  } else if(met2 == "BY") {
    perm$byfdr <- p.adjust(perm[,11], method = "BY")
  } else if(met2 == "Bonferroni") {
    perm$byfdr <- p.adjust(perm[,11], method = "bonferroni")
  }
  # Empirical pvalue thresholds
  pthreshold <- getThreshold(perm[,11:12], opt_fdr)
  cat("    Permutation threshold at the 0.05 q-value level =", pthreshold, "\n")
  # Calculate nominal pvalue thresholds based on empirical pvalue threshold
  pval0 <- qbeta(pthreshold, perm[,3], perm[,4], ncp = 0, lower.tail = TRUE, log.p = FALSE)
  # quantile function for the F distribution; the 5th: Dummy
  test0 = qf(pval0, 1, perm[,5], ncp = 0, lower.tail = FALSE, log.p = FALSE)
  corr0 = sqrt(test0 / (perm[,5] + test0))
  test1 = nindiv * corr0 * corr0 / (1 - corr0 * corr0)
  pval1 = pf(test1, 1, nindiv, ncp = 0, lower.tail = FALSE, log.p = FALSE)
  cat("  * pval0 = ", mean(pval0), " +/- ", sd(pval0), "\n")
  cat("  * test0 = ", mean(test0), " +/- ", sd(test0), "\n")
  cat("  * corr0 = ", mean(corr0), " +/- ", sd(corr0), "\n")
  cat("  * test1 = ", mean(test1), " +/- ", sd(test1), "\n")
  cat("  * pval1 = ", mean(pval1), " +/- ", sd(pval1), "\n")
  # Write the thresholds (after correction, pval1 is more conservative than pval0)
  nthresholds <- data.frame(gene = perm[,1], nthresholds = pval1)
  return(nthresholds)
}


# -------------------- 3. Calculate sensitivity and FDR --------------------
# Calculate TPR - how many true eGenes could be identified
# "resultsdata": every method has 4 columns; R2_causalSNP; each row is one scenairo or one replicate
# "methodlist": list of method names, or "all"
calculate_TPR <- function(resultsdata, methodlist = "all") {
  if(methodlist == "all") {
    jj <- which(names(resultsdata) == "n.CP") + 1
    methodlist <- unlist(strsplit(names(resultsdata)[seq(jj, ncol(resultsdata)-1, by = 4)], split = "_nasso"))
  }
  # TPR
  TPR_eGene <- resultsdata[, paste0(methodlist, "_TPgene")]/200
  colnames(TPR_eGene) <- methodlist
  return(TPR_eGene)
}


# Calculate FDR - how many significant results are actually false positives
# "resultsdata": every method has 4 columns; R2_causalSNP; each row is one scenairo or one replicate
# "methodlist": list of method names, or "all"
calculate_FDR <- function(resultsdata, methodlist = "all") {
  if(methodlist == "all") {
    jj <- which(names(resultsdata) == "n.CP") + 1
    methodlist <- unlist(strsplit(names(resultsdata)[seq(jj, ncol(resultsdata)-1, by = 4)], split = "_nasso"))
  }
  # FDR
  FDR_eGene <- (resultsdata[, paste0(methodlist, "_ngene")] - resultsdata[, paste0(methodlist, "_TPgene")]) / resultsdata[, paste0(methodlist, "_ngene")]
  for(ii in 1: ncol(FDR_eGene)) {FDR_eGene[,ii][is.na(FDR_eGene[,ii])] <- 0}
  colnames(FDR_eGene) <- methodlist
  return(FDR_eGene)
}


# Calculate the average values (FDR/TPR)
# "n_reps": number of replicates, default is 10
# "data_reps": each row indicates one rep for each scenario. Columns contain "[Rep]", "[SS]", and "[MAF]"; "[effectsize]" may also exist. Each of the resting columns indicate the performance of one method.
calculate_average_values <- function(data_reps, n_reps) {
  # Whether a fixed eQTL effect size was simulated in each scenairo
  if("effectsize" %in% names(data_reps)) {
    # Put replicates together
    data_reps <- data_reps[order(data_reps$SS, data_reps$MAF, data_reps$effectsize, data_reps$Rep),]
    # Unique scenarios
    ave_data <- unique(data_reps[,c("SS","MAF","effectsize")])
    # Take average
    for(cc in 5:ncol(data_reps)) {
      ave_data$newc <- sapply(1:nrow(ave_data), function(rr) { mean(data_reps[1:n_reps+(rr-1)*n_reps, cc]) })
      names(ave_data)[ncol(ave_data)] <- names(data_reps)[cc]
    } # End of loop for all methods
  } else {
    # Put replicates together
    data_reps <- data_reps[order(data_reps$SS, data_reps$MAF, data_reps$Rep),]
    # Unique scenarios
    ave_data <- unique(data_reps[,c("SS","MAF")])
    # Take average
    for(cc in 4:ncol(data_reps)) {
      ave_data$newc <- sapply(1:nrow(ave_data), function(rr) { mean(data_reps[1:n_reps+(rr-1)*n_reps, cc]) })
      names(ave_data)[ncol(ave_data)] <- names(data_reps)[cc]
    } # End of loop for all methods
  }
  
  return(ave_data)
}

  
# Sum the results of replicates for each scenario and then calculate TPR/FDR based on all discoveries from all replicates 
# "resultsdata": each row is one replicate of a scenario; first 4 columns: [Rep], [SS], [MAF], and [n.CP] and every method has four columns.
# "methodlist": default is "all", all the method in the table. Or give a list of methods.
# "FDRorTPR": calculate "FDR" or "TPR".

calculateTPRFDR_allreps <- function(resultsdata, methodlist = "all", FDRorTPR) {
  # 1. Sum the number of TP and FP for each scenario
  resultsdata <- resultsdata[order(resultsdata$SS, resultsdata$MAF, resultsdata$Rep),]
  sumallreps <- unique(resultsdata[,c("SS","MAF")])

  # Which methods
  if(methodlist == "all") {
    jj <- which(names(resultsdata) == "n.CP") + 1
    methodlist <- unlist(strsplit(names(resultsdata)[seq(jj, ncol(resultsdata))], split = "_"))
    methodlist <- unique(setdiff(methodlist, c("nasso","ngene","TP","TPgene")))
  }

  # Keep the column method_ngene and method_TPgene
  for(mm in methodlist) {
    for(whichone in c("_ngene","_TPgene")) {
      sumallreps$newcol <- sapply(1:nrow(sumallreps), function(ii) {
        sum(resultsdata[which(resultsdata$SS == sumallreps$SS[ii] & resultsdata$MAF == sumallreps$MAF[ii]), paste0(mm, whichone)])
      })
      names(sumallreps)[ncol(sumallreps)] <- paste0(mm, whichone)
    } # End of loop 2
  } # End of loop for methods
  
  # Calculate TPR or FDR
  if(FDRorTPR == "TPR") {
    returndd <- sumallreps[,paste0(methodlist, "_TPgene")] / (200*nrow(resultsdata)/36)
  } else if(FDRorTPR == "FDR") {
    returndd <- (sumallreps[,paste0(methodlist, "_ngene")] - sumallreps[,paste0(methodlist, "_TPgene")]) / sumallreps[,paste0(methodlist, "_ngene")]
  }
  colnames(returndd) <- methodlist
  returndd <- cbind(sumallreps[,c("SS","MAF")], returndd)

  return(returndd)
}


# -------------------- 4. Calculate average/sum for each scenario --------------------
# Calculate the average for each scenario, usually across all 100 replicates, after removing NAs if any exists.
# "data_reps": data for each scenario occupy a few rows. First Column: [Rep], contains column [SS] and [MAF], maybe [effectsize], the resting columns need to be taken average.
# "scenario_parameters": a vector indicating which parameters each scenario has. Default is c("SS","MAF"), or c("SS","MAF","effectsize")

calculate_ave_ignoreNAs <- function(data_reps, scenario_parameters = c("SS","MAF")) {
  # Unique scenarios (also the return table)
  ave_data <- unique(data_reps[,scenario_parameters])
  
  # Take average for each of the remaining columns
  if(identical(scenario_parameters, c("SS","MAF"))) {
    for(cc in setdiff(names(data_reps), c("Rep", scenario_parameters))) {
      ave_data$newc <- sapply(1:nrow(ave_data), function(rr) {
        dd <- data_reps[which(data_reps$SS == ave_data$SS[rr] & data_reps$MAF == ave_data$MAF[rr]), cc]
        # Report which scenarios contain NAs
        if(sum(is.na(dd))>0) {cat(cc,"SS",ave_data$SS[rr],"MAF",ave_data$MAF[rr],":",sum(is.na(dd)),"NA out of",length(dd),"variables \n")}
        dd <- dd[which(complete.cases(dd))]
        return(mean(dd))
      })
      names(ave_data)[ncol(ave_data)] <- cc
    } # End of loop for columns
  } else if(identical(scenario_parameters, c("SS","MAF","effectsize"))) {
    for(cc in setdiff(names(data_reps), c("Rep", scenario_parameters))) {
      ave_data$newc <- sapply(1:nrow(ave_data), function(rr) {
        dd <- data_reps[which(data_reps$SS == ave_data$SS[rr] & data_reps$MAF == ave_data$MAF[rr] & data_reps$effectsize == ave_data$effectsize[rr]), cc]
        # Report which scenarios contain NAs
        if(sum(is.na(dd))>0) {cat(cc,"SS",ave_data$SS[rr],"MAF",ave_data$MAF[rr],"effect size",ave_data$effectsize[rr],":",sum(is.na(dd)),"NA out of",length(dd),"variables \n")}
        dd <- dd[which(complete.cases(dd))]
        return(mean(dd))
      })
      names(ave_data)[ncol(ave_data)] <- cc
    } # End of loop for columns
  }
  
  return(ave_data)
}


# Calculate the sum for each scenario, usually across all 100 replicates, after removing NAs if any exists.
# "data_reps": data for each scenario occupy a few rows. First Column: [Rep], contains column [SS] and [MAF], maybe [effectsize], the resting columns need to be taken average.
# "scenario_parameters": a vector indicating which parameters each scenario has. Default is c("SS","MAF"), or c("SS","MAF","effectsize")

calculate_sum_ignoreNAs <- function(data_reps, scenario_parameters = c("SS","MAF")) {
  # Unique scenarios (also the return table)
  sum_data <- unique(data_reps[,scenario_parameters])
  
  # Take sum for each of the remaining columns
  if(identical(scenario_parameters, c("SS","MAF"))) {
    for(cc in setdiff(names(data_reps), c("Rep", scenario_parameters))) {
      sum_data$newc <- sapply(1:nrow(sum_data), function(rr) {
        dd <- data_reps[which(data_reps$SS == sum_data$SS[rr] & data_reps$MAF == sum_data$MAF[rr]), cc]
        # Report which scenarios contain NAs
        if(sum(is.na(dd))>0) {cat(cc,"SS",sum_data$SS[rr],"MAF",sum_data$MAF[rr],":",sum(is.na(dd)),"NA out of",length(dd),"variables \n")}
        dd <- dd[which(complete.cases(dd))]
        return(sum(dd))
      })
      names(sum_data)[ncol(sum_data)] <- cc
    } # End of loop for columns
  } else if(identical(scenario_parameters, c("SS","MAF","effectsize"))) {
    for(cc in setdiff(names(data_reps), c("Rep", scenario_parameters))) {
      sum_data$newc <- sapply(1:nrow(sum_data), function(rr) {
        dd <- data_reps[which(data_reps$SS == sum_data$SS[rr] & data_reps$MAF == sum_data$MAF[rr] & data_reps$effectsize == sum_data$effectsize[rr]), cc]
        # Report which scenarios contain NAs
        if(sum(is.na(dd))>0) {cat(cc,"SS",sum_data$SS[rr],"MAF",sum_data$MAF[rr],":",sum(is.na(dd)),"NA out of",length(dd),"variables \n")}
        dd <- dd[which(complete.cases(dd))]
        return(sum(dd))
      })
      names(sum_data)[ncol(sum_data)] <- cc
    } # End of loop for columns
  }
  
  return(sum_data)
}








