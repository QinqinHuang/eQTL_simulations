############################################
# Which Replicates:
whichreps <- 11:100
howmanyreps <- length(whichreps)
############################################

library(foreach)
library(doMC); registerDoMC(cores = 40)
options(stringsAsFactors = F)
source("~/Simulation_eQTL_replicates/Scripts/001_all_functions_for_eQTL_simulations.R")


#----------------------------------------------------------------------
# Simulate genotype data for 6 sample sizes 
#----------------------------------------------------------------------

nothing <- foreach(IdxRep = rep(whichreps, each = 6), SS = rep(c(100,200,500,1000,2000,5000), howmanyreps)) %dopar% {

  # Create a folder for current replicate/sample size and a folder under that for genotypes
  geno_dir <- paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/Genotype/") 
  dir.create(geno_dir, recursive = T)
  # Working directory
  setwd(geno_dir)
  
  # Use HAPGEN2 to simulate genotypes
  system(paste("hapgen2 -h ~/Simulation_eQTL_replicates/1000GP_phase3_hg19_GRCh37/FIN_chr22.hap -l ~/Simulation_eQTL_replicates/1000GP_phase3_hg19_GRCh37/FIN_chr22.legend -m ~/Simulation_eQTL_replicates/1000GP_phase3_hg19_GRCh37/HapMap_phase2_b37_genetic_map/genetic_map_chr22_combined_b37.txt -o hapgen2_out -dl 16050840 1 2 4 -n", SS, 1))
  
  # Use Plink to get VCF genotype file (for FastQTL)
  system("plink --gen hapgen2_out.controls.gen --sample hapgen2_out.controls.sample --oxford-single-chr 22 --out chr22.genotype")
  system("plink --bfile chr22.genotype --maf 0.005 --hwe 1e-6 --make-bed --out filtered.genotype")
  system("plink --bfile filtered.genotype --recode vcf --out filtered")  
  system("bgzip filtered.vcf && tabix -p vcf filtered.vcf.gz")
  
  # Use Plink to get MatrixEQTL input genotype file
  system("plink --bfile filtered.genotype --recode A-transpose --out filtered.data")
  system("cut -f 2,7- filtered.data.traw > SNP.txt && cut -f 1,2,4 filtered.data.traw > snpsloc.txt")
  snpsloc <- read.table("snpsloc.txt", header = T)
  snpsloc <- snpsloc[,c(2,1,3)]
  write.table(snpsloc, "snpsloc.txt", quote = F, sep = "\t", row.names = F)
  
  # Calculate MAFs
  system("plink --bfile filtered.genotype --freq --out MAF_filtered")
  system("cat MAF_filtered.frq | awk '{print$2,$5}' > SNPs_MAF.txt")

  # Remove some files
  system("rm chr22.genotype.* filtered.log filtered.data.* *.nosex MAF_filtered.* hapgen2_out.*")
  
  # Use Plink to do LD pruning (r2 threshold = 0.3) -- for determining true eQTLs
  system("plink --bfile filtered.genotype --indep-pairwise 1000kb 5 0.3 --out LD_pruned_window1000kb")
  system("rm LD_pruned_window1000kb.prune.out LD_pruned_window1000kb.log LD_pruned_window1000kb.nosex")
  
  return(NULL)
} # End of loop for 45 replicates and 6 sample sizes (45*6 = 270)



#----------------------------------------------------------------------
# Determine true eQTLs
#----------------------------------------------------------------------
#----- My function -----
# Log file when determining true causal eSNPs, under current folder, one file for each simulation
putlog <- function(newline, IdxRep) {
  write(newline, file = paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/LogFIles/Rep", IdxRep, "_log_Choose_true_eQTLs.txt"), append = T)
}


# Read gene location file
geneloc <- read.table("~/Simulation_eQTL_replicates/hg19_gene_annotation/geneloc_for_ME.txt", header = T)

# Determine true eQTLs, run in parallel for each replicate
nothing <- foreach(IdxRep = whichreps) %dopar% {
  # 36 scenarios in sequence
  for(SS in c(100,200,500,1000,2000,5000)) {
    putlog(paste0("# ----- Rep = ", IdxRep, " & Sample size = ", SS, " -----"), IdxRep = IdxRep)
    
    # Directory of this sample size
    setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS))
    # Read pruned list of SNPs
    pruned_snps <- read.table("Genotype/LD_pruned_window1000kb.prune.in", header = F)
        # Rep1-10 doesnt have this file, we used 
        # pruned_snps <- read.table("Genotype/LD_pruned_03/LD_pruned_window1000kb.prune.in", header = T)
    # MAFs
    pruned_snps_MAF <- read.table("Genotype/SNPs_MAF.txt", header = T)
    pruned_snps_MAF <- pruned_snps_MAF[which(pruned_snps_MAF$SNP %in% pruned_snps$V1),]
    # SNP location
    snpspos <- read.table("Genotype/snpsloc.txt", header = T)
    pruned_snps_MAF <- merge(snpspos[,-2], pruned_snps_MAF, by = "SNP")
    # Number of SNPs after LD pruning
    putlog(paste0("   Number of SNPs after LD pruning 0.3: ", nrow(pruned_snps_MAF)), IdxRep = IdxRep)
    
    # Loop for 6 MAFs
    for(MAF in c(50,25,10,5,1,0.5)) {
      # Create a new folder for current MAF
      dir.create(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
      setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
      
      #---------- Decide MAF boundaries ----------
      # Start from exact the value
      b1 <- MAF/100
      b2 <- MAF/100
      while(b1 > 0) {
        candidate <- which(pruned_snps_MAF$MAF >= b1 & pruned_snps_MAF$MAF <= b2)
        # If MAF is 50%, change the boundary
        if(MAF == 50) { b1 <- 0.49 }
        # If there are enough candidate SNPs
        if(length(candidate) >= 300) {break}
        # When MAF is 0.5%
        if(MAF == 0.5) {
          if(SS %in% c(100,200) | b2>=0.006) {break}
        }
        # When MAF is 1%
        if(MAF == 1) {
          if(SS == 100) {break}
        }
        b1 <- max((b1 - 0.5/SS), 0.005)
        b2 <- min((b2 + 0.5/SS), 0.5)
      }
      putlog(paste0("   # SS = ", SS, " MAF = ", MAF), IdxRep = IdxRep)
      putlog(paste("      Boundary:", b1, b2), IdxRep = IdxRep)
      #-----------------------------------------------------
      
      # Candidate SNPs
      pool <- pruned_snps_MAF[which(pruned_snps_MAF$MAF >= b1 & pruned_snps_MAF$MAF <= b2), ]
      putlog(paste0("      length of candidates: ", nrow(pool)), IdxRep = IdxRep)
      rownames(pool) <- 1:nrow(pool)
      
      # How many genes within the ±1Mb window:
      pool$n.cisgene <- sapply(1:nrow(pool), function(ii) {
        sum(geneloc$right < pool$POS[ii]+1e6-1 & geneloc$left > pool$POS[ii]-1e6+1)
      })
      # Distance from the next SNP
      pool$dist <- NA
      pool$dist[1:(nrow(pool)-1)] <- sapply(1:(nrow(pool)-1), function(ii) {
        pool$POS[ii+1] - pool$POS[ii]
      })
      # Remove SNPs that are too close
      if(length(which(pool$dist > 100)) >= 250) {
        pool <- pool[which(pool$dist > 100),]
        putlog(paste("      length of candidates after removing distance <= 100:", nrow(pool)), IdxRep = IdxRep)
      } else if(length(which(pool$dist > 10)) >= 250) {
        pool <- pool[which(pool$dist > 10),]
        putlog(paste("      length of candidates after removing distance <= 10:", nrow(pool)), IdxRep = IdxRep)
      } 
      
      # Randomly choose true associations
      while(TRUE) {
        # Subsample 200 eSNPs, or keep all eSNPs if there're less than 200.
        if(nrow(pool) >= 200) {
          eQTLsnp <- pool[sample(1:nrow(pool), 200), -5]
        } else {
          eQTLsnp <- pool[, -5]
        }
        # Order by number of genes in cis
        eQTLsnp <- eQTLsnp[order(eQTLsnp$n.cisgene),]
        # Choose one eGene for each true eQTL at random
        eQTLsnp$eGene <- ""
        for(i in 1:nrow(eQTLsnp)) {
          eQTLsnp$eGene[i] <- try( sample(setdiff(geneloc[which(geneloc$right < eQTLsnp$POS[i]+1000000-5 & geneloc$left > eQTLsnp$POS[i]-1000000+5), "geneid"], eQTLsnp$eGene[1:(i-1)]), 1) )
        }
        if(sum(eQTLsnp$eGene %in% geneloc$geneid) == nrow(eQTLsnp) & length(unique(eQTLsnp$eGene)) == nrow(eQTLsnp)) { break }
      }
      # TSS of eGenes
      eQTLsnp <- merge(eQTLsnp, geneloc[,c(1,3)], by.x = "eGene", by.y = "geneid")
      colnames(eQTLsnp)[ncol(eQTLsnp)] <- "TSS"
      eQTLsnp <- eQTLsnp[,c(2,3,4,5,1,6)]
      eQTLsnp <- eQTLsnp[order(eQTLsnp$POS),]
      rownames(eQTLsnp) <- 1:nrow(eQTLsnp)
      
      putlog(paste("   Number of true eSNPs:", nrow(eQTLsnp)), IdxRep = IdxRep)
      write.table(eQTLsnp, "True_eQTLs.txt", quote = F, sep = "\t", row.names = F)
      rm(pool)
      
    } # End of 6 MAFs
    
  } # End of 6 sample sizes
  
  return(NULL)
} # End of 45 replicates


#----------------------------------------------------------------------
# LD structure around the true causal eSNPs
#----------------------------------------------------------------------
# All 36 scenarios
nothing <- foreach(IdxRep = rep(whichreps, each = 36), SS = rep(rep(c(5000,2000,1000,500,200,100), each = 6), howmanyreps), MAF = rep(c(50,25,10,5,1,0.5), 6*howmanyreps) ) %dopar% {
  # Directory of current scenario
  setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
  
  # Read the true eQTLs
  eQTLsnp <- read.table("True_eQTLs.txt", header = T)
  
  # For each causal eSNP, calcualte the LD measure r2 for all SNPs within 2Mb from it and keep those with r2 > 0.01
  system("tail -n +2 True_eQTLs.txt | cut -f 1 > eSNP_list")
  #system("plink --bfile ../Genotype/filtered.genotype --r2 --ld-window 100000 --ld-window-kb 2000 --ld-window-r2 0.5 --ld-snp-list eSNP_list --out LD_eSNP_r2_05_window_2Mb")
  system("plink --bfile ../Genotype/filtered.genotype --r2 --ld-window 100000 --ld-window-kb 2000 --ld-window-r2 0.01 --ld-snp-list eSNP_list --out LD_eSNP_r2_001_window_2Mb")
  
  return(NULL)
}


#----------------------------------------------------------------------
# Simulate expression data
#----------------------------------------------------------------------

nothing <- foreach(IdxRep = rep(whichreps, each = 6), SS = rep(c(100,200,500,1000,2000,5000), howmanyreps)) %dopar% {
  
  # Read genotype data of the currect sample size
  setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS))
  sim_genotype <- read.table("./Genotype/SNP.txt", header = T)
  
  # Loop for 6 MAFs
  for(MAF in c(50,25,10,5,1,0.5)) {
    # Directory of current scenario
    setwd(paste0("~/Simulation_eQTL_replicates/",whichfolder(IdxRep),"/rep", IdxRep, "/SS", SS, "/SS", SS, "maf", MAF))
    
    # Read true causal eSNPs
    eQTLsnp <- read.table("True_eQTLs.txt", header = T)
    
    # eQTL effect sizes from a beta distribution
    eQTLsnp$sim_beta <- rgamma(nrow(eQTLsnp), shape = 3.508, rate = 6.257)
    
    # Expression data
    exprdata <- simulate_exp_data(geno_data = sim_genotype, gene_pos = geneloc, true_eQTLs_eGenes = eQTLsnp, needfastqtl = T)
    write.table(eQTLsnp, file = "True_eQTLs.txt", row.names = F, quote = F, sep = "\t")
    write.table(exprdata, "gene_expression.txt", quote = F, sep = "\t")
    
  } # End of 6 MAFs
  
  return(NULL)
}








# For Rep1-10, there is "LD_pruned_03" folder under each "Genotype" folder.
nothing <- foreach(IdxRep = rep(whichreps, each = 6), SS = rep(c(100,200,500,1000,2000,5000), howmanyreps)) %dopar% {
  
  # Use Plink to do LD pruning using r2 threshold 0.3 and 0.1
  dir.create("LD_pruned_03")
  dir.create("LD_pruned_01")
  for(r2cutoff in c(0.3, 0.1)) {
    setwd(paste("~/Simulation_eQTL_replicates/simulation_results/rep", IdxRep, "/SS", SS, "/Genotype/LD_pruned_0", r2cutoff*10, sep = ""))
    
    # LD pruning, window size = 1Mb
    system(paste("plink --bfile ../filtered.genotype --indep-pairwise 1000kb 5 ", r2cutoff, " --out LD_pruned_window1000kb", sep = ""))
    # Subset genotypes of these SNPs
    system("plink --bfile ../filtered.genotype --extract LD_pruned_window1000kb.prune.in --make-bed --out LD_pruned_geno")
  }
} 







