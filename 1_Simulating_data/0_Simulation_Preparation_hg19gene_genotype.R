# 2018-04-11
library(foreach)
library(doMC); registerDoMC(cores = 6)
library(ggplot2)
options(stringsAsFactors = F)

#----------------------------------------------------------------------
# Refseq gene annotation hg19/GRCh37
#----------------------------------------------------------------------
# Data downloaded
setwd("~/Simulation_eQTL_replicates/hg19_gene_annotation")

# Refseq gene data (There was a hash# at the beginning of the column names, removed.); 1284 rows
refgene <- read.table("RefSeq_Genes_chr22_hg19_GRCh37.txt", header = T, sep = "\t")

# Check data are on chr 22
if(!all(refgene$chrom == "chr22")) { cat("Not all genes are on chr22!\n") }

# Gene location
geneloc <- refgene[, c("name2","strand","txStart","txEnd")]
# Remove duplicated gene names
geneloc <- geneloc[!duplicated(geneloc$name2),]

# TSS
geneloc$TSS <- geneloc$txStart
geneloc$TSS[which(geneloc$strand == "-")] <- geneloc$txEnd[which(geneloc$strand == "-")]
# Unique TSS (618 genes)
geneloc <- geneloc[!duplicated(geneloc$TSS),]

# Write
write.table(geneloc, "refseq_gene_unique_hg19.txt", row.names = F, quote = F, sep = "\t")

# gene location file for the Matrix eQTL
d <- data.frame(geneid = geneloc$name2, chr = 22, left = geneloc$TSS, right = geneloc$TSS)
write.table(d, "geneloc_for_ME.txt", quote = F, sep = "\t", row.names = F)


#----------------------------------------------------------------------
# Generate input files of HAPGEN2
#----------------------------------------------------------------------
# Input files - .hap .legend and .map
# Convert 1000 Genomes project phase 3 VCF data to gen/legend/sample file.
# Download HapMap phase2 b37 file -> genetic map 
setwd("~/Simulation_eQTL_replicates/1000GP_phase3_hg19_GRCh37/")

# Read sample information
popinfo <- read.table("1000Genomes_phase3/integrated_call_samples_v3.20130502.ALL.panel", header = T)

# Subset CEU samples; 99 samples
CEU <- popinfo[which(popinfo$pop == "CEU"), c("sample")]
write.table(CEU, "CEU_samples_99.txt", quote = F, sep = "\t", col.names = F, row.names = F)

# Subset FIN samples; 99 samples
FIN <- popinfo[which(popinfo$pop == "FIN"), c("sample")]
write.table(FIN, "FIN_samples_99.txt", quote = F, sep = "\t", col.names = F, row.names = F)

# Use vcftools to convert vcf to IMPUTE input format (hap and legend files); IMPUTE requires phased data and only bi-allelic sites are included. 
# Data of CEU individuals; fitering - MAF and HWE
system("vcftools --gzvcf 1000Genomes_phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep CEU_samples_99.txt --IMPUTE --maf 0.005 --hwe 1e-6 --out CEU_chr22_filtered")
# Remove sample list file, which is identical to .hap.indv output.
file.remove("CEU_samples_99.txt")

# Data of FIN individuals; fitering - MAF and HWE
system("vcftools --gzvcf 1000Genomes_phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep FIN_samples_99.txt --IMPUTE --maf 0.005 --hwe 1e-6 --out FIN_chr22_filtered")
file.remove("FIN_samples_99.txt")
# Compress some files
system("gzip *_chr22_filtered.impute.hap && gzip *_chr22_filtered.impute.legend")

# hap/legend files contain SNP IDs that appear more than once.
foreach(pop = c("CEU","FIN")) %dopar% {
  chr22.hap <- read.table(paste(pop, "_chr22_filtered.impute.hap.gz", sep = ""), header = F)
  chr22.legend <- read.table(paste(pop, "_chr22_filtered.impute.legend.gz", sep = ""), header = T)
  
  # Remove SNPs that have multiple rows: 
  #(1) some of them have different alternative alleles - same pos and ID
  #(2) some of SNPs were mapped to the same location - different IDs but same pos).
  SNPskeep <- setdiff( 1:nrow(chr22.legend), which(chr22.legend$pos %in% chr22.legend$pos[which(duplicated(chr22.legend$pos))]) )
  chr22.hap <- chr22.hap[SNPskeep,]
  chr22.legend <- chr22.legend[SNPskeep,]
  write.table(chr22.hap, file = paste(pop, "_chr22.hap", sep = ""), row.names = F, col.names = F, quote = F, sep = " ")
  write.table(chr22.legend, file = paste(pop, "_chr22.legend", sep = ""), row.names = F, quote = F, sep = " ")
  
  return(NULL)
}



#----------------------------------------------------------------------
# Sanity check
#----------------------------------------------------------------------
# 1. MAF distribution
setwd("~/Simulation_eQTL_replicates/1000GP_phase3_hg19_GRCh37")

# CEU genotypes
CEU.hap <- read.table("CEU_chr22.hap", header = F)
CEU.leg <- read.table("CEU_chr22.legend", header = T)
# Calculate MAFs
CEU.leg$MAF <- apply(CEU.hap, 1, sum)/198
CEU.leg$MAF[which(CEU.leg$MAF > 0.5)] <- 1-CEU.leg$MAF[which(CEU.leg$MAF > 0.5)]

# FIN genotypes
FIN.hap <- read.table("FIN_chr22.hap", header = F)
FIN.leg <- read.table("FIN_chr22.legend", header = T)
# Calculate MAFs
FIN.leg$MAF <- apply(FIN.hap, 1, sum)/198
FIN.leg$MAF[which(FIN.leg$MAF > 0.5)] <- 1-FIN.leg$MAF[which(FIN.leg$MAF > 0.5)]

# MAF density
MAFs <- rbind( data.frame(MAF = CEU.leg$MAF, pop = "CEU"),
               data.frame(MAF = FIN.leg$MAF, pop = "FIN") )
pden <- ggplot(MAFs, aes(MAF)) + geom_density(aes(color = pop))
ggsave(filename = "MAFs_density_CEU.pdf", plot = pden, width = 6, height = 5)


#----- number of bases of alleles -----
library(stringr)
CEU.leg$lalleles <- sapply(1:nrow(CEU.leg), function(ii) { str_length(CEU.leg$allele0[ii])+str_length(CEU.leg$allele1[ii]) })
FIN.leg$lalleles <- sapply(1:nrow(FIN.leg), function(ii) { str_length(FIN.leg$allele0[ii])+str_length(FIN.leg$allele1[ii]) })

table(CEU.leg$lalleles)
table(FIN.leg$lalleles)

# Percentage of length of 2 in all variants
table(CEU.leg$lalleles)[1]/nrow(CEU.leg)
table(FIN.leg$lalleles)[1]/nrow(FIN.leg)
  # CEU: 175105 (90.97%); FIN: 154600 (90.17%).




