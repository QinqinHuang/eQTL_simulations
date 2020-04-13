library(reshape2)
library(ggplot2)
library(RColorBrewer)
options(stringsAsFactors = FALSE)

## Multiple testing correction method ##
localcormet <- "eigenMT"
globalcormet <- "BH"
# Use this method to determine TP eGenes.

# My Functions 
# 1. MAF groups
MAFasgourps <- function(dd, addMAF = F, decreasing = F) {
  if(!addMAF & !decreasing) { return(factor(paste0(dd,"%"), levels = rev(c("50%","25%","10%","5%","1%","0.5%")))) }
  if(!addMAF & decreasing) { return(factor(paste0(dd,"%"), levels = c("50%","25%","10%","5%","1%","0.5%"))) }
  if(addMAF & !decreasing) { return(factor(paste0("MAF: ",dd,"%"), levels = rev(paste0("MAF: ", c("50%","25%","10%","5%","1%","0.5%"))))) }
  if(addMAF & decreasing) { return(factor(paste0("MAF: ",dd,"%"), levels = paste0("MAF: ", c("50%","25%","10%","5%","1%","0.5%")))) }
}

# 2. effect size groups
ESasgourps <- function(dd, addES = F, decreasing = F) {
  if(!addES & !decreasing) { return(factor(dd, levels = rev(c(1.5,1,0.5,0.25)))) }
  if(!addES & decreasing) { return(factor(dd, levels = c(1.5,1,0.5,0.25))) }
  if(addES & !decreasing) { return(factor(paste0("Effect size: ",dd), levels = rev(paste0("Effect size: ", c(1.5,1,0.5,0.25))))) }
  if(addES & decreasing) { return(factor(paste0("Effect size: ",dd), levels = paste0("Effect size: ", c(1.5,1,0.5,0.25)))) }
}

# Colours (hex values)
# 1. Six MAFs:
RedtoBlue6 <- c("#b2182b", "#d6604d", "#f4a582", "#4393c3", "#2166ac", "#053061")

# 2. Multiple testing correction methods:
# Colors of six hierarchical correction procedures
# "ST-BH","BH-BH","BY-BH","Bonferroni-BH","eigenMT-BH","BPerm1k-BH"
hier6color <- c("#d2478f","#fdb462","#8c64cf", "#62a85b","#1f78b4","#c95b4a")

# Colors of four global correction methods
# "ST", "BH", "BY", "Bonferroni"
step2m4 <- c("#b98d3e","#9970c1","#64a860","#cc545e")

# 3. Bootstrap estimators:
# Naive: red, weighted estimator: green, outof: purple, shrinkage: blue
estimators4 <- c("#ef6548", "#33a02c", "#9970ab", "#1f78b4")
# eGenes that were not identified: grey
color_notidentified <- "#bdbdbd"


#--------------------------------------------------------------------------------------------
# 1. Comparison of multiple testing correction methods
#--------------------------------------------------------------------------------------------
# Read FDR & TPR of all multiple tetsing correction methods
FDRTPRallmet <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/2_Multiple_testing_correction_method_comparison/FDR_TPR_all_methods_r2threshold_0.8_UnfixedEffectSize.txt", header = T)
colnames(FDRTPRallmet)[3] <- "Method"

# 35 methods; Perm1k and APerm10k only have 10 simulations
unique(FDRTPRallmet$Method)
unique(FDRTPRallmet$Method[which(FDRTPRallmet$n_sims == 10)])

#----------- A. Inflated FDR of pooled FDR methods -----------
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/MultipleTestingCorrection_method_comparison/Pooled_FDR_methods/")

# My function - change method name (will be shown as each plot's title)
changename <- function(mydata, whichcol, oldname, newname) {
  for(ii in 1:length(oldname)) {
    mydata[,whichcol][which(mydata[,whichcol] == oldname[ii])] <- newname[ii]
  }
  mydata[,whichcol] <- factor(mydata[,whichcol], levels = newname)
  return(mydata)
}

# Each plot is a pooled method and MAFs have different colors
FDRTPRallmet$MAFgroup <- MAFasgourps(FDRTPRallmet$MAF, addMAF = F, decreasing = T)

# Compare 4 pooled methods
metnames <- c("pooledST","pooledBH","pooledBY","pooledBonf")
# Keep data for these methods, and change names of the methods
pooledFDRTPR <- changename(FDRTPRallmet[which(FDRTPRallmet$Method %in% metnames),], whichcol = "Method", oldname = metnames, newname = c("Pooled ST", "Pooled BH", "Pooled BY", "Pooled Bonf"))

# Line plots
# 1.FDR
pooledFDR <- ggplot(pooledFDRTPR, aes(x = SS, y = FDR, colour = MAFgroup)) +
  geom_point() + geom_line() +
  geom_hline(yintercept = 0.05, linetype = 2) + 
  scale_x_log10(breaks = c(100,200,500,1000,2000,5000)) + 
  facet_wrap(~Method, nrow = 1) + xlab("Sample Size") + 
  scale_color_manual(values = RedtoBlue6, name = "MAF")
ggsave(filename = "FDR_pooled_methods.pdf", pooledFDR, width = 10, height = 2.7)
# 2.TPR
pooledTPR <- ggplot(pooledFDRTPR, aes(x = SS, y = TPR, colour = MAFgroup)) + 
  geom_line() + geom_point() +
  scale_x_log10(breaks = c(100,200,500,1000,2000,5000)) + 
  facet_wrap(~Method, nrow = 1) + xlab("Sample Size") + ylab("TPR") +
  scale_color_manual(values = RedtoBlue6, name = "MAF")
ggsave(filename = "TPR_pooled_methods.pdf", pooledTPR, width = 10, height = 2.7)


#----------- B. Comparison of hierarchical correction procedures -----------
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/MultipleTestingCorrection_method_comparison/Hierarchical_correction_procedures/")

# My function - line plots
LinePlot <- function(alldata, metnames, whatisylab, chosencolor = NULL, namesforlegend = NULL) {
  # Extract data used in plotting 
  dd <- alldata[which(alldata$Method %in% metnames),]
  dd$Method <- factor(dd$Method, levels = metnames)
  # Line plots
  if(whatisylab == "FDR") {
    p <- ggplot(dd, aes(SS, FDR, color = Method)) + geom_hline(yintercept = 0.05, linetype = 2)
  } else if(whatisylab == "TPR") {
    p <- ggplot(dd, aes(SS, TPR, color = Method)) + expand_limits(y = c(0, 1))
  }
  p <- p + geom_point() + 
    xlab("Sample size") + ylab(whatisylab) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.text.x = element_text( angle = 45, hjust = 1) ) +
    scale_x_log10(breaks = c(100,200,500,1000,2000,5000)) + 
    geom_line() +
    facet_wrap(~MAFgroup, nrow = 1, scales = "free_y")
  # If other colors are provided
  if(length(chosencolor) == length(metnames) & length(namesforlegend) == length(metnames)) {
    p <- p + scale_color_manual(values = chosencolor, labels = namesforlegend)
  } else if(length(chosencolor) == length(metnames)) {
    p <- p + scale_color_manual(values = chosencolor)
  } else {
    p <- p + scale_color_brewer(palette = "Paired")
  }
  return(p)
}

# Each plot has a different MAF
FDRTPRallmet$MAFgroup <- MAFasgourps(FDRTPRallmet$MAF, addMAF = T, decreasing = F)

#----- 1. Compare hierarchical methods (global method: BH) -----
metnames <- paste0(c("ST","BH","BY","Bonferroni","eigenMT","BPerm1k"), "_", "BH")
namesforlegend <- paste0(c("ST","BH","BY","Bonferroni","eigenMT","BPerm1k"), "-", "BH")

# FDR
# Maximun FDR when MAF=10%: 0.0722
maxFDRMAF10 <- max(FDRTPRallmet[which(FDRTPRallmet$Method %in% metnames & FDRTPRallmet$MAF==10),"FDR"])
pFDR <- LinePlot(FDRTPRallmet, metnames, whatisylab = "FDR", chosencolor = my6color, namesforlegend) + expand_limits(y = c(0,0.07))
ggsave(filename = "FDR_6_hierarchical_met_r2threshold08.pdf", pFDR, width = 11, height = 2.9)
# Sensitivity
pTPR <- LinePlot(FDRTPRallmet, metnames, whatisylab = "TPR", chosencolor = my6color, namesforlegend)
ggsave(filename = "TPR_6_hierarchical_met_r2threshold08.pdf", pTPR, width = 11, height = 2.9)

#----- 2. Compare four global correction methods -----
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/MultipleTestingCorrection_method_comparison/Hierarchical_correction_procedures_global_methods/")
# Four global correction methods
globalmet4 <- c("ST","BH","BY","Bonferroni")

# Local correction methods:
for(localmet in c("ST","BH","BY","Bonferroni","eigenMT","BPerm1k")) {
  metnames <- paste0(localmet, "_", globalmet4)
  namesforlegend <- paste0(localmet, "-", globalmet4)
  # FDR
  pFDR <- LinePlot(FDRTPRallmet, metnames, whatisylab = "FDR", chosencolor = step2m4, namesforlegend) + expand_limits(y = c(0,0.08))
  ggsave(filename = paste0("FDR_4globalmet_",localmet,"_r2threshold08.pdf"), pFDR, width = 11, height = 2.9)
  # Sensitivity
  pTPR <- LinePlot(FDRTPRallmet, metnames, whatisylab = "TPR", chosencolor = step2m4, namesforlegend)
  ggsave(filename = paste0("TPR_4globalmet_",localmet,"_r2threshold08.pdf"), pTPR, width = 11, height = 2.9)
}

#----- 3. Compare permutation approaches (10 simulations) -----
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/MultipleTestingCorrection_method_comparison/Hierarchical_correction_procedures_permutations_10sim/")

# Read BPerm1k FDR/TPR in first 10 simulations
FDRTPRBPerm1k10sim <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/2_Multiple_testing_correction_method_comparison/FDR_TPR_BPerm1k_10sim_r2threshold_0.8_UnfixedEffectSize.txt", header = T)
colnames(FDRTPRBPerm1k10sim)[3] <- "Method"
# Other permutation approaches in 10 simulations
otherperm <- FDRTPRallmet[which(FDRTPRallmet$n_sims == 10),]
# Merge into one tabel
allperm <- rbind(FDRTPRBPerm1k10sim, otherperm[,1:10])

# Each plot has a different MAF
allperm$MAFgroup <- MAFasgourps(allperm$MAF, addMAF = T, decreasing = F)

# Methods
metnames <- paste0(rep(c("Perm1k","APerm10k","BPerm1k"), each=2), "_", rep(c("ST","BH"),3))
namesforlegend <- paste0(rep(c("Perm1k","APerm10k","BPerm1k"), each=2), "-", rep(c("ST","BH"),3))

# Colors: global ST - light colours; global BH - dark colours
# Perm1k: green; APerm10k: red; BPerm1k: purple
perm6colors <- brewer.pal(12,"Paired")[c(1,2,5,6,9,10)]
# FDR
maxpermFDR <- max(allperm[which(allperm$Method %in% metnames & allperm$MAF>=10),"FDR"])     # Maximum: 0.0838
pFDR <- LinePlot(allperm, metnames, whatisylab = "FDR", chosencolor = perm6colors, namesforlegend) + expand_limits(y = c(0,0.084))
ggsave(filename = paste0("FDR_permutations_10sim_",localmet,"_r2threshold08.pdf"), pFDR, width = 11, height = 2.9)
# Sensitivity
pTPR <- LinePlot(allperm, metnames, whatisylab = "TPR", chosencolor = perm6colors, namesforlegend)
ggsave(filename = paste0("TPR_permutations_10sim_",localmet,"_r2threshold08.pdf"), pTPR, width = 11, height = 2.9)

#----- Permutations as local, compare four global methods (10 simulations) -----
for(localmet in c("Perm1k","APerm10k","BPerm1k")) {
  metnames <- paste0(localmet, "_", globalmet4)
  namesforlegend <- paste0(localmet, "-", globalmet4)
  # FDR
  pFDR <- LinePlot(allperm, metnames, whatisylab = "FDR", chosencolor = step2m4, namesforlegend) + expand_limits(y = c(0,0.084))
  ggsave(filename = paste0("FDR_4globalmet_permutations_10sim_",localmet,"_r2threshold08.pdf"), pFDR, width = 11, height = 2.9)
  # Sensitivity
  pTPR <- LinePlot(allperm, metnames, whatisylab = "TPR", chosencolor = step2m4, namesforlegend)
  ggsave(filename = paste0("TPR_4globalmet_permutations_10sim_",localmet,"_r2threshold08.pdf"), pTPR, width = 11, height = 2.9)
}


#--------------------------------------------------------------------------------------------
# 2. Power estimation using constant effect size
#--------------------------------------------------------------------------------------------
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/Power_estimation/")

#powerdata <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/3_Multiple_eSNPs_per_eGene_top_SNP_causal/ConsEffectSize/Proportion_top_SNP_causal_TPeGenes_LDstructure_eigenMT_BH_ConsEffectSize.txt", header = T)
powerdata <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/2_Multiple_testing_correction_method_comparison/FDR_TPR_all_methods_r2threshold_0.8_ConsEffectSize.txt", header = T)
# large effect size should be plotted first
powerdata <- powerdata[order(powerdata$effectsize, decreasing = T),]

# MAF as factor
powerdata$MAF <- MAFasgourps(powerdata$MAF, addMAF = T, decreasing = F)
# Effect size as factor
powerdata$effectsize <- ESasgourps(powerdata$effectsize, addES = F, decreasing = T)
# Local correction using eigenMT
powerdata <- powerdata[which(powerdata$Method == "eigenMT_BH"),]

# Plot power against SS
p <- ggplot(powerdata, aes(x=SS, y=TPR, color = effectsize)) + 
  geom_point() + facet_wrap(~MAF, nrow = 2) + 
  geom_hline(yintercept = 0.8, linetype = "dashed", size = 0.8) + 
  scale_color_brewer(palette = "PuOr") +
  xlab("Sample size") + ylab("Power") + labs(color = "Effect size") +
  scale_x_log10(breaks=c(100,200,500,1000,2000,5000)) + 
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  theme(axis.text.x = element_text( angle = 45, hjust = 1) ) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE, size = 0.8)
ggsave("Power_Plot_eigenMT-BH_binomial_smooth_r2threshold08.pdf", p, width = 6.5, height = 4)

# Large plot for flowchart
#p <- ggplot(data_power[which(data_power$MAF == "MAF: 10%"),], aes(x=SS,y=Bonferroni.BH, color = effectsize)) + 
#  geom_point() + geom_hline(yintercept = 0.8, linetype = "dashed", size = 1) + 
#  scale_color_brewer(palette = "PuOr") +
#  xlab("Sample size") + ylab("Power") + 
#  theme(legend.position = "none") + 
#  scale_x_log10(breaks=c(100,200,500,1000,2000,5000)) + scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
#  theme(axis.text = element_blank(), axis.ticks = element_blank() ) +#, axis.title = element_text(face = "bold") ) +
#  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE)
#ggsave("Large_Power_Plot_Bonferroni-BH_binomial_smooth.pdf", p, width = 1.2, height = 1.2)


#--------------------------------------------------------------------------------------------
# 3. Multiple eSNPs per TP eGene & identification of causal eSNP
#--------------------------------------------------------------------------------------------
#----- The number of eSNPs per TP eGene -----
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/Multiple_eSNPs_per_eGene/")

## Constant effect size ##
ave_num_eSNPs_cons <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/3_Multiple_eSNPs_per_eGene_top_SNP_causal/ConsEffectSize_rethreshold0.8/Average_number_eSNPs_per_TPgeneORsiggene_100reps_eigenMT_BH_ConsEffectSize.txt", header = T)

# MAF as factor
ave_num_eSNPs_cons$MAF <- MAFasgourps(ave_num_eSNPs_cons$MAF, addMAF = F, decreasing = T)
# Effect size as factor
ave_num_eSNPs_cons$effectsize <- ESasgourps(ave_num_eSNPs_cons$effectsize, addES = T, decreasing = F)

# eigenMT-BH
p <- ggplot(ave_num_eSNPs_cons, aes(x=SS, y=n_SNP_TPeGene_eigenMT_BH, color = MAF)) + 
  geom_point() + geom_line() + 
  facet_wrap(~effectsize, nrow = 1) +
  scale_color_manual(values = RedtoBlue6) + 
  scale_x_log10(breaks = c(100,200,500,1000,2000,5000)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Sample size") + ylab("Average number of eSNPs")
ggsave("Number_eSNPs_per_TP_eigenMTBH_ConsEffectSize.pdf", p, width = 8, height = 2.4)

## Unfixed effect size ##
#ave_num_eSNPs_unf <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/3_Multiple_eSNPs_per_eGene_top_SNP_causal/UnfixedEffectSize/Average_number_eSNPs_per_TPgeneORsiggene_100reps_eigenMT_BH_UnfixedEffectSize.txt", header = T)

# MAF as factor
#ave_num_eSNPs_unf$MAF <- MAFasgourps(ave_num_eSNPs_unf$MAF, addMAF = F, decreasing = T)

# eigenMT-BH
#p <- ggplot(ave_num_eSNPs_unf, aes(x=SS, y=n_SNP_TPeGene_eigenMT_BH, color = MAF)) + 
#  geom_point() + geom_line() + 
#  scale_color_manual(values = RedtoBlue6) + 
#  scale_x_log10(breaks = c(100,200,500,1000,2000,5000)) +
#  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#  xlab("Sample size") + ylab("Average number of eSNPs")
#ggsave("Number_eSNPs_per_TP_eigenMTBH_UnfixedEffectSize.pdf", p, width = 3.3, height = 2.4)


#----- Proportion of top SNPs that are causal -----
# Plot the proportion against (1) power or (2) LD structure (color by MAF)
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/Multiple_eSNPs_per_eGene/")

## Constant effect size ##
# Read the proportion and LD structure
proportion_causal <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/3_Multiple_eSNPs_per_eGene_top_SNP_causal/ConsEffectSize_rethreshold0.8/Proportion_top_SNP_causal_TPeGenes_LDstructure_eigenMT_BH_ConsEffectSize.txt", header = T)

# Calculate power
proportion_causal$power <- proportion_causal$nTP/20000

# MAF as factor
proportion_causal$MAF <- MAFasgourps(proportion_causal$MAF, addMAF = F, decreasing = T)

# 1. Plot proportion against power
p <- ggplot(proportion_causal, aes(power, proportion, color = MAF)) + 
  geom_line() + geom_point() + 
  scale_color_manual(values = RedtoBlue6) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Power (eigenMT-BH) \n ") + ylab("Proportion of top eSNPs that were causal") +
  theme(axis.title = element_text(size=10))
ggsave("Proportion_top_eSNPs_causal_TP_VS_power_alldata_ConsEffectSize.pdf", p, width = 4, height = 3.2)

# Exclude scenarios where power is < 1% (200 TP eGenes across all 100 simulations)
p <- ggplot(proportion_causal[which(proportion_causal$power>=0.01),], aes(power, proportion, color = MAF)) +
  geom_line() + geom_point() + 
  scale_color_manual(values = RedtoBlue6) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Power (eigenMT-BH) \n ") + ylab("Proportion of top eSNPs that were causal") +
  theme(axis.title = element_text(size=10))
ggsave("Proportion_top_eSNPs_causal_TP_VS_power_power001_MAF_colors_ConsEffectSize.pdf", p, width = 4, height = 3.2)

# 2. Plot proportion against LD structure
p <- ggplot(proportion_causal, aes(ave_num_SNP05_2Mb, proportion, color = MAF)) + 
  geom_point() +
  scale_color_manual(values = RedtoBlue6) +
  scale_y_continuous(labels = scales::percent) +
  xlab("LD structure\n(Average number of SNPs in LD with the causal eSNP)") + 
  ylab("  ") +
  theme(axis.title = element_text(size=10))
ggsave("Proportion_top_eSNPs_causal_TP_VS_LD_structure_alldata_ConsEffectSize.pdf", p, width = 4, height = 3.2)
# Exclude scenarios where power is < 1% (200 TP eGenes across all 100 simulations) 
p <- ggplot(proportion_causal[which(proportion_causal$power>=0.01),], aes(ave_num_SNP05_2Mb, proportion, color = MAF)) + 
  geom_point() +
  scale_color_manual(values = RedtoBlue6) +
  scale_y_continuous(labels = scales::percent) +
  xlab("LD structure\n(Average number of SNPs in LD with the causal eSNP)") + 
  ylab("  ") + 
  theme(axis.title = element_text(size=10)) + expand_limits(x = 0)
ggsave("Proportion_top_eSNPs_causal_TP_VS_LD_structure_power001_ConsEffectSize.pdf", p, width = 4, height = 3.2)


#----- LD r2 between top SNP (not causal) and causal eSNP -----
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/Multiple_eSNPs_per_eGene/")

# Read LD r2 of top SNPs (TP eGenes) with causal eSNPs
LDr2 <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/3_Multiple_eSNPs_per_eGene_top_SNP_causal/ConsEffectSize_rethreshold0.8/LDr2_list_top_SNP_causal_TPeGenes_eigenMT_BH_ConsEffectSize.txt", header = T)

# NA means r2 ≤ 0.01
LDr2$LDr2[which(is.na(LDr2$LDr2))] <- 0

# Top SNPs that are not causal
LDr2_notcausal <- LDr2[which(LDr2$LDr2 != 1),]

# 10 classes based on r2
LDr2_notcausal$class <- sapply(1:nrow(LDr2_notcausal), function(ii) {
  b1 <- floor(LDr2_notcausal$LDr2[ii] *10)
  b2 <- b1 + 1
  return(paste0(b1/10,"-",b2/10))
})

# Histogram
p <- ggplot(LDr2_notcausal, aes(class)) + 
  geom_bar(width = 0.98, aes(y = (..count..)/sum(..count..))) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text( angle = 45, hjust = 1) ) +
  xlab("Correlation between non-causal top SNP and causal SNP") + ylab("Frequency") +
  theme(axis.title.x = element_text(size=8))
ggsave("R2_distribution_topSNP_not_causal_ConsEffectSize.pdf", p, width = 3.5, height = 3.5)

# How many ≥ 0.8
sum(LDr2_notcausal$LDr2 >= 0.8)/nrow(LDr2_notcausal)    # 83.53%




#--------------------------------------------------------------------------------------------
# 4. Winner's Curse in effect size estimation & Bootstrap
#--------------------------------------------------------------------------------------------
# Multiple testing correction method: eigenMT-BH
# Bootstrap method: same SNP in detection group, nominal thresholds to determine significant eGenes, 200 bootstraps

#----- A. Examples of Winner's Curse (Unfixed Effect Size) -----
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/WinnersCurse_Bootstrap_est/Examples/")

# Read estimates from unfixed effect size
bootsest_unfixed <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/4_WinnersCurse_Bootstrapping/UnfixedEffectSize_r2threshold0.8/BootstrapQTL_estimators_Naive_not_identified_eigenMTBH_UnfixedEffectSize_sim10.txt", header = T, sep = "\t")

# Examples to show Winner's Curse and bootstrap correction
# Sample size 200, MAF 5%, 10%, and 25%
smallsce <- bootsest_unfixed[which(bootsest_unfixed$SS %in% c(200) & bootsest_unfixed$MAF %in% c(5,10,25) & bootsest_unfixed$Estimator %in% c("Naive","shrinkage")),]
smallsce$estimates <- abs(smallsce$estimates)
smallsce$Estimator <- factor(smallsce$Estimator, levels = c("Naive","shrinkage"))
smallsce$scenario <- factor(paste0("Size=",smallsce$SS," MAF=",smallsce$MAF,"%"), levels = paste0("Size=200 MAF=", c(5,10,25),"%"))
# Shuffle the rows
smallsce <- smallsce[sample(1:nrow(smallsce), nrow(smallsce), replace = F),]
# Plotting
p <- ggplot(smallsce, aes(sim_beta, estimates, color = Estimator)) + 
  geom_point(size = 0.5, alpha = 0.6) +
  facet_wrap(~scenario, nrow = 1) + xlim(c(0,2)) + ylim(c(0,2)) +
  scale_color_manual(values = estimators4[c(1,4)]) +
  xlab("True effect size") + ylab("Estimated effect size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  geom_smooth(method = "lm", formula = y ~ x, se=FALSE, fullrange = T) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, size = 1)
ggsave(p, filename = "Examples_low_power.pdf", width = 5, height = 2)


#----- B. Comparison of estimators (Constant Effect Size) ----- 
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/WinnersCurse_Bootstrap_est/Comparison_estimators/")

# Read MSE, Mean Ratio, and Median Error for constant effect sizes
bootsperformance <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/4_WinnersCurse_Bootstrapping/ConsEffectSize_r2threshold0.8/BootstrapQTL_estimators_Naive_MSE_MedianError_MeanRatio_eigenMTBH_ConsEffectSize_sim10.txt", header = T, sep = "\t")
# Remove scenarios where power was 0 and data is NA
bootsperformance <- bootsperformance[which(bootsperformance$power > 0),]
# Effect size as factor
bootsperformance$effectsizeclass <- factor(paste0("True effect size=",bootsperformance$effectsize))

# 1) MSE
MSE <- bootsperformance[,c("SS","MAF","effectsizeclass","power","MSE_nai","MSE_shri","MSE_oos","MSE_wei")]
MSEdd <- melt(MSE, measure.vars = colnames(MSE)[5:8])
# Estimator as factor
MSEdd$variable <- factor(MSEdd$variable, levels = c("MSE_nai","MSE_wei","MSE_oos","MSE_shri"))

# Plot MSE, Power ≥0.01
p <- ggplot(MSEdd[which(MSEdd$power>=0.01),], aes(power, value, color = variable)) + 
  geom_point(size = 1) + facet_wrap(~effectsizeclass, scales = "free", nrow = 1) +
  scale_color_manual(values = estimators4, labels = c("Naive estimator","Bootstrap weighted estimator","Bootstrap out-of-sample estimator","Bootstrap shrinkage estimator"), name = "Estimator") +
  ylab("Mean Squared Error") + xlab("Power") + theme(legend.position = "none") +
  geom_line()
ggsave("MSE_power001.pdf", p, width = 8, height = 2.3)
#ggsave("color_legend.pdf", p+theme(legend.position = "right"), width = 11, height = 2.3)

# 2) Mean Ratio 
MRatio <- bootsperformance[,c("SS","MAF","effectsizeclass","power","MR_nai","MR_shri","MR_oos","MR_wei")]
MRatiodd <- melt(MRatio, measure.vars = colnames(MRatio)[5:8])
# Estimator as factor
MRatiodd$variable <- factor(MRatiodd$variable, levels = c("MR_nai","MR_wei","MR_oos","MR_shri"))

# Plot Mean Ratio, Power ≥0.01
p <- ggplot(MRatiodd[which(MRatiodd$power>=0.01),], aes(power, value, color = variable)) + 
  geom_point(size = 1) + facet_wrap(~effectsizeclass, scales = "free", nrow = 1) +
  scale_color_manual(values = estimators4) + 
  expand_limits(y = range(MRatiodd$value[which(MRatiodd$power>=0.01)])) + #min and max
  ylab("Mean Ratio") + xlab("Power") + theme(legend.position = "none") +
  geom_hline(yintercept = 1, linetype = 2) + geom_line()
ggsave("MR_power001.pdf", p, width = 8, height = 2.3)

# If don't plot different effect sizes in differet plots
p <- ggplot(MRatiodd[which(MRatiodd$power>=0.01),], aes(power, value, color = variable)) + 
  geom_point(size = 1) + #facet_wrap(~effectsizeclass, scales = "free", nrow = 1) +
  scale_color_manual(values = estimators4) + 
  ylab("Mean Ratio") + xlab("Power") + theme(legend.position = "none") +
  geom_hline(yintercept = 1, linetype = 2) + geom_line()
ggsave("MR_power001_oneplot.pdf", p, width = 2.5, height = 2.3)

# 3) Median Error
ME <- bootsperformance[,c("SS","MAF","effectsizeclass","power","ME_nai","ME_shri","ME_oos","ME_wei","effectsize")]
MEdd <- melt(ME, measure.vars = colnames(ME)[5:8])
# Estimator as factor
MEdd$variable <- factor(MEdd$variable, levels = c("ME_nai","ME_wei","ME_oos","ME_shri"))

# Plot Median Error, Power ≥0.01
p <- ggplot(MEdd[which(MEdd$power>=0.01),], aes(power, value, color = variable)) + 
  geom_point(size = 1) + facet_wrap(~effectsizeclass, scales = "free", nrow = 1) +
  scale_color_manual(values = estimators4) + 
  ylab("Median Error") + xlab("Power") + theme(legend.position = "none") +
  geom_line()
ggsave("MedianError_power001.pdf", p, width = 8, height = 2.3)

# Median Error/true effect size
MEdd$scaledME <- MEdd$value/MEdd$effectsize
p <- ggplot(MEdd[which(MEdd$power>=0.01),], aes(power, scaledME, color = variable)) + 
  geom_point(size = 1) + geom_line() +
  facet_wrap(~effectsizeclass, nrow = 1) +#, scales = "free") +
  scale_color_manual(values = estimators4) + 
  #expand_limits(y = range(MEdd$scaledME[which(MEdd$power>=0.01)])) + #min and max
  ylab("Median Error/True effect size") + xlab("Power") + theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("MedianError_power001_scaled_by_effectsize.pdf", p, width = 7, height = 2.3)

# If don't plot different effect sizes in differet plots
p <- ggplot(MEdd[which(MEdd$power>=0.01),], aes(power, scaledME, color = variable)) + 
  geom_point(size = 1) + #facet_wrap(~effectsizeclass, scales = "free", nrow = 1) +
  scale_color_manual(values = estimators4) + 
  ylab("Median Error/True effect size") + xlab("Power") + theme(legend.position = "none") +
  geom_line()
ggsave("MedianError_power001_scaled_by_effectsize_oneplot.pdf", p, width = 2.5, height = 2.3)


#----- Supplementary Figures of all scenarios (Unfixed Effect Size) -----
# Point plot for all 36 scenarios
setwd("~/Data/Simulation_cis_eQTLs_2017_10_19/Plots_eigenMT_BH/WinnersCurse_Bootstrap_est/Examples/")

# Data of one simulation of 36 scenarios using unfixed effect sizes
# Three bootstrap estimators and naive estimator for TP eGenes;
# nominal betas for those that were not identified.
est_simall <- read.table("~/Data/Simulation_cis_eQTLs_2017_10_19/Simulation_Data_eigenMT_BH/4_WinnersCurse_Bootstrapping/UnfixedEffectSize_r2threshold0.8/BootstrapQTL_estimators_Naive_not_identified_eigenMTBH_UnfixedEffectSize_sim10.txt", header = T, sep = "\t")
est_sim1 <- est_simall[which(est_simall$Rep == 1),]

# 36 scenarios -> factor
est_sim1$scenario <- factor(paste0("Size=", est_sim1$SS," MAF=", est_sim1$MAF, "%"), levels = paste0("Size=", rep(c(5000,2000,1000,500,200,100),each=6), " MAF=", rep(c(0.5,1,5,10,25,50),6), "%"))

# 1). Naive estimator only
naive <- est_sim1[c(which(est_sim1$Estimator == "not identified"), which(est_sim1$Estimator == "Naive")),]
# Absolute value
naive$estimates <- abs(naive$estimates)

# Nominal betas of TP eGenes as well as those that were not identified
naive$Estimator <- factor(naive$Estimator, levels = c("Naive","not identified"))
p <- ggplot(naive, aes(sim_beta, estimates, color = Estimator)) + 
  geom_point(size = 0.3) + facet_wrap(~scenario, nrow = 6) +
  scale_color_manual(values = c(estimators4[1], color_notidentified), 
                     labels = c("Naive estimator", "Not identified")) +
  xlab("True effect size") + ylab("Estimated effect size") +
  theme(strip.text = element_text(size = 7),   # smaller plot labels
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 15),
        legend.title = element_blank(), # remove the legend title
        legend.position = "bottom") + # Put the legend at the bottom
  geom_smooth(data = subset(naive, Estimator=="Naive"),   # linear regression for only red dots
              method = "lm", formula = y ~ x, se = FALSE, size = 0.5, alpha = 0.6,
              fullrange = T) + # make the fit span full range of the plot, not just significant estimator
  ylim(c(0,2.6)) + # If we don't put limit on y axis, the liner regression would make y axis large
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.3)   # diagonal lines
# 1.Original y lim
ggsave("Pointplot_sim1_scenario36_naive_notidentified_originalylim.pdf", p + ylim(c(0,4)), width = 7, height = 8)
# 2. Smaller y lim
ggsave("Pointplot_sim1_scenario36_naive_notidentified.pdf", p, width = 7, height = 8)
# 3. Smaller y lim and no panel labels
ggsave("Pointplot_sim1_scenario36_naive_notidentified_nolabelfacet.pdf", 
       p+theme(strip.background = element_blank(), strip.text = element_blank()), 
       width = 7, height = 7)

# 2). Naive estimator vs bootstrap estimator
for(ee in c("shrinkage","outof","weighted")) {
  # Put grey/not detected dots in the back and shuffle naive and bootstrap estiator
  notdetected <- est_sim1[which(est_sim1$Estimator == "not identified"),]
  dd <- est_sim1[which(est_sim1$Estimator %in% c("Naive",ee)),]
  dd <- dd[sample(1:nrow(dd), nrow(dd)),]
  dd <- rbind(notdetected,dd)
  dd$Estimator <- factor(dd$Estimator, levels = c("Naive",ee,"not identified"))
  if(ee == "shrinkage") {eelab <- "Bootstrap shrinkage estimator";whatcolors <- c(estimators4[c(1,4)],color_notidentified)}
  if(ee == "outof") {eelab <- "Bootstrap out-of-sample estimator";whatcolors <- c(estimators4[c(1,3)],color_notidentified)}
  if(ee == "weighted") {eelab <- "Bootstrap weighted estimator";whatcolors <- c(estimators4[c(1,2)],color_notidentified)}
  
  pboot <- ggplot(dd, aes(sim_beta, abs(estimates), color = Estimator)) + 
    geom_point(size = 0.3) + facet_wrap(~scenario, nrow = 6) +
    scale_color_manual(values = whatcolors, labels = c("Naive estimator",eelab,"Not identified")) +
    xlab("True effect size") + ylab("Estimated effect size") +
    theme(strip.text = element_text(size = 7),   # smaller plot labels
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 15),
          legend.title = element_blank(), # remove the legend title
          legend.position = "bottom") + # Put the legend at the bottom
    geom_smooth(data = subset(dd, Estimator != "not identified"),   # linear regression for only red dots
                method = "lm", formula = y ~ x, se = FALSE, size = 0.5, alpha = 0.6,
                fullrange = T) + # make the fit span full range of the plot, not just significant estimator
    ylim(c(0,2.6)) + # If we don't put limit on y axis, the liner regression would make y axis large
    geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.3)   # diagonal lines
  
  ggsave(paste0("Pointplot_sim1_scenario36_naive_",ee,".pdf"), pboot, width = 7, height = 8)
}

# 3). Examples of scenarios with high power
# Examples to show that both corrected and nominal estimator estimate the effect size accurately in scenarios with sufficient power.
# Sample size 5k and 2k, MAF 25% and 50%
bigsce <- est_sim1[which(est_sim1$SS %in% c(2000,5000) & est_sim1$MAF %in% c(50,25) & est_sim1$Estimator %in% c("Naive","shrinkage")),]
bigsce$Estimator <- factor(bigsce$Estimator, levels = c("Naive","shrinkage"))
bigsce$scenario <- factor(bigsce$scenario, levels = paste0("Size=",c(2000,2000,5000,5000)," MAF=", c(25,50,25,50),"%"))
# Shuffle the rows
bigsce <- bigsce[sample(1:nrow(bigsce), nrow(bigsce), replace = F),]
# Plotting
p <- ggplot(bigsce, aes(sim_beta, abs(estimates), colour = Estimator)) + 
  geom_point(alpha = 0.8) +
  facet_wrap(~scenario, nrow = 1) + xlim(c(0,2)) + ylim(c(0,2)) +
  scale_color_manual(values = estimators4[c(1,4)], 
                     labels = c("Naive estimator", "Bootstrap shrinkage estimator")) +
  xlab("True effect size") + ylab("Estimated effect size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, size = 0.8, alpha = 0.8, fullrange = T) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, size = 0.8)
ggsave(p, filename = "Examples_high_power.pdf", width = 6, height = 2.6)




#--------------------------------------------------------------------------------------------
#----- The inflated FDR of pooled methods vs LD spectrum -----
# Plotting FDR of pooled methods in LD pruning datasets and no LD dataset
# Directory of Plots
setwd("~/Data/Simulation_Replicates_2017.03.27_hg19_GRCh37/Plots/Fig2_Pooled_methods_infalted_FDR/LDpruning_noLD/")

# LD pruning 0.3, LDp pruning 0.1, no LD
for(whichLD in c("LDpruned03","LDpruned01","noLD")) {
  # Read TPR and FDR
  TPR <- read.table(paste0("~/Data/Simulation_Replicates_2017.03.27_hg19_GRCh37/Simulation_Data/MTC_methods_10reps/LDpruning_noLD/Rep10_whole_TPR_36scenarios_",whichLD,".txt"), header = T, stringsAsFactors = F)
  FDR <- read.table(paste0("~/Data/Simulation_Replicates_2017.03.27_hg19_GRCh37/Simulation_Data/MTC_methods_10reps/LDpruning_noLD/Rep10_whole_FDR_36scenarios_",whichLD,".txt"), header = T, stringsAsFactors = F)
  # Plot large SS and MAF first
  TPR <- TPR[rev(order(TPR$SS, TPR$MAF)),]
  FDR <- FDR[rev(order(FDR$SS, FDR$MAF)),]
  # Add "%" to MAF values and make it as factor
  TPR$MAF <- MAFtoFactor(TPR$MAF, F)
  FDR$MAF <- MAFtoFactor(FDR$MAF, F)
  # Methods
  name.method <- c("pooledST","pooledBH","pooledBY","pooledBonf")
  
  # FDR
  FDR.pooled <- FDR[, c("SS", "MAF", name.method)]
  names(FDR.pooled)[3:6] <- c("Pooled ST", "Pooled BH", "Pooled BY", "Pooled Bonferroni")
  FDR.pooled <- melt(FDR.pooled, id.vars = c("SS","MAF"))
  names(FDR.pooled)[3:4] <- c("Method", "FDR")
  # TPR
  TPR.pooled <- TPR[,c("SS", "MAF", name.method)]
  names(TPR.pooled)[3:6] <- c("Pooled ST", "Pooled BH", "Pooled BY", "Pooled Bonferroni")
  TPR.pooled <- melt(TPR.pooled, id.vars = c("SS","MAF"))
  names(TPR.pooled)[3:4] <- c("Method", "TPR")
  
  # FDR
  p.pooled <- ggplot(FDR.pooled, aes(x = SS, y = FDR, colour = MAF)) + geom_line() + geom_point() +
    geom_hline(yintercept = 0.05, linetype = 2) + 
    scale_x_log10(breaks = c(100,200,500,1000,2000,5000)) + 
    facet_wrap(~Method, nrow = 1) + xlab("Sample Size") + ylab("FDR") +
    scale_color_manual(values = my6color)
  ggsave(filename = paste0("FDR_pooled_FDRmethods_pooledBonf_lineplots",whichLD,".pdf"), p.pooled, width = 10, height = 2.7)
  # TPR
  p.pooled <- ggplot(TPR.pooled, aes(x = SS, y = TPR, colour = MAF)) + geom_line() + geom_point() +
    scale_x_log10(breaks = c(100,200,500,1000,2000,5000)) + 
    facet_wrap(~Method, nrow = 1) + xlab("Sample Size") + ylab("TPR") +
    scale_color_manual(values = my6color)
  ggsave(filename = paste0("TPR_pooled_FDRmethods_pooledBonf_lineplots_",whichLD,".pdf"), p.pooled, width = 10, height = 2.7)
}










