############################################
############# Basic parameters #############
# Which Reps:
whichreps <- 1:10
howmanyreps <- length(whichreps)
typeeffectsize <- "Cons"
localcormet <- "eigenMT"
globalcormet <- "BH"
############################################
# 2017-10-20 [FINISHED]
# Aim: to compare estimators using MSE, MR, Median Error, ect

library(foreach)
library(doMC); registerDoMC(cores = 20)
library(reshape)
options(stringsAsFactors = F)

# Version of BootstrapQTL:
versionbootsQ <- "0.7.3"
n_bootstraps <- 200

# Data are under:
mydir <- "~/Simulation_eQTL_replicates/scenario_files/Bootstrap_reestimation_effect_size/BootstrapQTL_v0.7.3_eigenMT-BH_nomthres_200_ConsEffectSize/"
setwd(mydir)

# Read all estimators from 10 simulations
bootsest <- read.table( paste0("BootstrapQTL_estimators_Naive_not_identified_",localcormet,globalcormet,"_",typeeffectsize,"EffectSize_sim",howmanyreps,".txt"), header = T, sep = "\t")

# Calculate Mean Squared Error, Mean Ratio, and Median Error
sce <- unique(bootsest[,4:6])
rownames(sce) <- 1:nrow(sce)
bootsperformances <- foreach(ii = 1:nrow(sce), .combine = rbind) %dopar% {
  dd <- bootsest[which(bootsest$SS == sce$SS[ii] & bootsest$MAF == sce$MAF[ii] & bootsest$effectsize == sce$effectsize[ii]),]
  dd$estimates <- abs(dd$estimates)
  return(data.frame(power = sum(dd$Estimator == "Naive")/2000,
                    MSE_nai = mean((dd[which(dd$Estimator == "Naive"),"estimates"]-dd$effectsize)^2),
                    MSE_shri = mean((dd[which(dd$Estimator == "shrinkage"),"estimates"]-dd$effectsize)^2),
                    MSE_oos = mean((dd[which(dd$Estimator == "outof"),"estimates"]-dd$effectsize)^2),
                    MSE_wei = mean((dd[which(dd$Estimator == "weighted"),"estimates"]-dd$effectsize)^2),
                    ME_nai = median(dd[which(dd$Estimator == "Naive"),"estimates"]-dd$effectsize),
                    ME_shri = median(dd[which(dd$Estimator == "shrinkage"),"estimates"]-dd$effectsize),
                    ME_oos = median(dd[which(dd$Estimator == "outof"),"estimates"]-dd$effectsize),
                    ME_wei = median(dd[which(dd$Estimator == "weighted"),"estimates"]-dd$effectsize),
                    MR_nai = mean(dd[which(dd$Estimator == "Naive"),"estimates"]/dd$effectsize),
                    MR_shri = mean(dd[which(dd$Estimator == "shrinkage"),"estimates"]/dd$effectsize),
                    MR_oos = mean(dd[which(dd$Estimator == "outof"),"estimates"]/dd$effectsize),
                    MR_wei = mean(dd[which(dd$Estimator == "weighted"),"estimates"]/dd$effectsize) ) )
}
bootsperformances <- cbind(sce, bootsperformances)
# Save the results
write.table(bootsperformances, paste0("BootstrapQTL_estimators_Naive_MSE_MedianError_MeanRatio_",localcormet,globalcormet,"_",typeeffectsize,"EffectSize_sim",howmanyreps,".txt"), quote = F, sep = "\t", row.names = F)



