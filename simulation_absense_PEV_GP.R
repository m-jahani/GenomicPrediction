#train a rrblup genomic prediction model with PAVs and predict the phenotype for a range of absence allele from 0 to 100
#Mojtaba Jahani
#November 2021
#######################################################Load Packages########################################################

library(rrBLUP)
library(data.table)
library(tidyverse)

#########################################################Input data#########################################################

args = commandArgs(trailingOnly = TRUE)
MARKER <- args[1] #genotype data set "/data/users/mjahani/JOON_PAV/pav_gwas/SAM.5maf.noheader.sed12.rrblup.csv"
PHENOTYPE <- args[2] #phenotypes"/data/users/mjahani/JOON_PAV/pav_gwas/phenotype.csv"
trait <-  args[3] #column of the trait (1:53)
SAVE_DIR <- args[4] #directory to save result "/data/users/mjahani/JOON_PAV/pav_gwas/result"

##########################################################Read data#########################################################

Markers_impute <- fread(MARKER,header = F)#read genotype file
phenotype <- read.csv(PHENOTYPE,
                      header = T,
                      row.names = 1) %>%
  scale(center = T,scale = T) %>%
  as.data.frame()#read scale phenotypic data
##################################################Fit the prediction Model##################################################

GP_MODEL <- mixed.solve(phenotype[,as.numeric(trait)], #phenotyoic data
                        Z=Markers_impute, #marker data
                        K=NULL, 
                        method="REML",
                        SE=FALSE, 
                        return.Hinv=FALSE)#calculation of the model
BLUPS <- as.matrix(GP_MODEL$u) #extract matix of BLUPS
BLUE <- as.vector(GP_MODEL$beta) #extract the BLUE value
HER <- as.vector(GP_MODEL$Vu)/(as.vector(GP_MODEL$Vu)+(as.vector(GP_MODEL$Ve)))
GEBVS <- (as.matrix(Markers_impute) %*% (BLUPS))+BLUE#calculate breeding values
NMARKER <- ncol(Markers_impute)  
rm(Markers_impute,GP_MODEL)
fit=lm(phenotype[,as.numeric(trait)]~GEBVS)
############################################explained variance prediction model accuracy##############################################

data.frame(phenotype = colnames(phenotype)[as.numeric(trait)],#trait name
           rsq = summary(fit)$r.squared,
           beta1 = summary(fit)$coefficients[2,1],
           pbeta1 = summary(fit)$coefficients[2,4]) %>% #save BLUE value
  fwrite(paste0(SAVE_DIR,"/NONCDS_explained_variance"),
           append = T,
           sep = "\t",
           col.names = F)


rm(GEBVS)

#############################################absence allele range simulation###############################################

  for (absence in seq(0,100,0.5)) { #loop through precent absent allele simulation
    matrix(1, 
           100, 
           NMARKER) -> X #matrix with 100 row and NMARKER columns with all elemnts set to 1 (presence)
    for(DUMMY in 1:nrow(X)) { #loop through each row of matrix to assign absence alleles
      X[DUMMY,
        sample(1:NMARKER,
               floor(absence*0.01*NMARKER))] <- -1 
    }#matrix of dummy markers
    ((X  %*% (BLUPS)) + BLUE) %>%
      as.data.frame() %>% 
      mutate(absense_percentage=absence,
             phenotype = colnames(phenotype)[as.numeric(trait)]) %>%
      fwrite(paste0(SAVE_DIR,"/NONCDS_absense_range_simulation_GP"),
             append = T,
             sep = "\t",
             col.names = F)
    rm(X)
  }#percentage range loop
