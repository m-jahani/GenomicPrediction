#Mojtaba Jahani
#November 2021
#######################################################Load Packages########################################################

library(rrBLUP)
library(data.table)
library(tidyverse)
library(modelr)
#########################################################Input data#########################################################

args = commandArgs(trailingOnly = TRUE)
MARKER <- args[1] #genotype data set "/data/users/mjahani/JOON_PAV/pav_gwas/SAM.5maf.noheader.sed12.rrblup.csv"
PHENOTYPE <- args[2] #phenotypes"/data/users/mjahani/JOON_PAV/pav_gwas/phenotype.csv"
trait <-  args[3] #column of the trait (1:53)
SAVE_DIR <- args[4] #directory to save result "/data/users/mjahani/JOON_PAV/pav_gwas/result"
fold <- args[5]

#MARKER <-"/data/users/mjahani/JOON_PAV/pav_gwas/SAM.5maf.noheader.sed12.rrblup.csv"
#PHENOTYPE <-"/data/users/mjahani/JOON_PAV/pav_gwas/phenotype.csv"
#trait <- 1
#SAVE_DIR <- "/data/users/mjahani/JOON_PAV/pav_gwas/result"

##########################################################Read data#########################################################

Markers_impute <- fread(MARKER,header = F)#read genotype file
phenotype <- read.csv(PHENOTYPE,
                      header = T,
                      row.names = 1) %>%
  scale(center = T,scale = T) %>%
  as.data.frame() %>%
  select(as.numeric(trait))#read phenotypic data

##########################################################Options###########################################################
#number of iterations for cross
iteration=20
##########################################################Test the accuracy########################################################
for (r in 1:iteration) {
  folds <- crossv_kfold(phenotype, k = as.numeric(fold)) 
  for (k in 1:as.numeric(fold)) {
    
    #define a subset as training population
    pheno_train <- phenotype[as.integer(folds$train[[k]]),]
    m_train <- Markers_impute[as.integer(folds$train[[k]]),]

    #define a subset as validation population
    pheno_valid <- phenotype[as.integer(folds$test[[k]]),]
    m_valid <- Markers_impute[as.integer(folds$test[[k]]),]
    
  answer_train <- mixed.solve(pheno_train, Z=m_train, K=NULL, method="REML",SE=FALSE, return.Hinv=FALSE)
  BLUE <- as.vector(answer_train$beta)
  GEBV_valid <- (as.matrix(m_valid) %*%  (as.matrix(answer_train$u))) + BLUE
 
  data.frame(phenotype = colnames(phenotype),#trait name
             Kfold = as.numeric(fold),
             it = iteration,
             Nfold = k,
             corr = cor(GEBV_valid,pheno_valid,use="complete"),
             Rsq = summary(lm(pheno_valid~GEBV_valid))$r.squared) %>% #calculate and save correlation between all GEBV and phenotype in validation population with all iterations
    fwrite(paste0(SAVE_DIR,
                  "/Kfold_cross_validation_cor_rsq_data"),
           append = T,
           sep = "\t",
           col.names = F)
  
  rm(pheno_train,
     m_train,
     pheno_valid,
     m_valid,
     answer_train,
     GEBV_valid,
     BLUE)
  }#folds
  rm(folds)
}#iteration

