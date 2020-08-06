#!/usr/bin/env Rscript

setwd('/notebook/Projects/clustering/')

source('clust_lca_aux.R')
source('clust_mclust_nbclust_aux.R')
source('clust_info_table.R')

setwd('/notebook/Desktop/catch me - new biomarkers/')
#setwd('/rds/homes/v/vxr610/gkoutosg-variant-prediction/BBCAF - new biomarkers')

#####################
# load packages

library(readr)
#library(mclust)
library(poLCA)
library(ggplot2)
library(dplyr)
library(stats)
library(dataPreparation)
library(xlsx)

# xlsx was tricky to install
# in the docker image
# apt-get install r-cran-rjava
# apt-get install openjdk-8-jdk
# R CMD javareconf -e
# R: install.packages('xlsx')
# restart R after installing


#####################
#load the data and print summary information
df <- read_csv('/notebook/Desktop/catch me - new biomarkers/BBC_AF_data_11Nov19.csv')
df <- df[!is.na(df$Filter),]

BIO_vars <- c('IL_6', 'hsCRP', 'ANG2', 'BMP10', 'ESM1', 'IGFBP7', 'NTproBNP', 'hsTnT', 'GDF15', 'FABP3', 'CA125', 'FGF23')
CONT_vars <- c('eGFR_CKDEPI', 'Age', 'BMI', BIO_vars)
CAT_vars <- c('Sex', 'Stroke_TIA', 'HF', 'HTN', 'AF', 'BMI_cat')
EXTRA_CAT_vars <- c(CAT_vars, 'n_eGFR') # <------------------------- these will be used
ALL_CAT <- c(EXTRA_CAT_vars, 'PAD', 'CAD', 'Diabetes', 'MOCA_GE_26')
COMORB_LIST <- c('Stroke_TIA', 'HF', 'HTN', 'AF', 'PAD', 'CAD', 'Diabetes')

colnames(df)[colnames(df) == 'MOCA≥26'] <- 'MOCA_GE_26'

df$n_eGFR <- as.factor(sapply(df$eGFR_CKDEPI, function(x) {
  if (x >= 90) {
    return('G1')
  } else if (x >= 60) {
    return('G2')
  } else if (x >= 45) {
    return('G3a')
  } else if (x >= 30) {
    return('G3b')
  } else if (x >= 15) {
    return('G4')
  } else {
    return('G5')
  }
}))

df$BMI_cat <- as.factor(sapply(df$BMI, function(x) {
  if (x < 18.5) {
    return('underweight')
  } else if (x < 25) {
    return('normal')
  } else if (x < 30) {
    return('overweight')
  } else if (x < 35) {
    return('obese')
  } else if (x < 40) {
    return('severe obese')
  } else {
    return('very severe obese')
  }
}))
df$AF <- df$Outcome_AF_SR
df$Outcome_AF_SR <- NULL

# set the categorical variables as factor: please check if they start with 0 just to avoid issues. 
# it is also important to check how the bmi categorical behaves; the order might be odd and should be verified in the end
fix_categorical <- function(df) {
  for(i in c(ALL_CAT)) { 
    if(i %in% colnames(df)) {
      df[[i]] <- as.factor(as.character(df[[i]]))
      print(paste(i, paste(levels(df[[i]]), collapse=' ')))
    }
  }
  return(df)
}
df <- fix_categorical(df)

bio_scales <- build_scales(dataSet = df[df$Filter == 1, ], cols = BIO_vars, verbose = T)
df <- data.frame(fastScale(dataSet = df, scales = bio_scales, verbose=T))

df_pca <- prcomp(df[, BIO_vars])
df <- cbind(df, df_pca$x)
PCA_vars <- colnames(df_pca$x)

print(summary(df))
#####################


#####################
# load parameters
if(F) {
  args = commandArgs(trailingOnly=TRUE)
  if (length(args) < 2) {
    stop('requires seed and number of groups')
  }
  
  current_seed <- args[1]
  groups <- args[2]
} else {
  current_seed <- 1
  groups <- 1:8
}
#####################

stop('The block ahead does LCA on the categorical variables')

#####################
# lca clustering

variables_used <- EXTRA_CAT_vars
formula <- create_formula(variables_used)
FILE_FORMAT <- '20191127_bootstrapping/20191127_LCA_bootstrapping_AF_all_patients_'
loop_operation_LCA(1:100, 1:8, df[df$Filter == 1, variables_used], formula, FILE_FORMAT, nrep=5, redo=F)
results <- load_results_LCA(1:100, 1:8, FILE_FORMAT)


variables_used <- setdiff(EXTRA_CAT_vars, c('AF'))
formula <- create_formula(variables_used)
FILE_FORMAT <- '20191127_bootstrapping/20191127_LCA_bootstrapping_noAF_all_patients_'
loop_operation_LCA(1:100, 1:8, df[df$Filter == 1, variables_used], formula, FILE_FORMAT, nrep=5, redo=F)
results <- load_results_LCA(1:100, 1:8, FILE_FORMAT)


variables_used <- setdiff(EXTRA_CAT_vars, c('AF'))
formula <- create_formula(variables_used)
FILE_FORMAT <- '20191127_bootstrapping/20191127_LCA_bootstrapping_noAF_noAF_patients_'
loop_operation_LCA(1:100, 1:8, df[df$AF == "0" & df$Filter == 1, variables_used], formula, FILE_FORMAT, nrep=5, redo=F)
results <- load_results_LCA(1:100, 1:8, FILE_FORMAT)


variables_used <- setdiff(EXTRA_CAT_vars, c('AF'))
formula <- create_formula(variables_used)
FILE_FORMAT <- '20191127_bootstrapping/20191127_LCA_bootstrapping_noAF_AF_patients_'
loop_operation_LCA(1:100, 1:8, df[df$AF == "1" & df$Filter == 1, variables_used], formula, FILE_FORMAT, nrep=1, redo=F)
results <- load_results_LCA(1:100, 1:8, FILE_FORMAT)


variables_used <- c(setdiff(EXTRA_CAT_vars, c('n_eGFR', 'AF')), 'Diabetes')
formula <- create_formula(variables_used)
FILE_FORMAT <- '20191127_bootstrapping/20191127_LCA_bootstrapping_only_AF_no_eGFR_w_Diabetes'
loop_operation_LCA(1:100, 1:8, df[df$AF == "1" & df$Filter == 1, variables_used], formula, FILE_FORMAT, nrep=1, redo=F)
results <- load_results_LCA(1:100, 1:8, FILE_FORMAT)

#####################
# compile the results for LCA

# boxplot with each class
ggplot(results, aes(x=nclasses, group=nclasses, y=bic)) + geom_boxplot()

t_test_results <- t.test(results[results$nclasses == 3,]$bic, results[results$nclasses == 4,]$bic)
t_test_results

# lines with the different bics
ggplot(results, aes(x=nclasses, group=seed, y=bic, color=seed)) + geom_line() + theme(legend.position = "none")

# frequency of the values
results %>% group_by(seed) %>% slice(which.min(bic)) -> a

frequency <- data.frame(table(a$nclasses))
colnames(frequency) <- c('nclasses', 'frequency')

ggplot(frequency, aes(x=nclasses, y=frequency)) + geom_bar(stat="identity")


#####################
# compute new model and get the distribution of the values
stop('need the best cluster option')
BEST_N <- 4

model <- run_single_LCA(df, formula, BEST_N, seed=1, nrep=1)

predicted_df <- df
predicted_df$predclass <- model$predclass

cat_table <- table_cat_values(predicted_df, c(ALL_CAT)) # <------------------ CHECK
comorb <- table_n_comorb(predicted_df, c(COMORB_LIST)) # <------------------ CHECK
cont_table <- table_continuous_values(predicted_df, setdiff(colnames(predicted_df), c(ALL_CAT, 'predclass', 'ID', 'predicted_df', "Pred_Prob", "BS_Prob", "Filter", "AgeSexBMI", "AgeSexBMI_NTproBNP", "AgeSexBMI_BMP10", "AgeSexBMI_BMP10ANG2", "MOCA_GE_26")))

#####################


stop('The block ahead does mclust/nbclust on the Biomarkers')

#####################
# mclust clustering (biomarkers)

FILE_FORMAT <- '20191127_bootstrapping/20191127_bootstrapping_all_patients_mclust_'
loop_operation_MCLUST_NBCLUST(1:100, 1:30, df[df$Filter==1, PCA_vars], FILE_FORMAT=FILE_FORMAT, redo=F, do_mclust=T, do_nbclust=F)
results_mclust <- load_results_MCLUST_NBCLUST(1:100, 1:30, FILE_FORMAT, do_mclust=T, do_nbclust=F)

results_mclust$groups <- as.numeric(as.character(results_mclust$groups))

mdf <- melt(results_mclust, measure.vars=colnames(results_mclust)[1:14], na.rm=T)
ggplot(mdf, aes(x=groups, group=groups, y=value)) + geom_boxplot()

mdf$group_variable <- paste(mdf$groups, mdf$variable)
ggplot(mdf, aes(x=groups, group=group_variable, fill=variable, y=value)) + geom_boxplot()

ggplot(mdf, aes(x=groups, group=variable, color=variable, y=value)) + geom_point() # VVI is the highest
ggplot(mdf, aes(x=groups, group=variable, color=variable, y=value)) + geom_smooth() # this is an approximation, but can be seem that VVI is higher

ggplot(mdf[mdf$variable == 'VVI',], aes(x=groups, group=groups, y=value)) + geom_boxplot()

##########
# nbclust clustering (biomarkers)

FILE_FORMAT <- '20191127_bootstrapping/20191127_bootstrapping_all_patients_nbclust_PCA_'
loop_operation_MCLUST_NBCLUST(1:100, 1:30, df[df$Filter==1, PCA_vars], FILE_FORMAT=FILE_FORMAT, redo=F, do_mclust=F, do_nbclust=T)
results_nbclust <- load_results_MCLUST_NBCLUST(1:100, 1:30, FILE_FORMAT, do_mclust=F, do_nbclust=T)

ggplot(results_nbclust, aes(x=Number_clusters, group=Number_clusters)) + geom_bar()
ggplot(results_nbclust, aes(y=Number_clusters, x=Number_clusters, group=seed)) + geom_boxplot()

results_nbclust %>% group_by(seed, Number_clusters) %>% count() -> grouped
grouped$percent <- grouped$n/26
ggplot(grouped, aes(x=Number_clusters, y=percent, group=Number_clusters)) + geom_boxplot()
ggplot(grouped, aes(x=Number_clusters, y=n, group=Number_clusters)) + geom_boxplot()

##########
# compiling results

set.seed(1)
c12 <- run_mclust(df[df$Filter == 1, BIO_vars],
                  groups=12,
                  modelNames=c('VVI'))
c12_df <- df[df$Filter == 1,]
c12_df$predclass <- c12$classification


set.seed(1)
c17 <- run_mclust(df[df$Filter == 1, BIO_vars],
                  groups=17,
                  modelNames=c('VVI'))
c17_df <- df[df$Filter == 1,]
c17_df$predclass <- c17$classification


n2 <- kmeans(df[df$Filter == 1, BIO_vars],
             centers=2,
             algorithm='Hartigan-Wong',
             nstart=1,
             iter.max = 100)
n2_df <- df[df$Filter == 1,]
n2_df$predclass <- n2$cluster
print(table(n2_df$predclass))

##########
# 

dataset <- n2_df
output_file <- '20191128_bbcaf_kmeans_2_df.xlsx'

comparison_df <- data.frame(fastScale(dataSet = dataset, scales = bio_scales, verbose=T, way='unscale'))
r <- compile_results_to_xlsx(comparison_df, 
                        continuous_variables = CONT_vars, 
                        categorical_variables = ALL_CAT,
                        comorbidity_variables = COMORB_LIST,
                        subgroup_cases = c(1, 2, 3),
                        output_file = output_file)

