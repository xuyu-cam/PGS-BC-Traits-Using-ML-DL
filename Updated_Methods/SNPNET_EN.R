library("snpnet")
library(stringr)
library(dplyr)
library(tidyr)
library("data.table")

args = commandArgs(trailingOnly=TRUE)


# trait id
trait = args[1]

# the folder hosts a expanded variant set (.pgen format)
vars_sele_folder =  args[2]

# one of the five sample partitions in UKB 
pheno_folder =  args[3]

#folder saves intermediate results
results_dir = paste(".../temp_results/",vars_sele_folder,"_",trait,pheno_folder,sep='')

# configure model settings
configs <- list(
  results.dir = results_dir, # needed when saving intermediate results
  save = TRUE, # save intermediate results per iteration (default FALSE)
  nCores=20, # number of cores available (default 1)
  alpha=0.5,
  mem = 290000,
  use.glmnetPlus = TRUE, # recommended for faster computation
  plink2.path = ".../plink_2.0/plink2" # path to plink2 program
)

# genotype file - including all chromosomes
genotype.pfile=paste(".../",vars_sele_folder,'/',trait,'_ukb_imp_v3_SNPs_maf0.01_ldthin0.5_full',sep='')

# phenotype file 
phenotype.file=paste(".../",trait,"_trait_levels_with_sample_split_",pheno_folder,'.txt',sep='')

# trait
phenotype <- trait

# covariates - [all the input phenotypes were already adjusted for covariates]
covariates <- c()


# run EN to find the best lamada
fit_snpnet <- snpnet(
  genotype.pfile = genotype.pfile,
  phenotype.file = phenotype.file,
  phenotype = phenotype,
  covariates = covariates,
  split.col = "split", # split column name in phenotype.file with train/val/test labels (pre-processed)
  configs = configs
)


# index of the best lamada (identified from the validaion data)
max_idx=which.max(fit_snpnet$metric.val)


#predict trait values of all models on testing samples (interval testing)
pred_snpnet <- predict_snpnet(
  fit = fit_snpnet,
  new_genotype_file = genotype.pfile, 
  new_phenotype_file = phenotype.file, 
  phenotype = phenotype, 
  covariate_names = covariates,   
  split_col = "split",
  split_name = c("test")
  )


# get the predicted PGSs of the selected optimal lamada for testings samples 
df <-melt(pred_snpnet$prediction$test)
col_name= paste('s',as.character(max_idx-1),sep='')
colnames(df)= c("IIID", "lamada", "value") 
df=df[df$lamada==col_name,]

# get te actual trait level data of all testing samples for a trait
df_pheno=fread(phenotype.file)
df_pheno$IIID = paste(df_pheno$FID,df_pheno$IID,sep='_')
df_pheno=df_pheno[df_pheno$split=='test']
df_merge = merge(x = df_pheno, y = df, by = "IIID", all = TRUE)

# get test metrics 
r_score = cor.test(df_merge[[trait]], df_merge$value, method = "pearson")
r_value = unname(r_score$estimate)
conf_int_low = r_score$conf.int[1]
conf_int_high = r_score$conf.int[2]
pval = r_score$p.value
r2_value = r_value*r_value

# for writing
test_r2 = r2_value
test_r = r_value


# create PGS model file with the trained model
df_beta=stack(fit_snpnet$beta[[max_idx]])
df_beta$rsid= df_beta$ind
df_beta=df_beta %>%
  separate(ind, c("chr", "pos","col3","other_allele","effect_allele"), "_")
colnames(df_beta)[1]='effect'
setcolorder(df_beta, c("rsid","chr","pos", "effect_allele","other_allele",'effect'))
df_beta$col3=NULL
df_beta=df_beta[order(df_beta$chr),]
df_beta$rsid= paste(df_beta$chr,df_beta$pos,df_beta$effect_allele,df_beta$other_allele,sep='_')



# create results dataframe for internal validation results
df_results <- data.frame(trait= character(0), N_vars= numeric(0),folder=numeric(0),r2=numeric(0),r=numeric(0),pval=numeric(0),confinf_low=numeric(0),confinf_high=numeric(0),stringsAsFactors=FALSE)
df_results[nrow(df_results)+1, ] <- c(trait, dim(df_beta)[1],pheno_folder,test_r2,test_r,pval,conf_int_low,conf_int_high)


# write PGS model - beta list
save_file_folder_beta=paste(".../PGS_models/",vars_sele_folder,sep='')
dir.create(save_file_folder_beta, showWarnings = FALSE)
save_file_beta = paste(save_file_folder_beta,'/',trait,'_PGS_model_',pheno_folder,sep='')
write.table(df_beta, file=save_file_beta,col.names=TRUE,row.names=FALSE,quote = FALSE,sep='\t')

# write internal validation results
save_file_folder_results=paste(".../results/",vars_sele_folder,sep='')
dir.create(save_file_folder_results, showWarnings = FALSE)
save_file_results = paste(save_file_folder_results,'/',trait,'_results_',pheno_folder,sep='')
write.table(df_results, file=save_file_results,col.names=TRUE,row.names=FALSE,quote = FALSE,sep='\t')



