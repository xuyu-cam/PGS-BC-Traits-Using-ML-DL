library(bigsnpr)


args <- commandArgs(trailingOnly = TRUE)

# the folder hosts data for a expanded variant set
vars_sele_folder = as.character(args[1])

# trait id
trait = as.character(args[2])



print(paste0("Start working with variant set:", vars_sele_folder,' , Trait:',trait,' -- ' ,Sys.time()))


# read summary stats data for a trait given a variant folder
sumstat_file = paste(".../genotype/ukb/ldthined_r2_0.5/",vars_sele_folder,'_gwas/',trait,"_gwas_normalised_imputed_r2_0.5_vars.out",sep='')
sumstats <- bigreadr::fread2(sumstat_file,nThread = 10)   
sumstats$tmp=sumstats$a0
sumstats$a0=sumstats$a1
sumstats$a1=sumstats$tmp
sumstats['tmp']<-NULL


# read all genotype data (including variants of all variant sets and all traits)
genofile_rds =".../ukb_imp_v3_SNPs_maf0.01_ldthin0.5_full.rds"
obj.bigSNP <- snp_attach(genofile_rds)
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))  - 1



#match variants in saummary stats with variants in genotype data
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map,strand_flip=FALSE)



# remove unqualified varaints based on the recommendation in LDpred2 paper
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_y_est = min( sd_ldref * info_snp$beta_se * sqrt(info_snp$n_eff))
sd_ss = with(info_snp, sd_y_est / sqrt(n_eff * beta_se^2))
is_bad <-  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05
info_snp = info_snp[!is_bad,]


# folder stores the correlation matrix between all the variants in Gb by chromosome (pre-trained using 10,000 EU samples in UKB)
corr_dir = ".../full_bed_corr/"

#### create a tempfolder to save temp files
## use ramdisk to speed up I/Os on a cluster node
mainDir="/ramdisks/"
subDir = paste0(trait,"_temp")
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
temp_folder = paste0(mainDir,subDir)
tmp <- tempfile(tmpdir = temp_folder)


# correlation matrices were stored by chromosome.
# Extract correlation matrices of variants by chromosome (given the trait and the variant set) and combine them to a single correlation matrix
for (chr in 1:22) {
  print(paste0("Reading corr matrix of chr:", chr,' -- ' ,Sys.time()))
  
  ## indices in 'sumstats'
  ind.chr <- which(info_snp$chr == chr)
  
  # make sure there is at least one variant for a chromosome
  if(length(ind.chr)>0) {
  
      ## indices in 'G'
      ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
      ## indices in 'corr'
      ind.chr3 <- match(ind.chr2, which(CHR == chr))
      
      corr0 <- readRDS(paste0(corr_dir,"chr", chr, ".rds"))[ind.chr3, ind.chr3]
      #print(paste0("finished rds reading:",Sys.time()))
      
      
      # the first chromosome
      if (chr == 1) {
        df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
        df_model_beta = info_snp[ind.chr, c('rsid','chr','pos','a0','a1')] 
        
        # if only one variant
        if (length(ind.chr)==1){
          ld = c(1)
        }
        else {
        ld <- Matrix::colSums(corr0^2)
        }
     
        # if only one variant
        if (length(ind.chr)==1){
          corr = as_SFBM(as.matrix(corr0), tmp)
        }
        else {
          corr = as_SFBM(corr0, tmp)
        }
        
      } else {
        df_beta <- rbind(df_beta, info_snp[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")])
        df_model_beta <- rbind(df_model_beta, info_snp[ind.chr, c('rsid','chr','pos','a0','a1')])  
        
        if (length(ind.chr)==1){
          ld = c(ld, 1)
        }
        else {
          ld <- c(ld, Matrix::colSums(corr0^2))
        }
        
        if (length(ind.chr)==1){
          corr$add_columns(as(as.matrix(corr0), "dgCMatrix"), nrow(corr))
        }
        else {
          corr$add_columns(corr0, nrow(corr))
        }
      }
  
  }
}

#intercept=1 - need to set default intercept (to be estimated in the function) for variants sets with highly related variants!
#but you will need to set intercept=1 if you are using highly independent variant sets, e.g. CA variant sets
(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,sample_size = n_eff,blocks = NULL)))
h2_est <- ldsc[["h2"]]


# LDpred2-grid parameters
h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

# running ldpred2 with genome-wide variants
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)


#create save folder for the given variant set and given trait
save_mainfolder = ".../training_models/"
dir.create(file.path(save_mainfolder, vars_sele_folder), showWarnings = FALSE)
save_mainfolder = paste0(save_mainfolder,vars_sele_folder)
dir.create(file.path(save_mainfolder, trait), showWarnings = FALSE)

# save tested parameter sets 
save_file_paras= paste0(save_mainfolder,'/',trait,'/','params.txt')
write.csv(params,save_file_paras,row.names = FALSE, quote=FALSE)

# loop through each set of betas trained and save as a pgs model
for(i in 1:dim(params)[1]) {
  pgs_model = data.frame(df_model_beta)
  colnames(pgs_model)[4:5]= c('effect_allele','other_allele')
  col = beta_grid[,i]
  pgs_model$effect = col
  save_file_model = paste0(save_mainfolder,'/',trait,'/',"pgs_model_",as.character(i),".txt")
  write.table(pgs_model,save_file_model,sep="\t",row.names = FALSE, quote=FALSE)
}

### LDpred2 dosen't support genotype data with missing variants when calculating genetic scores with a PGS model (snp_ldpred2_grid method);
### Also, the filesystem on our cluster doesn't support frequent matadata operations (memory mapping in LDpred2), i.e. extremely slow with snp_fastImputeSimple();
### so we seperately used Plink2 to calculate PGS scores with these trained PGS models above,and selected the best model on validation data (then do internal and external validation)

unlink(temp_folder,recursive = TRUE)





