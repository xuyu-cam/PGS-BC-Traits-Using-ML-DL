

from Methods.ElasticNet import full_fit_ElasticNetCV_GD
from Methods.BayesianRidge import full_fit_BayesianRidge
from Methods.Traditional_GRS import traditional_GRS_selected_vars
import datetime,os
from sklearn.model_selection import KFold
from sklearn.externals import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from scipy.stats import pearsonr


def get_trait_abbrs(trait_abbr_file):
    trait_abbr =[]
    with open(trait_abbr_file) as f:
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            trait_abbr.append(line[0])
    return trait_abbr


def read_trait_genotypes(trait,xdata_path,var_ids):
    geno_file = xdata_path + trait + '_xdata.csv.gz'
    df = pd.read_csv(geno_file, sep=',',compression='gzip')
    sample_ids = list(df["sample_ids"])
    df = df[var_ids]
    return sample_ids,np.array(df)


def read_variant_ids(trait,variants_path):
    variants_file = variants_path + trait + '_condsig'
    df = pd.read_csv(variants_file,header=None)
    return list(df[0])

def read_related_sampleIDs(file):
    sampleIDs = []
    with open(file) as f:
        for line in f:
            sampleIDs.append(int(line.strip()))
    return sampleIDs

# function that gets the pheno data by a vector
def get_trait_vector(pheno_name,pheno_file_path):
    col_name_index_map = {}
    pheno_vec = []
    with open(pheno_file_path) as f:
        col_name_line = f.readline().strip().split()
        for i in range(len(col_name_line)):
            col_name_index_map[col_name_line[i]]=i
        target_index = col_name_index_map[pheno_name+ '_gwas_normalised']
    #    f.readline() # skip the null line
        for line in f:
            line = line.strip().split()
            pheno_vec.append(line[target_index])
    return pheno_vec  #still a string list


def get_valid_sample_index(ydata):
   valid_index = []
   for i in range(len(ydata)): #Filtering NA samples in xdata and ydata
       if ydata[i] != 'NA':
           valid_index.append(i)
   return valid_index

def get_valid_sample_relevant_id_filtering(ydata,sample_ids,relevant_ids):
   valid_index = []
   for i in range(len(ydata)): #Filtering NA samples in xdata and ydata;
       if ydata[i] != 'NA' and sample_ids[i] not in relevant_ids:
           valid_index.append(i)
   return valid_index

def get_prediction_measures(model,X,y):
    y_pred = model.predict(X)
    return r2_score(y, y_pred), pearsonr(y, y_pred)[0]


def run_experiments_5_folders(trait_abbr_file,xdata_path,variants_path,beta_path,model_path,results_file,ukb_traits_value_file,relevant_sampleIDs_file,Remove_relavated_samples):

    trait_abbrs = get_trait_abbrs(trait_abbr_file) # Trait list that includes all the trait names
    relevant_sample_ids = read_related_sampleIDs(relevant_sampleIDs_file)
    kf = KFold(n_splits=5, shuffle=True, random_state=21)

    with open(results_file,'w') as f:
        f.write("TRAIT\tTot_Samples\tFolder\tN_of_Vars\tPRS_R2\tPRS_R\tElasticNet_R2\tElasticNet_R\tBayesianRidge_R2\tBayesianRidge_R\n")
       # f.write("TRAIT\tTot_Samples\tFolder\tN_of_Vars\tPRS_R2\tPRS_R\n")
        for trait in trait_abbrs:

            if trait not in ['mscv','mrv','rdw_cv']:
                print("{} - Start processing trait: {}".format(datetime.datetime.now(),trait))
                var_ids = read_variant_ids(trait, variants_path)  # read variant ids in a list - in the order of them listed in the file
                sample_ids_pre,X = read_trait_genotypes(trait,xdata_path,var_ids) # read samples IDs in a list, and sample-genotype X data matrix
                y_pre = get_trait_vector(trait,ukb_traits_value_file) # read a vector of trait values against the X data matrix

                if Remove_relavated_samples == False:
                    valid_index = get_valid_sample_index(y_pre)   # non-"NA" samples are kept only
                else:
                    valid_index = get_valid_sample_relevant_id_filtering(y_pre,sample_ids_pre,relevant_sample_ids)  # only non-"NA" and non-revelant samples are kept

                X = X[valid_index,:]
                y = np.array([float(y_pre[i]) for i in valid_index])
                sample_ids = [sample_ids_pre[i] for i in valid_index]


                folder_count = 0
                for train_index, test_index in kf.split(X):
                    folder_count += 1
                    print("folder-{}".format(folder_count))


                    model_file = model_path + trait + "_EN_model_" + str(folder_count) + ".pkl"
                    print("model file: {}".format(model_file))
                    en_model = joblib.load(model_file)
                    EN_r2, EN_r = get_prediction_measures(en_model, X, y)
                    print("R2: {}, r: {}".format(EN_r2, EN_r))


                    print("{}-folder-{} Running BayesianRidge...".format(datetime.datetime.now(), folder_count))
                    model_file = model_path + trait + "_BR_model_" + str(folder_count) + ".pkl"
                    print("model file: {}".format(model_file))
                    br_model = joblib.load(model_file)
                    BR_r2, BR_r = get_prediction_measures(br_model, X, y)
                    print("R2: {}, r: {}".format(BR_r2, BR_r))


                    print("{}-folder-{} Running PRS...".format(datetime.datetime.now(), folder_count))
                    beta_file = beta_path + trait + "_beta"
                    GRS_r2, GRS_r = traditional_GRS_selected_vars(beta_file, X, y,var_ids)
                    print("R2: {}, r: {}".format(GRS_r2, GRS_r))

                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(trait,X.shape[0],folder_count,X.shape[1], GRS_r2,GRS_r, EN_r2, EN_r, BR_r2, BR_r))
                    #f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(trait, X.shape[0],folder_count, X.shape[1], GRS_r2,GRS_r))
                    f.flush()


if __name__ == "__main__":

    # File lists all the trait names
    trait_abbr_file = "/TRAIT_MAP.tsv"

    # Folder path under which the INTERVAL genotype data of each trait is stored by file
    data_path = "/condsig_variants_com/"

    # Folder path under which the GWAS betas of each trait (conditonal analysis variants) is stored by file
    beta_path = "/condsig_variants_com_beta/"

    # File that store the results
    results_file = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/interval_results_UNI_only_condsig_com_PC_adjusted.txt"
    model_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/models/EN_BR_condsig_variants_com_PC_adjusted/"

    variants_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/variants/condsig_variants_com/"
    ukb_traits_value_file = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/data/interval_blood_cell_trait_PC_adjusted/int_170k.sample"
    relevant_sampleIDs_file = "/home/yx322/tests/interval_related_UKBB.txt"

    #if os.path.isdir(model_path) == False:
    #       os.mkdir(model_path)
    # run_experiments_5_folders(trait_abbr_file, data_path, beta_path, model_path, results_file)

    run_experiments_5_folders(trait_abbr_file,data_path,variants_path,beta_path,model_path,results_file,ukb_traits_value_file,relevant_sampleIDs_file,Remove_relavated_samples=False)