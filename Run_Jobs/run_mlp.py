

from Methods.MLP import mlp_train
import datetime,os,sys
from sklearn.model_selection import KFold
from sklearn.externals import joblib
import numpy as np
import pandas as pd
from Utls.data_coding_converter import convert_to_one_hot


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

# function that gets the pheno data by a vector
def get_trait_vector(pheno_name,pheno_file_path):
    col_name_index_map = {}
    pheno_vec = []
    with open(pheno_file_path) as f:
        col_name_line = f.readline().strip().split()
        for i in range(len(col_name_line)):
            col_name_index_map[col_name_line[i]]=i
        target_index = col_name_index_map[pheno_name+ '_gwas_normalised']
   #     f.readline() # skip the null line
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


def run_experiments_5_folders(trait_abbr_file,xdata_path,variants_path,networks_path,model_path,results_file,ukb_traits_value_file,trait,one_hot_coding):

    #trait_abbrs = get_trait_abbrs(trait_abbr_file) # Trait list that includes all the trait names
    kf = KFold(n_splits=5, shuffle=True, random_state=21)

    with open(results_file,'w') as f:
        f.write("TRAIT\tTot_Samples\tFolder\tN_of_Vars\tMLP_R2\tMLP_R\n")

        print("{} - Start processing trait: {}".format(datetime.datetime.now(),trait))
        var_ids = read_variant_ids(trait, variants_path)
        sample_ids_pre,X = read_trait_genotypes(trait,xdata_path,var_ids)
        y_pre = get_trait_vector(trait,ukb_traits_value_file)

        valid_index = get_valid_sample_index(y_pre)
        X = X[valid_index,:]

        if one_hot_coding==True:
            X = convert_to_one_hot(X)

        y = np.array([float(y_pre[i]) for i in valid_index])
        sample_ids = [sample_ids_pre[i] for i in valid_index]

        folder_count = 0
        for train_index, test_index in kf.split(X):
            folder_count += 1
            print("folder-{}".format(folder_count))
            x_train, x_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            print("{}-folder-{} MLP...".format(datetime.datetime.now(), folder_count))
            network_file = networks_path + trait + "_cond_maf_mlp.txt"
            model, mlp_r2, mlp_r = mlp_train(network_file,x_train,y_train,x_test,y_test)
            print("R2: {}, r: {}".format(mlp_r2, mlp_r))
            model_file = model_path + trait + "_mlp_model_" + str(folder_count) + ".h5"
            model.save(model_file)

            f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(trait,X.shape[0],folder_count,x_train.shape[1], mlp_r2,mlp_r))
            f.flush()


if __name__ == "__main__":

    # trait_abbr_file = "E:\\Blood Cell Trait Data\\TRAIT_MAP.tsv"
    # data_path = "E:\\Blood Cell Trait Data\\xdata_ydata\\"
    # results_file = "E:\\PyCharmProjects\\GenoImpute\\Results\\blood_traits_EN_BR_results_SNPs_condig.txt"
    # run_experiments(trait_abbr_file, data_path, results_file)

    trait_index = int(sys.argv[1])


    trait_abbr_file = "/home/yx322/tests/TRAIT_MAP.tsv"
    trait_abbrs_all = get_trait_abbrs(trait_abbr_file)  # Trait list that includes all the trait names
    trait_abbrs = [trait for trait in trait_abbrs_all if trait not in ['mscv', 'mrv', 'rdw_cv']]
    trait = trait_abbrs[trait_index]

    data_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/condsig_xdata/condsig_xdata_com/"

    result_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/results/ukb_results_MLP_condsig_com_PC_adjusted/"

    results_file = result_path + trait + "_ukb_results_MLP_condsig_com_PC_adjusted.txt"
    model_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/models/MLP_condsig_variants_com_PC_adjusted/"
    variants_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/variants/condsig_variants_com/"
    ukb_traits_value_file ="/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/data/ukb_blood_cell_trait_eu_PC_adjusted/ukbb_500k.sample"
    networks_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/nn_structures/MLP_Results/"

    if os.path.isdir(model_path) == False:
           os.mkdir(model_path)
    if os.path.isdir(result_path) == False:
        os.mkdir(result_path)
    # run_experiments_5_folders(trait_abbr_file, data_path, beta_path, model_path, results_file)

    run_experiments_5_folders(trait_abbr_file,data_path,variants_path,networks_path,model_path,results_file,ukb_traits_value_file,trait,False)