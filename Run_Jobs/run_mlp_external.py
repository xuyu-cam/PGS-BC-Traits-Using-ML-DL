

from Methods.ElasticNet import full_fit_ElasticNetCV_GD
from Methods.BayesianRidge import full_fit_BayesianRidge
from Methods.Traditional_GRS import traditional_GRS
import datetime,os
from sklearn.model_selection import KFold
from sklearn.externals import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from keras.models import load_model
from Utls.data_coding_converter import convert_to_one_hot
from Utls.read_traits import get_trait_abbrs
from Utls.read_trait_genotypes import read_trait_genotypes
from Utls.read_variants_ids import read_variant_ids
from Utls.read_adjusted_trait_levels import get_trait_vector
from Utls.get_QCed_sample_index import get_valid_sample_index
from Utls.read_revelant_sample_ids import read_related_sampleIDs
from Utls.get_QCed_non_relevant_sample_ids import get_valid_sample_relevant_id_filtering


def get_prediction_measures(model,X,y):
    r = pearsonr(y,model.predict(X).ravel())[0]
    r2 = r2_score(y, model.predict(X).ravel())
    return r2,r


def run_experiments_5_folders(trait_abbr_file,xdata_path,variants_path,model_path,results_file,ukb_traits_value_file,relevant_sampleIDs_file,Remove_relavated_samples,one_hot_coding):
    # Testing the performance of the trained MLP model (internal test) on the INTERVAL data

    trait_abbrs = get_trait_abbrs(trait_abbr_file) # Trait list that includes all the trait names
    relevant_sample_ids = read_related_sampleIDs(relevant_sampleIDs_file)
    kf = KFold(n_splits=5, shuffle=True, random_state=21)

    with open(results_file,'w') as f:
        f.write("TRAIT\tTot_Samples\tFolder\tN_of_Vars\tMLP_R2\tMLP_R\n")
        for trait in trait_abbrs:
            print("{} - Start processing trait: {}".format(datetime.datetime.now(),trait))
            var_ids = read_variant_ids(trait, variants_path)  # read variant ids in a list - in the order of them listed in the file
            sample_ids_pre,X = read_trait_genotypes(trait,xdata_path,var_ids) # read samples IDs in a list, and sample-genotype X data matrix
            y_pre = get_trait_vector(trait,ukb_traits_value_file) # read a vector of trait values against the X data matrix

            if Remove_relavated_samples == False:
                valid_index = get_valid_sample_index(y_pre)   # non-"NA" samples are kept only
            else:
                valid_index = get_valid_sample_relevant_id_filtering(y_pre,sample_ids_pre,relevant_sample_ids)  # only non-"NA" and non-revelant samples are kept

            X = X[valid_index,:]
            if one_hot_coding == True:
                X = convert_to_one_hot(X)

            y = np.array([float(y_pre[i]) for i in valid_index])
            sample_ids = [sample_ids_pre[i] for i in valid_index]

            folder_count = 0
            for train_index, test_index in kf.split(X):
                folder_count += 1
                print("folder-{}".format(folder_count))
                model_file = model_path + trait + "_mlp_model_" + str(folder_count) + ".h5"
                print("model file: {}".format(model_file))
                mlp_model = load_model(model_file)
                mlp_r2, mlp_r = get_prediction_measures(mlp_model, X, y)
                print("R2: {}, r: {}".format(mlp_r2, mlp_r))

                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(trait,X.shape[0],folder_count,X.shape[1], mlp_r2,mlp_r))
                f.flush()


if __name__ == "__main__":

    trait_abbr_file = "/TRAIT_MAP.tsv"
    data_path = "/condsig_variants_com/"
    results_file = "/interval_results_MLP_condsig_com.txt"

    # Folder path under which all the trained CNN models were stored for each trait
    model_path = "/models/MLP_condsig_variants_com/"
    variants_path = "/variants/condsig_variants_com/"
    ukb_traits_value_file = "/int_170k.sample"
    relevant_sampleIDs_file = "/interval_related_UKBB.txt"

    run_experiments_5_folders(trait_abbr_file,data_path,variants_path,model_path,results_file,ukb_traits_value_file,relevant_sampleIDs_file,Remove_relavated_samples=True,one_hot_coding=False)