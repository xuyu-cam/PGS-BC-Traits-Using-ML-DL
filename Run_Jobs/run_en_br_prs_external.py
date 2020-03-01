

from Methods.Traditional_GRS import traditional_GRS_selected_vars
import datetime,os
from sklearn.model_selection import KFold
from sklearn.externals import joblib
import numpy as np
import pandas as pd

from Utls.read_traits import get_trait_abbrs
from Utls.read_trait_genotypes import read_trait_genotypes
from Utls.read_variants_ids import read_variant_ids
from Utls.read_adjusted_trait_levels import get_trait_vector
from Utls.get_QCed_sample_index import get_valid_sample_index
from Utls.read_revelant_sample_ids import read_related_sampleIDs
from Utls.get_QCed_non_relevant_sample_ids import get_valid_sample_relevant_id_filtering
from Utls.get_model_prediction import get_prediction_measures


def run_experiments_5_folders(trait_abbr_file,xdata_path,variants_path,beta_path,model_path,results_file,ukb_traits_value_file,relevant_sampleIDs_file,Remove_relavated_samples):
    # Testing the performance of the trained EN and BR model (internal test) and the univariant method on the INTERVAL data (external test)
    # Remove_relavated_samples: FALSE OR TRUE

    trait_abbrs = get_trait_abbrs(trait_abbr_file) # Trait list that includes all the trait names

    #read sample ids that are related to samples in UKB
    relevant_sample_ids = read_related_sampleIDs(relevant_sampleIDs_file)

    kf = KFold(n_splits=5, shuffle=True, random_state=21)
    with open(results_file,'w') as f:
        f.write("TRAIT\tTot_Samples\tFolder\tN_of_Vars\tPRS_R2\tPRS_R\tElasticNet_R2\tElasticNet_R\tBayesianRidge_R2\tBayesianRidge_R\n")
        for trait in trait_abbrs:

            print("{} - Start processing trait: {}".format(datetime.datetime.now(),trait))

            # read conditional analysis variants ids in a list - in the order of them listed in the file
            var_ids = read_variant_ids(trait, variants_path)

            # read samples IDs in a list, and sample-genotype X data matrix
            sample_ids_pre,X = read_trait_genotypes(trait,xdata_path,var_ids)

            # read a vector of adjusted trait values against the X data matrix
            y_pre = get_trait_vector(trait,ukb_traits_value_file)

            if Remove_relavated_samples == False:  # choose whether remove these relevant samples
                valid_index = get_valid_sample_index(y_pre)   # non-"NA" samples are kept only
            else:
                valid_index = get_valid_sample_relevant_id_filtering(y_pre,sample_ids_pre,relevant_sample_ids)  # only non-"NA" and non-revelant samples are kept

            # get the final valid genotype data X and and trait value data y
            X = X[valid_index,:]
            y = np.array([float(y_pre[i]) for i in valid_index])
            sample_ids = [sample_ids_pre[i] for i in valid_index]


            folder_count = 0

            #run external test for 5 times using the 5 trained models from UKB
            for train_index, test_index in kf.split(X):
                folder_count += 1
                print("folder-{}".format(folder_count))

                model_file = model_path + trait + "_EN_model_" + str(folder_count) + ".pkl"
                print("model file: {}".format(model_file))
                en_model = joblib.load(model_file) # Load a EN model
                EN_r2, EN_r = get_prediction_measures(en_model, X, y)
                print("R2: {}, r: {}".format(EN_r2, EN_r))


                print("{}-folder-{} Running BayesianRidge...".format(datetime.datetime.now(), folder_count))
                model_file = model_path + trait + "_BR_model_" + str(folder_count) + ".pkl"
                print("model file: {}".format(model_file))
                br_model = joblib.load(model_file) # Load a BR model
                BR_r2, BR_r = get_prediction_measures(br_model, X, y)
                print("R2: {}, r: {}".format(BR_r2, BR_r))


                print("{}-folder-{} Running PRS...".format(datetime.datetime.now(), folder_count))
                beta_file = beta_path + trait + "_beta"
                GRS_r2, GRS_r = traditional_GRS_selected_vars(beta_file, X, y,var_ids)
                print("R2: {}, r: {}".format(GRS_r2, GRS_r))

                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(trait,X.shape[0],folder_count,X.shape[1], GRS_r2,GRS_r, EN_r2, EN_r, BR_r2, BR_r))
                f.flush()


if __name__ == "__main__":

    # File lists all the trait names
    trait_abbr_file = "/TRAIT_MAP.tsv"

    # Folder path under which the INTERVAL genotype data of each trait is stored by file
    data_path = "/condsig_variants_com/"

    # Folder path under which the GWAS betas of each trait (conditonal analysis variants) is stored by file
    beta_path = "/condsig_variants_com_beta/"

    # File that store the results of external validation of EN, BR and UNI
    results_file = "/interval_results_UNI_EN_BR_condsig_com.txt"

    # Folder path under which all the traied EN and BR models (using UNB data) were stored
    model_path = "/EN_BR_condsig_variants_com/"

    # Folder path under which conditonal analysis variants of each trait is stored by file
    variants_path = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/ukb_blood_traits_genetics/5k_vars_com_ukb_interval/variants/condsig_variants_com/"

    # Folder path under which adjusted trait values for each trait in INTERVAL is stored by file
    interval_traits_value_file = "/int_170k.sample"

    # Identified samples in INTERVAL that are releted to samples in UKB
    relevant_sampleIDs_file = "/home/yx322/tests/interval_related_UKBB.txt"

    run_experiments_5_folders(trait_abbr_file,data_path,variants_path,beta_path,model_path,results_file,interval_traits_value_file,relevant_sampleIDs_file,Remove_relavated_samples=False)