

from Methods.ElasticNet import full_fit_ElasticNetCV_GD
from Methods.BayesianRidge import full_fit_BayesianRidge
from Methods.Traditional_GRS import traditional_GRS_selected_vars
import datetime,os
from sklearn.model_selection import KFold
from sklearn.externals import joblib
import numpy as np
from Utls.read_traits import get_trait_abbrs
from Utls.read_trait_genotypes import read_trait_genotypes
from Utls.read_variants_ids import read_variant_ids
from Utls.read_adjusted_trait_levels import get_trait_vector
from Utls.get_QCed_sample_index import get_valid_sample_index



def run_experiments_5_folders(trait_abbr_file,xdata_path,variants_path,beta_path,model_path,results_file,ukb_traits_value_file):
    # Partition the UKB data into 5 folders, and train 5 BR/EN PGS models on any 4 folders of the data,
    # and testing the performance of the trained model (internal test) and the univariant method on the remaining folder


    trait_abbrs = get_trait_abbrs(trait_abbr_file) # Trait list that includes all the trait names
    kf = KFold(n_splits=5, shuffle=True, random_state=21)

    with open(results_file,'w') as f:
        f.write("TRAIT\tTot_Samples\tFolder\tN_of_Vars\tPRS_R2\tPRS_R\tElasticNet_R2\tElasticNet_R\tBayesianRidge_R2\tBayesianRidge_R\n")

        for trait in trait_abbrs:
            print("{} - Start processing trait: {}".format(datetime.datetime.now(),trait))

            #read conditional variants of a trait
            var_ids = read_variant_ids(trait, variants_path)

            #read genotype data (matrix: Samples X variants), and the corresponding sample ids
            sample_ids_pre,X = read_trait_genotypes(trait,xdata_path,var_ids)

            #get adjusted trait values of a trait
            y_pre = get_trait_vector(trait,ukb_traits_value_file)

            #Get valid indexes of samples (trait values of excluded samples are set as 'NA' in the stored file)
            valid_index = get_valid_sample_index(y_pre)

            #get the final valid genotype data X and and trait value data y
            X = X[valid_index,:]
            y = np.array([float(y_pre[i]) for i in valid_index])
            sample_ids = [sample_ids_pre[i] for i in valid_index]

            folder_count = 0
            # partioning the samples to 5 folders
            for train_index, test_index in kf.split(X):
                folder_count += 1
                print("folder-{}".format(folder_count))

                #Use 4 folders as training data and the remaining one folder as testing data (internal test)
                x_train, x_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]

                print("{}-folder-{} Running ElasticNet...".format(datetime.datetime.now(), folder_count))
                EN_model, EN_r2, EN_r = full_fit_ElasticNetCV_GD(x_train, x_test, y_train, y_test)
                print("R2: {}, r: {}".format(EN_r2, EN_r))
                model_file = model_path + trait + "_EN_model_" + str(folder_count) + ".pkl"
                joblib.dump(EN_model, model_file)  # save trained EN model

                print("{}-folder-{} Running BayesianRidge...".format(datetime.datetime.now(), folder_count))
                BR_model, BR_r2, BR_r = full_fit_BayesianRidge(x_train, x_test, y_train, y_test)
                print("R2: {}, r: {}".format(BR_r2, BR_r))
                model_file = model_path + trait + "_BR_model_" + str(folder_count) + ".pkl"
                joblib.dump(BR_model, model_file) # save trained BR model

                print("{}-folder-{} Running traditional PRS...".format(datetime.datetime.now(), folder_count))
                beta_file = beta_path + trait + "_beta"
                GRS_r2, GRS_r = traditional_GRS_selected_vars(beta_file, x_test, y_test,var_ids)
                print("R2: {}, r: {}".format(GRS_r2, GRS_r))

                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(trait,X.shape[0],folder_count,x_train.shape[1], GRS_r2,GRS_r, EN_r2, EN_r, BR_r2, BR_r))
                f.flush()


if __name__ == "__main__":

    # File lists all the trait names
    trait_abbr_file = "TRAIT_MAP.tsv"

    #Folder path under which the UKB genotype data of each trait is stored by file
    data_path = "/condsig_xdata_com/"

    #Folder path under which the GWAS betas of each trait (conditonal analysis variants) is stored by file
    beta_path = "/condsig_variants_com_beta/"

    #File that store the results
    results_file = "/results/ukb_results_UNI_BR_EN_condsig_com.txt"

    #Folder path under which all the traied EN and BR models will be stored
    model_path = "/EN_BR_condsig_variants_com/"

    #Folder path under which conditonal analysis variants of each trait is stored by file
    variants_path = "/condsig_variants_com/"

    #Folder path under which adjusted trait values of each trait in UKB is stored by file
    ukb_traits_value_file ="/ukbb_500k.sample"

    if os.path.isdir(model_path) == False:
           os.mkdir(model_path)

    run_experiments_5_folders(trait_abbr_file,data_path,variants_path,beta_path,model_path,results_file,ukb_traits_value_file)