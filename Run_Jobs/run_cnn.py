

from Methods.CNN import cnn_train
import datetime,os,sys
from sklearn.model_selection import KFold
import numpy as np
from Utls.data_coding_converter import convert_to_one_hot
from Utls.read_traits import get_trait_abbrs
from Utls.read_trait_genotypes import read_trait_genotypes
from Utls.read_variants_ids import read_variant_ids
from Utls.read_adjusted_trait_levels import get_trait_vector
from Utls.get_QCed_sample_index import get_valid_sample_index


def run_experiments_5_folders(xdata_path,variants_path,networks_path,model_path,results_file,ukb_traits_value_file,trait,one_hot_coding):
    # Partition the UKB data into 5 folders, and train 5 CNN PGS models on any 4 folders of the data,
    # and testing the performance of the trained model (internal test) on the remaining folder


    # 5 folder partition
    kf = KFold(n_splits=5, shuffle=True, random_state=21)
    with open(results_file,'w') as f:
        f.write("TRAIT\tTot_Samples\tFolder\tN_of_Vars\tMLP_R2\tMLP_R\n")
        print("{} - Start processing trait: {}".format(datetime.datetime.now(),trait))

        # read conditional analysis variants of a trait
        var_ids = read_variant_ids(trait, variants_path)

        # read genotype data (matrix: Samples X variants), and the corresponding sample ids
        sample_ids_pre,X = read_trait_genotypes(trait,xdata_path,var_ids)

        # get adjusted trait values of a trait
        y_pre = get_trait_vector(trait,ukb_traits_value_file)

        # Get valid indexes of samples (trait values of excluded samples are set as 'NA' in the stored file)
        valid_index = get_valid_sample_index(y_pre)

        X = X[valid_index,:]

        # one-hot-encoding converting
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

            print("{}-folder-{} CNN...".format(datetime.datetime.now(), folder_count))
            network_file = networks_path + trait + "_cond_maf_cnn.txt"
            model, cnn_r2, cnn_r = cnn_train(network_file,x_train,y_train,x_test,y_test)
            print("R2: {}, r: {}".format(cnn_r2, cnn_r))
            model_file = model_path + trait + "_cnn_model_" + str(folder_count) + ".h5"
            model.save(model_file)


            f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(trait,X.shape[0],folder_count,x_train.shape[1], cnn_r2,cnn_r))
            f.flush()


if __name__ == "__main__":

    #choose a trait by an index
    trait_index = int(sys.argv[1])

    # File lists all the trait names
    trait_abbr_file = "/TRAIT_MAP.tsv"
    trait_abbrs_all = get_trait_abbrs(trait_abbr_file)  # Trait list that includes all the trait names

    trait = trait_abbrs_all[trait_index]

    # Folder path under which the UKB genotype data of each trait is stored by file
    data_path = "/condsig_xdata_com/"

    # Path where the result files stored
    result_path = "/ukb_results_CNN_condsig_com/"

    #results file
    results_file = result_path + trait + "_ukb_results_CNN_condsig_com.txt"

    #Folder path under which all the traied CNN models will be stored
    model_path = "/CNN_condsig_variants_com/"

    #Folder path under which conditonal analysis variants of each trait is stored by file
    variants_path = "/condsig_variants_com/"

    # Folder path under which adjusted trait values of each trait in UKB is stored by file
    ukb_traits_value_file = "/home/yx322/rds/rds-jmmh2-projects/inouye_lab_other/impute_genomics/yx322/data/ukb_blood_cell_trait_eu_PC_adjusted/ukbb_500k.sample"

    #Folder path under which all the selected top 10 CNN network structures of each trait were stored by txt files
    #Note that, based on the protocal described in https://github.com/paubellot/DL-Biobank, we selected the top 10 CNN models via a validation step on the training data.
    networks_path = "/CNN_structures/"

    if os.path.isdir(model_path) == False:
           os.mkdir(model_path)
    if os.path.isdir(result_path) == False:
        os.mkdir(result_path)

    # one_hot_coding: choose whether to convert the genotype X data to one-hot encoding representation
    run_experiments_5_folders(data_path,variants_path,networks_path,model_path,results_file,ukb_traits_value_file,trait,one_hot_coding=False)