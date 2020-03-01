import pandas as pd
import numpy as np




def read_trait_genotypes(trait,xdata_path,var_ids):
    #get genotype data of a trait
    #trait: trait name
    #xdata_path: path where the genotype data is stored; genotype data of a trait is stored with a gzipped csv file in which the data entry is seperated using ','
    #var_ids: the list of variants ids used
    #returning the used sample ids and a matrix of samples X variants (snps)
    geno_file = xdata_path + trait + '_xdata.csv.gz'
    df = pd.read_csv(geno_file, sep=',',compression='gzip')
    sample_ids = list(df["sample_ids"])
    df = df[var_ids]
    return sample_ids,np.array(df)



