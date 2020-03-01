import pandas as pd


def read_variant_ids(trait,variants_path):
    #read the full list of conditional variants ids of a given trait
    variants_file = variants_path + trait + '_condsig'
    df = pd.read_csv(variants_file,header=None)
    return list(df[0])


