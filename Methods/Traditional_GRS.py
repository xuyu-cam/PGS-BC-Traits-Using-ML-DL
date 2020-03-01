
import numpy as np
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from sklearn.model_selection import train_test_split
import pandas as pd

def get_beta_vec(beta_file):
    #get GWAS effect sizes of all the conditional analysis variants

    beta_vec = []
    with open(beta_file) as f:
        for line in f:
            beta = float(line.strip().split(' ')[-1])
            beta_vec.append(beta)
    return np.array(beta_vec)


def get_beta_vec_by_vars_ids(beta_file,vars):
    # get GWAS effect sizes of the selected variants

    beta_vec = []
    df = pd.read_csv(beta_file,header=None,sep=' ',index_col=0)
    betas = df.loc[vars][1]
    return np.array(betas)


def traditional_GRS(beta_file,X,y):
    # Univariate PGS method (all the conditional analysis variants)
    # beta_file: stores the GWAS effect sizes of variants by row (Order of the variants are set the same as test data X in the previous analysis)
    # X: testing genotype data
    # y: testing trait value data

    beta_vec = get_beta_vec(beta_file)
    beta_vec = -beta_vec  # flip the sign of the effect sizes because the GWAS were conducted against the reference alleles
    y_pred = X.dot(beta_vec)
    # print(y_pred[:10],len(y_pred))

    y_pred = y_pred/X.shape[1]
    return r2_score(y, y_pred),pearsonr(y, y_pred)[0]

def traditional_GRS_selected_vars(beta_file,X,y,vars):
    # Univariate PGS method (only using given variants)
    beta_vec = get_beta_vec_by_vars_ids(beta_file,vars)
    beta_vec = -beta_vec
    y_pred = X.dot(beta_vec)
    # print(y_pred[:10],len(y_pred))
    y_pred = y_pred/X.shape[1]
    return r2_score(y, y_pred),pearsonr(y, y_pred)[0]

def traditional_GRS_nk(beta_file, k,X, y):
    # Univariate PGS method (Uing k*1000 variants)
    beta_vec = get_beta_vec(beta_file)[:k*1000]
    beta_vec = -beta_vec
    y_pred = X.dot(beta_vec)
    # print(y_pred[:10],len(y_pred))
    y_pred = y_pred / X.shape[1]
    return r2_score(y, y_pred), pearsonr(y, y_pred)[0]






# if __name__ == "__main__":
#
#
#     X = np.load("mpv_xdata.npy")
#     y = np.load("mpv_ydata.npy")
#     beta_file = "mpv_gwas_variants"
#     x_train, x_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=6)
#     print(traditional_GRS2(beta_file, x_test, y_test))
#
