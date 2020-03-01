from sklearn.model_selection import train_test_split
import random
from sklearn.linear_model import SGDRegressor
import numpy as np
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from sklearn.linear_model import BayesianRidge
from sklearn.linear_model import BayesianRidge




def full_fit_BayesianRidge(x_train, x_test, y_train, y_test):
    # Bayesian Ridge Regression Method
    #
    # Note prior Gamma distribution are set as (100000, -10) and (-10, 0) which were selected via cross-validation step using training data (see codes below: full_fit_BayesianRidge_para_turing;
    # all traits shared the same)
    # X_train: training genotype data  (Numpy matrix - samples X viriants)
    # X_test: testing genotype data
    # y_train: training trait value data (Numpy vector - 1 X N)
    # y_test: testing trait value data
    # return the learned model, r and r2 performance

    model, r2 = get_BayesianRidge_prediction(x_train, y_train, x_test, y_test, 100000, -10, -10, 0)
    y_pred = model.predict(x_test)
    return model,r2_score(y_test, y_pred),pearsonr(y_test, y_pred)[0]



#Bayesian Ridge Regression Method, given a set of prior Gamma distribution
def get_BayesianRidge_prediction(x_train, y_train, x_val, y_val, alpha_1, alpha_2, lambda_1, lambda_2):
    model = BayesianRidge(alpha_1=alpha_1, alpha_2=alpha_2, lambda_1=lambda_1, lambda_2=lambda_2)
    model.fit(x_train, y_train)
    y_pred = model.predict(x_val)
    return model,r2_score(y_val, y_pred)


#hyper-parameter Tuning - finding the best prior Gamma distributions on the training set of a trait
#Grid search on (1e10, 1e5, 1e1, 0, -1e1, -1e5, -1e10)
#return the best 'alpha_1', 'alpha_2' 'lambda_1' 'lambda_2'
def full_fit_BayesianRidge_para_turing(x_train, x_val, y_train, y_val):
    nums = (1e10, 1e5, 1e1, 0, -1e1, -1e5, -1e10)
    alpha_1 = nums
    alpha_2 = nums
    lambda_1 = nums
    lambda_2 = nums
    best_model = None
    best_r2 = 0
    for a1 in alpha_1:
        for a2 in alpha_2:
            for l1 in lambda_1:
                for l2 in lambda_2:
                    # print("Training BayesianRidge with alpha_1: {}, alphs_2: {}, lambda_1: {}, lambda_2:{}".format(a1,a2,l1,l2))
                    model,r2 = get_BayesianRidge_prediction(x_train,y_train,x_val,y_val,a1,a2,l1,l2)
                    if best_model == None or r2 > best_r2:
                        best_model = model
                        best_r2 = r2
                        best_params =  {'alpha_1':a1, 'alpha_2': a2, 'lambda_1': l1, 'lambda_2': l2}
    print("Best Para: {}".format(best_params))
    return best_model



#if __name__ == "__main__":
   # X = np.load("mpv_xdata.npy")
   # y = np.load("mpv_ydata.npy")
   # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=5)
   # X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, random_state=42)
   # model, r2 = get_BayesianRidge_prediction(X_train, y_train, X_test, y_test, 100000, -10, -10, 0)

