
from sklearn.model_selection import train_test_split
import random
from sklearn.linear_model import SGDRegressor
import numpy as np
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from sklearn.linear_model import ElasticNetCV
import datetime
import logging



def full_fit_ElasticNetCV_GD(X_train, X_test, y_train, y_test):
    # Elastic Net Regression with 10 folder cross validation

    # Note that primilary analyses were conducted by setting l1_ratio as [0.1,0.2, 0.5, 0.7, 0.9, 0.95, 0.99, 1] on 20% the training data (for the consideration of runing time)
    # and 0.5 was found the best or close to the best perfomrnace on the remaing traitning data (all the traits). Thus, the defult 0.5 was applied for l1_ratio in this study.
    # X_train: training genotype data  (Numpy matrix - samples X viriants)
    # X_test: testing genotype data
    # y_train: training trait value data (Numpy vector - 1 X N)
    # y_test: testing trait value data
    # return the learned model, r and r2 performance of the model on testing data

    regr = ElasticNetCV(cv=10, random_state=0,n_jobs=-1,max_iter =1500)
    regr.fit(X_train, y_train)
    y_pred = regr.predict(X_test)
    # print("ElasticNet with full training data load - R2 result: ", r2_score(y_test, y_pred))
    # print("ElasticNet with full training data load - Pearson R result: ", pearsonr(y_test, y_pred))
    return regr,r2_score(y_test, y_pred),pearsonr(y_test, y_pred)[0]


# if __name__ == "__main__":
#     X = np.load("mpv_xdata.npy")
#     y = np.load("mpv_ydata.npy")
#     x_train, x_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=5)
#     full_fit_ElasticNetCV_GD(x_train, x_test, y_train, y_test)
