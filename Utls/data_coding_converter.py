
import numpy as np
from sklearn.preprocessing import OneHotEncoder


def convert_to_one_hot(X):

    X = np.round(X)
    enc = OneHotEncoder()
    enc.fit(X)
    return enc.transform(X).toarray()

