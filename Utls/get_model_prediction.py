from sklearn.metrics import r2_score
from scipy.stats import pearsonr

def get_prediction_measures(model,X,y):
    y_pred = model.predict(X)
    return r2_score(y, y_pred), pearsonr(y, y_pred)[0]
