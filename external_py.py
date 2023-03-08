import numpy as np
from scipy import stats
from sklearn import linear_model as lm


def _get_model_info(model, Z, y):
    """
    Helper function to get fitted coefficients, R2, and log-likelihood
    """
    n, p = Z.shape
    coef = model.coef_

    r2 = model.score(Z, y)
    ystar = model.predict(Z)
    s2e = sum((y - ystar) ** 2) / (n - 1)

    logl = sum(stats.norm.logpdf(y, loc=ystar, scale=np.sqrt(s2e)))

    return coef, r2, logl


def fit(Z, y, h2g, b_qtls=None, args=None):
    n, p = Z.shape

    model = lm.LinearRegression()
    model.fit(Z, y)

    coef, r2, logl = _get_model_info(model, Z, y)

    return coef, r2, logl
