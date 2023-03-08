"""
This file is an example of how to define an external/custom function to fit a
predictive model of gene expression from genotype to be used by `twas_sim`. Here
we are using the `LinearRegression` class from `sklearn` to illustrate external
python functionality.

External modules -must- include a function named `fit` that takes as arguments:
    Z: numpy matrix of genotype
    y: numpy vector of gene expression/phenotype
    h2g: the true h2g of gene expression
    b_qtls: the true beta/effect-sizes for gene expression (i.e. eQTL)
    args: the argparse object from twas_sim; useful for pulling `args.output`
        as a prefix for temp data.

Similarly, it must return a tuple containing (coef, r2, logl):
    coef: the numpy vector for estimated eQTL weights
    r2: the predictive r2 (optional; None)
    logl: the log likelihood of the model (optional; None)

"""
import numpy as np
from scipy import stats
from sklearn import linear_model as lm


def fit(Z, y, h2g, b_qtls=None, args=None):
    n, p = Z.shape

    model = lm.LinearRegression()
    model.fit(Z, y)

    coef = model.coef_

    r2 = model.score(Z, y)
    ystar = model.predict(Z)
    s2e = sum((y - ystar) ** 2) / (n - 1)

    logl = sum(stats.norm.logpdf(y, loc=ystar, scale=np.sqrt(s2e)))

    return coef, r2, logl
