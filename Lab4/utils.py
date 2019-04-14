import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl


import statsmodels.formula.api as smf
import statsmodels.tsa.api as smt
import statsmodels.api as sm
import scipy.stats as scs

def multiplot(y, lags=None, figsize=(16, 12)):
    """
    Takes a Pandas serie and plots the Autocorrelation, partial autocorrelation, normal qq plot and probability plot.
    """
    with plt.style.context('bmh'):    
        fig, (ax1, ax2) = plt.subplots(2,2, figsize=figsize)
        smt.graphics.plot_acf(y, lags=lags, ax=ax1[0], alpha=0.05)
        smt.graphics.plot_pacf(y, lags=lags, ax=ax1[1], alpha=0.05)
        sm.qqplot(y, line='s', ax=ax2[0])
        ax2[0].set_title('Normal QQ Plot')        
        scs.probplot(y, sparams=(y.mean(), y.std()), plot=ax2[1])

        plt.tight_layout()
    return 