import os
import re
from functools import partial
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, mannwhitneyu
from sklearn.metrics import r2_score


def save_plt_fig(name, format):
    if format == 'tiff':
        kwargs = {'compression': 'tiff_lzw'} if format == 'tiff' else None
        plt.savefig(name, format=format, pil_kwargs=kwargs, dpi=350)
    else:
        plt.savefig(name, format=format, dpi=350)


def getattr_with_kwargs(module, method):
    if isinstance(method, dict):
        return partial(getattr(module, method['name']), **method.get('kwargs', {}))

    return getattr(module, method)


def intersect_dfs(dfs):
    common_index = set(dfs[0].index).intersection(*map(lambda df: df.index, dfs[1:]))
    return [df.loc[common_index] for df in dfs]


def make_sure_dir_exists(dir_path):
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)

    return dir_path


def logit(x):
    return np.log(x / (1 - x))


def inlogit(y):
    return np.exp(y) / (1 + np.exp(y))


def find_substring_occurrences(substring, string):
    return [m.start() for m in re.finditer(substring, string)]


def predict(df, coefs, logit=True):
    cols = set(df.columns) & set(coefs.index)
    df_pred = df[cols].mul(coefs.loc[cols, 'Estimate'], axis=1).sum(axis=1)
    if '(Intercept)' in coefs.index:
        df_pred += coefs.loc['(Intercept)']['Estimate']
    return df_pred if logit else inlogit(df_pred)


def get_accuracy(coefs, df):
    predictions = predict(df, coefs)
    return get_scores(predictions, df['fraq'])


def get_scores(pred, true):
    cor = pearsonr(pred, logit(true))[0]
    r2 = max(0, r2_score(pred, logit(true)))
    mann_w = mannwhitneyu(inlogit(pred), true)[1]
    #
    mean_pred, mean_true = np.median(inlogit(pred)), np.median(true)
    uplift = (mean_pred - mean_true) / mean_true
    uplift_score = max(0, 1 - abs(uplift))
    #
    return {
        'cor': cor,
        'r2': r2,
        'mann-w': mann_w,
        'uplift': uplift,
        'uplift.score': uplift_score,
        'mean_pred': mean_pred,
        'mean_true': mean_true,
    }
