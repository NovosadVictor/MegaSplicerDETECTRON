import os
import re
from functools import partial
import matplotlib.pyplot as plt
import numpy as np


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
