import os
import re
import numpy as np


def logit(x):
    return np.log(x / (1 - x))


def inlogit(y):
    return np.exp(y) / (1 + np.exp(y))


def find_substring_occurrences(substring, string):
    return [m.start() for m in re.finditer(substring, string)]


def make_sure_dir_exists(dir_path):
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
