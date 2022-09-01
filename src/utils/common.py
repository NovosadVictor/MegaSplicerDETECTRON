import sys
import os
import re
import json
from functools import partial
from shutil import copyfile
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


from src.helpers.pipeline import \
    load_rbp_data, \
    load_isoforms, \
    load_rbps, \
    filter_columns_by_expression, \
    load_data, \
    set_variable_exons
from src.helpers.plots import plot_gene_isoforms


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


def load_config_and_input_data(config_path):
    """Load configuration file and input data
    Parameters
    ----------
    config_path : string
        Path to config file (json).
    Returns
    -------
    dict, dict, pd.DataFrame, pd.DataFrame, pd.DataFrame
    """

    print('Loading config...')
    try:
        config_file = open(config_path, 'r')
        config = json.load(config_file)
    except:
        print('Cannot open configuration file', file=sys.stderr)
        sys.exit(1)

    # Paths are absolute or relative to config file
    config_dirname = os.path.dirname(config_path)

    gene = config['gene']
    if os.path.isfile(config.get('rbp_data_path')):
        rbp_df = pd.read_csv(os.path.join(config_dirname, config['rbp_data_path']), index_col=0)
    else:
        rbp_df = load_rbp_data()
    rbp_df = filter_columns_by_expression(
        rbp_df,
        tresh_mean=config.get('rbps_tresh_mean', 1),
        tresh_var=config.get('rbps_tresh_var', 3),
    )
    if os.path.isfile(config.get('isoforms_data_path')):
        isoforms_df = pd.read_csv(os.path.join(config_dirname, config['isoforms_data_path']), index_col=0)
    else:
        isoforms_df = load_isoforms(gene)
    isoforms_df = filter_columns_by_expression(
        isoforms_df,
        tresh_mean=config.get('isoforms_tresh_mean', 1),
        tresh_var=config.get('isoforms_tresh_var', 10),
    )
    if os.path.isfile(config.get('rbps_path')):
        rbps = pd.read_csv(os.path.join(config_dirname, config['rbps_path']), index_col=0)
    else:
        rbps = load_rbps()

    rbp_df, isoforms_df = intersect_dfs([rbp_df, isoforms_df])

    gene_data = load_data(gene)
    gene_data['transcripts'] = [t for t in gene_data['transcripts'] if t['transcript_id'] in isoforms_df.columns]
    gene_data = set_variable_exons(gene_data)
    gene_data['sequence'] = gene_data['sequence'].replace('T', 'U')

    output_dir = make_sure_dir_exists(os.path.join(config_dirname, config['output_dir']))
    copyfile(config_path, os.path.join(output_dir, 'config.json'))
    plot_gene_isoforms(gene_data, output_dir=output_dir)

    print('Loaded config...')

    return config, gene_data, rbp_df, isoforms_df, rbps


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
