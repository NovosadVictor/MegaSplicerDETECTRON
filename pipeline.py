from multiprocessing import Pool
import pandas as pd
import numpy as np
import traceback

from helpers.pipeline import \
    load_data, \
    map_motifs_to_exons, \
    make_exons_sf_df, \
    load_isoforms, \
    set_variable_exons
from consts import genes_data, tcga_sfs
from lr import elastic_net
from utils import make_sure_dir_exists


def process(gene):
    try:
        make_sure_dir_exists(f'results/{gene}')
        gene_data = load_data(gene)

        gene_isoforms = load_isoforms(gene)
        gene_data['transcripts'] = [t for t in gene_data['transcripts'] if t['transcript_id'] in gene_isoforms.columns]

        gene_data = set_variable_exons(gene_data)

        gene_data['sequence'] = gene_data['sequence'].replace('T', 'U')

        if 0 < len(gene_data['variable_exons']) <= 20:
            exons_motifs = map_motifs_to_exons(gene_data)
            tree = make_exons_sf_df(
                gene_data,
                tcga_sfs, gene_isoforms,
                gene_exon_motifs=exons_motifs,
            )
            nodes = [tree.left_child, tree.rightt_child]
            while len(nodes):
                for node in nodes[::2]:
                    df = node.df
                    if len(df.columns) > 2:
                        res = elastic_net(df)
                        tree.res = res

                cur_nodes = []
                for node in nodes:
                    if node.left_child is not None:
                        cur_nodes += [node.left_child, node.right_child]
                nodes = cur_nodes

            return tree
        return None
    except Exception:
        print('Exception for ', gene)
        traceback.print_exc()

        return None


def process_chunk(gene_names_chunk, process_ind):
    results = pd.DataFrame()
    for i, gene in enumerate(gene_names_chunk):
        res = process(gene)
        if res is not None and not res.empty:
            results = pd.concat([results, res], axis=0)

        if i % 100 == 0:
            results.to_csv(f'{process_ind}_coefs_results.csv')

    return results


n_processes = 50
gene_names = list(genes_data.keys())
# gene_names = ['CD44']
gene_names_chunks = np.array_split(gene_names, n_processes)
with Pool(n_processes) as p:
    results = p.starmap(
        process_chunk,
        [(gene_names_chunks[process_ind], process_ind) for process_ind in range(n_processes)],
    )
    results = pd.concat(results, axis=0)
    results.to_csv('coefs_results.csv')
    print(results)
