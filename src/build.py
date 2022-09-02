import pickle
import sys

from src.helpers.plots import plot_isoforms_tree
from src.lr import elastic_net
from src.helpers.pipeline import map_motifs_to_exons, make_exons_sf_df, load_config_and_input_data


def main(config_path):
    # Load config and input data
    config, gene_data, rbp_df, isoforms_df, rbps = load_config_and_input_data(config_path)
    print(rbps)

    exons_motifs = map_motifs_to_exons(gene_data, rbps)
    tree = make_exons_sf_df(
        gene_data,
        rbp_df, isoforms_df,
        gene_exon_motifs=exons_motifs,
    )
    plot_isoforms_tree(tree, config['output_dir'])

    nodes = [tree.left_child, tree.right_child]
    while len(nodes):
        for node in nodes[::2]:
            df = node.df
            if len(df.columns) > 2:
                res = elastic_net(df)
                tree.res = res
                if config['tissue_specific'] and res is not None:
                    tree.tissue_res = {}
                    for tissue in set(df['Tissue']):
                        tissue_df = df[df['Tissue'] == tissue]
                        tissue_df['Freq'] = 1
                        res = elastic_net(tissue_df)
                        tree.tissue_res[tissue] = res

        cur_nodes = []
        for node in nodes:
            if node.left_child is not None:
                cur_nodes += [node.left_child, node.right_child]
        nodes = cur_nodes

    with open(f'{config["output_dir"]}/tree.wb', 'wb') as res_file:
        pickle.dump(tree, res_file)

    return tree


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Please specify configuration file', file=sys.stderr)
        sys.exit(1)

    main(sys.argv[1])
