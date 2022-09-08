import pickle
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split

from src.helpers.pipeline import map_motifs_to_exons, make_exons_sf_df
from src.helpers.plots import plot_isoforms_tree, plot_gene_isoforms
from src.lr import elastic_net
from src.utils.common import predict, get_accuracy, get_scores, add_freq_to_df


class Pipeline:
    def __init__(self, config, gene_data, rbp_df, isoforms_df, rbps):
        self.config = config
        self.gene_data = gene_data
        self.rbp_df = rbp_df
        self.isoforms_df = isoforms_df
        self.rbps = rbps
        self.tree = None
        self.tissue_specific = self.config.get('tissue_specific', False) and 'Tissue' in self.rbp_df
        if self.tissue_specific:
            self.rbp_df = add_freq_to_df(self.rbp_df)

        if 'Dataset.Type' in self.rbp_df.columns:
            self.train_index = self.rbp_df[self.rbp_df['Dataset.Type'] == 'Training'].index
            self.val_index = self.rbp_df[self.rbp_df['Dataset.Type'] == 'Validation'].index
        else:
            stratify = self.rbp_df['Tissue'] if self.tissue_specific else None
            self.train_index, self.val_index = train_test_split(self.rbp_df.index, test_size=.25, stratify=stratify)

    def run(self):
        exons_motifs = map_motifs_to_exons(self.gene_data, self.rbps)
        tree = make_exons_sf_df(
            self.gene_data,
            self.rbp_df, self.isoforms_df,
            gene_exon_motifs=exons_motifs,
        )

        nodes = [tree.left_child, tree.right_child]
        while len(nodes):
            for node in nodes[::2]:
                df = node.df.loc[self.train_index]
                if len(df.columns) > 2:
                    res = elastic_net(df)
                    node.res = res
                    if self.tissue_specific:
                        node.tissue_res = {}
                        for tissue in set(df['Tissue']):
                            print(tissue)
                            tissue_df = df[df['Tissue'] == tissue]
                            tissue_df['Freq'] = 1
                            res = elastic_net(tissue_df)
                            node.tissue_res[tissue] = res
            cur_nodes = []
            for node in nodes:
                if node.left_child is not None:
                    cur_nodes += [node.left_child, node.right_child]
            nodes = cur_nodes

        with open(f'{self.config["output_dir"]}/pipeline.wb', 'wb') as res_file:
            pickle.dump(self, res_file)

        plot_gene_isoforms(self.gene_data, output_dir=self.config['output_dir'])
        plot_isoforms_tree(tree, output_dir=self.config['output_dir'])

        self.tree = tree
        self.predict()
        self.accuracy()
        self.plot()

    @staticmethod
    def load_from_file(path_to_file):
        with open(path_to_file) as class_file:
            return pickle.load(class_file)

    def predict(self):
        if self.tree is None:
            return

        parents = [self.tree.left_child, self.tree.right_child]
        while len(parents):
            cur_parents = []
            for node_id, parent in enumerate(parents):
                if node_id % 2 == 0:
                    parent.res['predictions.train'] = predict(parent.df.loc[self.train_index], parent.res['coefs'])
                    parent.res['predictions.validation'] = predict(parent.df.loc[self.val_index], parent.res['coefs'])

                    for tissue in getattr(parent, 'tissue_res', {}):
                        parent.tissue_res[tissue]['predictions.train'] = predict(
                            parent.df.loc[self.train_index].query(f'Tissue == "{tissue}"'),
                            parent.tissue_res[tissue]['coefs'],
                        )
                        parent.tissue_res[tissue]['predictions.validation'] = predict(
                            parent.df.loc[self.val_index].query(f'Tissue == "{tissue}"'),
                            parent.tissue_res[tissue]['coefs'],
                        )
                else:
                    parent.res['predictions.train'] = 1 - parents[node_id - 1].res['predictions.train']
                    parent.res['predictions.validation'] = 1 - parents[node_id - 1].res['predictions.validation']

                    for tissue in getattr(parent, 'tissue_res', {}):
                        parent.tissue_res[tissue]['predictions.train'] = 1 - parents[node_id - 1].tissue_res[tissue]['predictions.train']
                        parent.tissue_res[tissue]['predictions.validation'] = 1 - parents[node_id - 1].tissue_res[tissue]['predictions.validation']

                parent.res['predictions.train.accumulative'] = parent.parent.res.get('predictions.train.accumulative', 1) * parent.res['predictions.train']
                parent.res['predictions.validation.accumulative'] = parent.parent.res.get('predictions.validation.accumulative', 1) * parent.res['predictions.validation']

                for tissue in getattr(parent, 'tissue_res', {}):
                    parent.tissue_res[tissue]['predictions.train.accumulative'] = parent.parent.tissue_res.get(tissue, {}).get('predictions.train.accumulative', 1) * parent.tissue_res[tissue]['predictions.train']
                    parent.tissue_res[tissue]['predictions.validation.accumulative'] = parent.parent.tissue_res.get(tissue, {}).get('predictions.validation.accumulative', 1) * parent.tissue_res[tissue]['predictions.validation']

                if parent.left_child is not None:
                    cur_parents += [parent.left_child, parent.right_child]

            parents = cur_parents

    def accuracy(self):
        if self.tree is None:
            return

        parents = [self.tree.left_child, self.tree.right_child]
        while len(parents):
            cur_parents = []
            for node_id, parent in enumerate(parents):
                print(parent.res['predictions.train'], parent.df.loc[self.train_index]['fraq'])
                train_acc = get_scores(parent.res['predictions.train'], parent.df.loc[self.train_index]['fraq'])
                val_acc = get_scores(parent.res['predictions.validation'], parent.df.loc[self.val_index]['fraq'])
                parent.res['accuracy'] = {
                    'train': train_acc,
                    'validation': val_acc,
                }
                for tissue in getattr(parent, 'tissue_res', {}):
                    train_acc = get_scores(parent.tissue_res[tissue]['predictions.train'], parent.df.loc[self.train_index].query(f'Tissue == "{tissue}"')['fraq'])
                    val_acc = get_scores(parent.tissue_res[tissue]['predictions.validation'], parent.df.loc[self.val_index].query(f'Tissue == "{tissue}"')['fraq'])
                    parent.tissue_res[tissue]['accuracy'] = {
                        'train': train_acc,
                        'validation': val_acc,
                    }
                if parent.left_child is not None:
                    cur_parents += [parent.left_child, parent.right_child]
                else: # leaf node (isoform)
                    parent.res['accuracy']['train.accumulative'] = get_scores(
                        parent.res['predictions.train.accumulative'],
                        parent.df.loc[self.train_index]['fraq'],
                    )
                    parent.res['accuracy']['validation.accumulative'] = get_scores(
                        parent.res['predictions.validation.accumulative'],
                        parent.df.loc[self.val_index]['fraq'],
                    )

                    for tissue in getattr(parent, 'tissue_res', {}):
                        parent.tissue_res[tissue]['accuracy']['train.accumulative'] = get_scores(
                            parent.tissue_res[tissue]['predictions.train.accumulative'],
                            parent.df.loc[self.train_index].query(f'Tissue == "{tissue}"')['fraq'],
                        )
                        parent.tissue_res[tissue]['accuracy']['validation.accumulative'] = get_scores(
                            parent.tissue_res[tissue]['predictions.validation.accumulative'],
                            parent.df.loc[self.val_index].query(f'Tissue == "{tissue}"')['fraq'],
                        )

            parents = cur_parents

    def plot(self):
        if self.tree is None:
            return

        parents = [self.tree.left_child, self.tree.right_child]
        isoform_accuracies = {}
        while len(parents):
            cur_parents = []
            for node_id, parent in enumerate(parents):
                if parent.left_child is not None:
                    cur_parents += [parent.left_child, parent.right_child]
                else:
                    isoform_accuracies[parent.kwargs[0]['transcript_id']] = {
                        'var.train': parent.df.loc[self.train_index]['fraq'],
                        'var.validation': parent.df.loc[self.val_index]['fraq'],
                        'train': parent.res['accuracy']['train.accumulative'],
                        'validation': parent.res['accuracy']['validation.accumulative'],
                        'tissue': {
                            tissue: {
                                'train': parent.tissue_res[tissue]['accuracy']['train.accumulative'],
                                'validation': parent.tissue_res[tissue]['accuracy']['validation.accumulative'],
                                'var.train': parent.df.loc[self.train_index].query(f'Tissue == "{tissue}"')['fraq'],
                                'var.validation': parent.df.loc[self.val_index].query(f'Tissue == "{tissue}"')['fraq'],
                            } for tissue in parent.tissue_res
                        }
                    }
            parents = cur_parents

        for score in ['cor', 'uplift.score']:
            plt.figure(figsize=(8, 6))
            ax = sns.scatterplot(
                x=isoform_accuracies.keys(), y=[isoform_accuracies[iso]['validation'][score] for iso in isoform_accuracies],
                color='blue', s=200,
            )
            ax.legend(['Combined', 'Trained'])
            plt.xticks(rotation=90)
            ax.set_ylim(0, 1)
            plt.grid(visible=True)
            plt.tight_layout()
            plt.savefig(f"{self.config['output_dir']}/{score}.png", dpi=300)
