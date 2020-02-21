"""
Take gene sets and put through neural networks (individual predictions)
to see resultant values
"""

from setup import *
from neuralnet import *
from sparse_x_generator import *

import seaborn as sns
from keras.utils import np_utils
import multiprocessing as mp
from functools import partial


class Perturbation:
    def __init__(self, dir, data_subdir, model_subdir, ae, pretrain, exp_name, label, pos_class_name, pred_type,
                 order_by, y_colname, plot_these_classes, data_type, cpus):
        self.num_cpu = cpus
        self.dir = dir
        self.datadir = data_subdir
        self.model_subdir = model_subdir
        self.ae = ae
        self.pretrain = pretrain
        self.exp_name = exp_name
        self.label = label
        self.pos_class_name = pos_class_name
        self.pred_type = pred_type
        self.order_by = order_by
        self.y_colname = y_colname
        self.plot_these_classes = plot_these_classes
        self.data_type = data_type
        if self.data_type == 'verhaak':
            self.classes = plot_these_classes

        self.exp = self.label + '/neuralnets/' + self.exp_name
        self.model_fn = os.path.join(self.dir, self.model_subdir, self.exp, 'model_retrained.h5')
        self.le = os.path.join(self.dir, self.model_subdir, self.exp, 'labelencoder.pkl')
        if self.pretrain is None:
            self.scaleFN_in = os.path.join(self.dir, self.model_subdir, self.exp, 'geneTrainScalers.pkl')
        else:
            self.scaleFN_in = os.path.join(self.pretrain, 'geneTrainScalers.pkl')

        # prep model input
        self.genes, self.y_labels, self.ids = get_gbm_data(data_dir=self.datadir,
                                                           data_type=self.data_type,
                                                           label_name=self.label)  # read in rg data
        self.gene_ids = list(self.genes.columns.values)
        _, input_dim = self.genes.shape

        scaler = pickle.load(open(self.scaleFN_in, 'rb'))  # scale data using model's training scaler
        genes = scaler.transform(self.genes)
        genes = genes.reshape(len(self.genes), input_dim, 1)  # reshape input into 3-D shape
        for i in self.y_labels.iloc[:, 0].unique():
            print(list(self.y_labels.iloc[:, 0]).count(i), '\t', list(self.y_labels.iloc[:, 0]).count(i) / self.y_labels.shape[0],
                  '\t', i)
        print(self.le)
        if os.path.isfile(self.le):
            self.le = pickle.load(open(self.le, 'rb'))
        self.genes = genes
        self.y_colname = self.y_labels.columns[0]
        self.label_dict = {i: self.y_labels.index.values[self.y_labels[y_colname] == i] for i in self.le.classes_}
        # transform labels
        self.y_labels['y_true'] = self.le.transform(self.y_labels[self.y_colname])  # positive class?

        print('original class names: ', self.le.classes_)  # check classes
        print('transformed class names: ', self.le.transform(self.le.classes_))
        print('order_by label encoding: ', self.le.transform([self.pos_class_name]))  # positive class?

    def get_nn_outputs(self, model, gene_set):
        # get sparse inputs
        num_cases, input_dim, _ = self.genes.shape
        sparse = [ get_sparse_x(genes=gene_set,
                                gene_ids=self.gene_ids,
                                ge_profile=self.genes[i, :, 0]) for i in range(self.genes.shape[0])]
        if sparse[0] is not None:
            sparse = np.reshape(sparse, (num_cases, input_dim, 1))  # reshape input
            output = model.predict(sparse)
        else:
            output = None
        return output

    def get_nn_outputs_all(self, model, gene_set_dict, gs_keys_list=None):
        # use verhaak genes as inputs to nn, and get output
        print('getting outputs')

        if gs_keys_list is None:
            gs_keys_list = list(gene_set_dict.keys())

        perturbations = {}
        for i in gs_keys_list:
            print(i)
            a = self.get_nn_outputs(model, gene_set=gene_set_dict[i])
            if a is None:
                continue
            else:
                a = pd.DataFrame(a)
                a.columns = [ c + ' probability' for c in self.plot_these_classes]  # rename to gene set
                a.index = self.ids
                perturbations[i] = a
        return perturbations

    def order_nn_outputs(self, df):
        # label_dict is a dictionary of true labels, each category in label contains tcga ids
        # e.g., label['PN'] = ids of patients who have PN subtype
        # order by label (e.g., subtype) and true y labels and recombine

        print('available colnames for reordering: ', df.columns.values)
        df_reorder = []
        for i in self.le.classes_:
            print(i)
            sub = df.loc[df.index.intersection(self.label_dict[i])]  # get patients with true label
            # get true y
            labs = self.y_labels.loc[self.y_labels.index.intersection(self.label_dict[i])]
            labs = labs.join(sub)
            labs_reorder = []
            for j in labs[self.y_colname].unique():
                lr = labs.loc[labs[self.y_colname] == j]
                lr = lr.sort_values(by=[self.order_by], ascending=False)
                labs_reorder.extend(lr.index.values)
            labs = labs.loc[labs_reorder]
            df_reorder.append(labs)
        df_reorder = pd.concat(df_reorder, axis=0)

        return df_reorder

    def get_all_plots(self, model,  gs, name, fig_dir):
        # get nn outputs
        mod_pert = self.get_nn_outputs_all(model=model,
                                           gene_set_dict=gs,
                                           gs_keys_list=list(gs.keys()))
        print('plotting outputs')
        t = {}
        for i in list(gs.keys()):
            t[i] = self.order_nn_outputs(df=mod_pert[i])
            g = get_heatmap(values=t[i], gene_set_name=i, row_label_name=self.y_colname)
            g.savefig(fig_dir + '/' + self.label + '_perturbations_' + name + '_' + i + '.png', bbox_inches='tight')
        return mod_pert

    def get_gs_roc(self, preds, y_true_colname, y_pred_colname):
        s = {}
        if self.pred_type == 'binaryClass':
            y_true = preds[[y_true_colname]]
            y_pred = preds[[y_pred_colname]]
            s['fpr'], s['tpr'], _ = roc_curve(y_true=y_true, y_score=y_pred)
            s['roc_auc'] = auc(s['fpr'], s['tpr'])
            s['precision'], s['recall'], _ = precision_recall_curve(y_true=y_true, probas_pred=y_pred)
            s['pr_auc'] = auc(s['recall'], s['precision'])

        if self.pred_type == 'multiClass':
            print('classes order: ', self.le.classes_)
            print('colname order: ', y_true_colname)

            y_true_one_hot = np_utils.to_categorical(preds[[y_true_colname]])  # ints to one-hot encoding
            if self.data_type == 'vasari':
                y_preds_colnames = [self.label + ' ' + i + ' probability' for i in self.le.classes_]
            elif self.data_type == 'verhaak':
                y_preds_colnames = [i + ' probability' for i in self.plot_these_classes]
            else:
                print('no valid data type')

            print(preds.columns.values)
            print(y_preds_colnames)
            y_preds_one_hot = preds[y_preds_colnames]

            s['roc_auc'], s['pr_auc'] = multiclass_metrics(classes=self.le.classes_,
                                                           y_truth=np.array(y_true_one_hot),
                                                           y_preds=np.array(y_preds_one_hot),
                                                           fn=None, t=None, plot=False)
        return s

    def get_scores(self, gs_names, gs_perts, name, tab_dir, save=False):
        # scoring
        if bool(gs_perts):
            pool = mp.Pool(processes=self.num_cpu)
            func = partial(self.get_gs_roc,
                           y_true_colname='y_true',
                           y_pred_colname=self.order_by)
            preds = [gs_perts[i].join(self.y_labels) for i in gs_names]

            # save preds
            if save:
                for i in gs_names:
                    pert_fn = os.path.join(tab_dir, self.label + '_' + name + '_' + i + '_preds.csv')
                    print('predictions saved: ', pert_fn)
                    gs_perts[i].join(self.y_labels).to_csv(pert_fn)
            scores = pool.map(func, preds)
            scores = dict(zip(gs_names, scores))

            if self.pred_type == 'binaryClass':
                scores = pd.DataFrame.from_dict(scores, orient='index')
                scores = scores.sort_values(by=['roc_auc'], ascending=False)
                scores.to_csv(os.path.join(tab_dir, self.label + '_' + name + '_scores.csv'))
                scores.sort_values(by=['roc_auc'], ascending=False)
            else:
                s = {}
                s['roc_auc'] = pd.DataFrame.from_dict({(i, 'roc_auc'): scores[i]['roc_auc']  # contains f1_scores
                                                       for i in scores.keys()}, orient='index')
                s['pr_auc'] = pd.DataFrame.from_dict({(i, 'pr_auc'): scores[i]['pr_auc']
                                                      for i in scores.keys()}, orient='index')

                s['roc_auc'].to_csv(os.path.join(tab_dir, self.label + '_' + name + '_roc_scores.csv'))
                s['pr_auc'].to_csv(os.path.join(tab_dir, self.label + '_' + name + '_pr_scores.csv'))
                scores = s
        else:
            scores = None

        return scores


def get_puchalski_gene_sets(gs_file):
    puch = pd.read_csv(gs_file, sep=',', header=0).drop(0, axis='index')
    puch_dict = {i: puch[i].dropna().values for i in puch.columns}
    return puch_dict


def get_verhaak_gene_sets(gs_file):
    # get verhaak gene set
    verhaak = pd.read_table(gs_file, sep='\t', header=1)
    verhaak = verhaak[['Subtype', 'GO genes']]
    verhaak.columns = ['go_genes', 'subtype']
    all_genes = verhaak['go_genes'].values
    subtypes = ['NL', 'PN', 'CL', 'MES']
    verhaak = {i: verhaak['go_genes'].loc[verhaak['subtype'] == i].values for i in subtypes}
    verhaak['all'] = all_genes
    return verhaak


def get_bi_gene_sets(gs_file):
    # read msigdb gmt file
    with open(os.path.join(gs_file), "r") as f:
        gs = {}
        for line in f:
            l = line.rstrip('\n')
            l = l.split('\t')
            gs[l[0]] = l[2:]
    names = list(gs.keys())
    sizes = np.array([len(gs[i]) for i in names])
    print(gs_file)
    print(len(sizes))
    print(sum(sizes))
    print(min(sizes))
    print(max(sizes))
    return gs, names, sizes


def get_random_gs(genes, gs_exclude, num=200):
    # random gene set
    # gs_exclude - genes to exclude because they're being looked at e.g., verhaak_gs['all']
    np.random.seed(504)
    rand_pool = [g for g in genes.columns.values if g not in gs_exclude]
    random_gs = {'random': np.random.choice(rand_pool, size=num)}
    return random_gs


def get_all_gs(datadir):
    # gene set files

    # mig sig gs
    gs_dir = datadir + '/msigdb_v6.2_GMTs'
    gs_files = [i for i in os.listdir(gs_dir) if i.endswith('symbols.gmt')]

    gs = {}
    gs_names = {}
    gs_sizes = {}
    for i in gs_files:
        k = re.sub('.v6.2.symbols.gmt', '', i)
        gs[k], gs_names[k], gs_sizes[k] = get_bi_gene_sets(os.path.join(gs_dir, i))

    # verhaak gene sets
    verhaak_file = datadir + '/TCGA_unified_CORE_ClaNC840.txt'
    verhaak_gs = get_verhaak_gene_sets(verhaak_file)

    # also need verhaak true labels
    # create dictionary of subtypes and pat ids
    genes, verhaak_labels, _ = get_gbm_data(data_dir=os.path.dirname(verhaak_file), data_type='verhaak',
                                            label_name='subtype')
    verhaak_labels = {k: v.index.values for k, v in verhaak_labels.groupby('subtype')}
    verhaak_labels['NL'] = verhaak_labels.pop('Neural')  # consistent naming
    verhaak_labels['PN'] = verhaak_labels.pop('Proneural')
    verhaak_labels['CL'] = verhaak_labels.pop('Classical')
    verhaak_labels['MES'] = verhaak_labels.pop('Mesenchymal')

    gs['verhaak_gs'] = verhaak_gs
    gs_names['verhaak_gs'] = list(verhaak_gs.keys())
    gs_sizes['verhaak_gs'] = np.array([len(verhaak_gs[i]) for i in verhaak_gs.keys()])

    # random gs
    rand = {'random_100': get_random_gs(genes=genes, gs_exclude=verhaak_gs['all'], num=100)['random'],
            'random_200': get_random_gs(genes=genes, gs_exclude=verhaak_gs['all'], num=200)['random']}

    gs['verhaak_rand_gs'] = rand
    gs_names['verhaak_rand_gs'] = ['random']
    gs_sizes['verhaak_rand_gs'] = np.array([100, 200])

    gs['test_gs'] = {'test': ['na', 'notagene', ' blah'],
                     'random': gs['verhaak_rand_gs']['random_100']}
    gs_names['test_gs'] = ['test', 'random']
    gs_sizes['test_gs'] = np.array([3, 100])

    # puchalski gene sets used
    puchalski_gs = get_puchalski_gene_sets(datadir + '/gene_sets_Puchalski/aaf2666_Table-S15.csv')
    gs['puchalski_gs'] = puchalski_gs
    gs_names['puchalski_gs'] = list(puchalski_gs.keys())
    gs_sizes['puchalski_gs'] = np.array([len(puchalski_gs[i]) for i in puchalski_gs.keys()])

    return gs, gs_names, gs_sizes


def get_heatmap(values, gene_set_name, row_label_name):
    s = values.pop(row_label_name)
    lut = dict(zip(s.unique(), "rbgy"))
    row_colors = s.map(lut)

    g = sns.clustermap(values, row_colors=row_colors, cmap='viridis',
                       row_cluster=False, col_cluster=False)
    g.fig.suptitle(gene_set_name + ' gene set perturbation')
    # Draw the legend bar for the classes
    for label in s.unique():
        g.ax_col_dendrogram.bar(0, 0, color=lut[label],
                                label=label, linewidth=0)
    g.ax_col_dendrogram.legend(loc="center left", ncol=1)
    return g


def get_all_perts(Pertrubation, model, gs, gs_names, gs_type,
                  tab_dir, fig_dir=None, save_preds=False):
    # todo : use get_all_plots instead when fig_dir is not none
    perts = {}
    scores = {}

    if gs_type == 'single_gene':
        single_gs = {i: np.array([i]) for i in Pertrubation.gene_ids}  # create gene dict with single genes

        if gs_names is None:
            gs_names = list(single_gs.keys())
        perts['single'] = Pertrubation.get_nn_outputs_all(model=model,
                                                          gene_set_dict=single_gs,
                                                          gs_keys_list=gs_names)  # not plotting
        scores['single'] = Pertrubation.get_scores(gs_names=gs_names,
                                                   gs_perts=perts['single'],
                                                   name='single',
                                                   tab_dir=tab_dir,
                                                   save=save_preds)
    else:
        gs_info = {'verhaak': 'verhaak_gs',
                   'random_verh': 'verhaak_rand_gs',
                   'puchalski': 'puchalski_gs',
                   'hallmark': 'h.all',
                   'chromosome': 'c1.all',
                   'reactome': 'c2.cp.reactome',
                   'biocarta': 'c2.cp.biocarta',
                   'kegg': 'c2.cp.kegg',
                   'canonical_paths': 'c2.cp',
                   'chem_gene_perts': 'c2.cgp',
                   'motif': 'c3.all',
                   'computational': 'c4.all',
                   'GO': 'c5.all',
                   'onco_sig': 'c6.all',
                   'immuno_sig': 'c7.all',
                   'test': 'test_gs'
                   }
        if fig_dir is None:
            k = gs_type
            v = gs_info[k]

            if gs_names is None:
                gs_names = list(gs[v].keys())
            # gs_names = ['chr2p14']  # testing
            perts = {k: Pertrubation.get_nn_outputs_all(model=model,
                                                        gene_set_dict=gs[v],
                                                        gs_keys_list=gs_names)}

            gs_names = list(perts[k].keys())
            scores = {k:  Pertrubation.get_scores(gs_names=gs_names,
                                                  gs_perts=perts[k],
                                                  name=k,
                                                  tab_dir=tab_dir,
                                                  save=save_preds)}
    return perts, scores
