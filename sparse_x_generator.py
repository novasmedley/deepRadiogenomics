"""
Functions to create sparse gene expression profiles, i.e., gene masking.
"""

import numpy as np
import pandas as pd
import re
import os
from random import shuffle
import pickle


def get_gene_ids(gene_names, symbols=True):
    if symbols:
        gene_ids = list(map(str, gene_names))
    else:
        gene_ids = [i.replace('geneid.', '') for i in gene_names]
        gene_ids = list(map(int, gene_ids))

    return gene_ids


def get_broadinst_gene_sets(gmt_file):
    gs = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            genes = line.strip('\n').split('\t')
            gs[genes[0]] = genes[2:]
    gs_keys = list(gs.keys())
    sizes = np.array([len(gs[i]) for i in gs_keys])
    return gs, gs_keys, sizes


def get_subset_gene_set(array_ids, gmt_file, min_size=100):
    gs, keys, sizes = get_broadinst_gene_sets(gmt_file=gmt_file)

    # get gs with at least 100 genes
    keys = np.array(keys)
    keys = keys[np.where(sizes >= min_size)]
    print('gene sets num: ', len(keys))
    gs = {k: gs.get(k) for k in keys}

    # get genes set in gene set
    genes_gs = {k: gs.get(k) for k in keys}  # genes in each list
    genes_gs = [j for i in genes_gs.values() for j in i]  # unlist
    gene_set = set(genes_gs)

    # get genes that are not in the gene set
    t = {k: 1 for k in array_ids if k not in gene_set}
    print('array genes not in gene set (added as a gene set): ', len(t))
    gs['others'] = list(t.keys())  # add other genes to gene set

    return gs


def get_broadinst_collection(gmt_files):
    gs = {}
    gs_names = {}
    gs_sizes = {}
    for i in gmt_files:
        k = re.sub('.v6.2.symbols.gmt', '', os.path.basename(i))
        gs[k], gs_names[k], gs_sizes[k] = get_broadinst_gene_sets(i)
    return gs, gs_names, gs_sizes


def format_chrome(gs, gs_keys):
    # further format chrome so each geneset is one chrome
    k = ['chr' + str(i) for i in np.arange(1, 23)]
    k.append('chry')
    k.append('chrx')
    gs = {t: sum([gs[i] for i in gs_keys if re.match(t + 'q', i)], []) for t in k}
    sizes = np.array([len(gs[i]) for i in gs_keys])
    return gs, k, sizes


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


def print_sizes(sizes):
    print(len(sizes))
    print(sum(sizes))
    print(min(sizes))
    print(max(sizes))


def get_gs_subset(gs, gs_sizes, min_genes=None):
    # remove gs that do not have at least x genes
    k = np.array(list(gs.keys()))
    k = k[gs_sizes >= min_genes]
    print(len(k))
    sub_gs = {i: gs.get(i) for i in k}
    return sub_gs


def get_sparse_idx(genes, gene_ids):
    # genes = list of genes in pathway
    # gene_ids = list of gene_ids of all genes
    gene_idx = [gene_ids.index(g) for g in genes if g in gene_ids]  # get the index of genes in pathway
    gene_idx = np.array(gene_idx)
    return gene_idx


def get_sparse_genes(gene_idx, ge_profile):
    gene_dim = ge_profile.shape[0]
    ge_profile = np.reshape(ge_profile, (gene_dim,))
    sparse_x = np.zeros(gene_dim)  # create empty vector
    sparse_x[gene_idx] = ge_profile[gene_idx]  # set values to patient values for each gene in pathway
    return sparse_x


def get_sparse_x(genes, gene_ids, ge_profile):
    # create sparse input to nn based on patient's values for pathway
    # genes = list of genes in pathway
    # gene_ids = list of gene_ids of all genes
    # ge_profile = array of gene expressions for one person
    gene_idx = get_sparse_idx(genes, gene_ids)  # get the index of genes in pathway
    if gene_idx.size !=0 :
        sparse_x = get_sparse_genes(gene_idx=gene_idx, ge_profile=ge_profile)
    else:
        sparse_x = None
    return sparse_x


def generate_data(y, gene_expr, gene_sets, gene_ids, batch, gs_idx_fn, sparse_batch=10, ae=False):
    # gene_expr = data frame, rows of patients, cols of gene expressions
    # batch = batch size
    # pathways = dictionary pathway key is name and values are list of genes

    while True:
        if ae:
            n, gene_dim = gene_expr.shape
        else:
            n, gene_dim, _ = gene_expr.shape

        path_names = list(gene_sets.keys())
        shuffle(path_names)

        # get sparse indices beforehand
        if os.path.isfile(gs_idx_fn):
            gs_indx = pickle.load(open(gs_idx_fn, 'rb'))
        else:
            gs_indx = {i: get_sparse_idx(genes=gene_sets.get(i), gene_ids=gene_ids) for i in path_names}
            pickle.dump(gs_indx,  open(gs_idx_fn, 'wb'))

        for i in range(0, n, batch):  # iterate through people
            iend = min(i+batch, n)

            for j in range(0, len(path_names), sparse_batch):  # iterate through gene set
                jend = min(j+sparse_batch, len(path_names))
                reps = jend - j

                if ae:
                    # inputs
                    paths = [get_sparse_genes(gene_idx=gs_indx.get(path_names[path]),
                                              ge_profile=gene_expr[person, :])
                             for path in range(j, jend)
                             for person in range(i, iend)]  # sparse paths
                    paths = np.array(paths)

                    # outputs
                    labels = np.tile(gene_expr[i:iend, :], (reps, 1))
                else:
                    # inputs
                    paths = [get_sparse_genes(gene_idx=gs_indx.get(path_names[path]),
                                              ge_profile=gene_expr[person, :, 0])
                             for path in range(j, jend)
                             for person in range(i, iend)]  # sparse paths
                    paths = np.array(paths)
                    paths = paths.reshape(len(paths), gene_dim, 1)  # reshape input into 3-D shape

                    # outputs  TODO check for group predictions, y is array
                    labels = y[i:iend]
                    labels = np.tile(labels, reps)

                yield (paths, labels)
