import argparse
import copy
import operator
import re

import numpy as np
import pandas as pd
import scanpy as sc
from SCCAF import *
from sklearn.metrics import f1_score

## functions

def get_top_features(adata: sc.AnnData, clusters=list, topn=3) -> pd.DataFrame:
    np.random.seed(42)
    y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF_assessment(adata.X, adata.obs[clusters])
    return get_topmarkers(clf, list(adata.var.index), topn=topn).set_index('gene')


def get_plot(ad: sc.AnnData, clusters: str, title: str, out_file: str) -> None:
    plot = sc.pl.umap(ad, color=[clusters], title=title, show=False)
    fig = plot.get_figure()
    fig.savefig(out_file)


def prepare_data(path: str, normalize_sc=False) -> sc.AnnData:
    adata = sc.read(path)
    adata.obs['louvain'] = adata.obs['seurat_clusters']
    if normalize_sc:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    sc.pp.neighbors(adata)
    return adata


def optimize_clusters(adata: sc.AnnData, res: int, k_added: str, prefix: str) -> sc.AnnData:
    np.random.seed(42)
    sc.tl.louvain(adata, resolution=res, key_added=k_added)
    SCCAF_optimize_all(ad=adata, plot=False, min_acc=0.9, prefix=prefix, use='pca')
    return adata


def loop_opt(ad, min_r, max_r, step):
    ad = copy.deepcopy(ad)
    opt_res = dict()
    for res in np.arange(min_r, max_r, step):
        np.random.seed(42)
        try:
            sc.tl.louvain(ad, resolution=np.round(res, 1), key_added=f'L{np.round(res, 1)}_Round0')
            y_prob, y_pred, y_test, clf, cvsm, acc = SCCAF_assessment(ad.X, ad.obs[ad.obs.keys()[-1]])
            opt_res[np.round(res, 1)] = f1_score(y_pred, y_test, average='micro')
        except:
            continue
    return opt_res


def find_opt_res(ad, min_r=1, max_r=3, step=1, thresh=0.9):
    opt_res_fir = max(loop_opt(ad, min_r, max_r, step).items(), key=operator.itemgetter(1))[0]
    opt_res_dict = loop_opt(ad, min_r=opt_res_fir - 0.1, max_r=opt_res_fir, step=0.2)
    try:
        opt_res_sec = max(k for k, v in opt_res_dict.items() if v > thresh)
    except:
        opt_res_sec = max(opt_res_dict.items(), key=operator.itemgetter(1))[0]
    return \
        max(loop_opt(ad, min_r=opt_res_sec - 0.05, max_r=opt_res_sec + 1, step=0.05).items(),
            key=operator.itemgetter(1))[0]

def pipeline(h5: str, sample_id: str, plt_before: str, plt_after: str, log: str, markers: str) -> None:
    ## load data
    adata = prepare_data(h5)
    ## draw plot with seurat clusters as labels
    get_plot(ad=adata, clusters='seurat_clusters', title=f'{sample_id}: seurat',
             out_file=plt_before)
    ## find optimal resolution based on F1 metrics
    opt_res = find_opt_res(adata)
    ## optimize clusters with given optimal resolution
    adata = optimize_clusters(adata, res=opt_res, k_added=f'L{opt_res}_Round0', prefix=f'L{opt_res}')
    ## Get name of all rounds, write some output to log file
    rounds_num = list(filter(lambda x: re.findall(f'L{opt_res}_Round.*self*', x), list(adata.obs.keys())))
    with open(log, 'a+') as log_file:
        log_file.write(f'Optimal resolution = {opt_res}\n')
        log_file.write(f'Number of rounds = {len(rounds_num)}\n')
    ## draw plot with optimized clusters as labels
    get_plot(ad=adata, clusters=adata.obs.keys()[-1],
             title=f'{sample_id}: sccaf, optimal resolution: {opt_res}, rounds: {len(rounds_num)}', out_file=plt_after)
    ## write markers
    try:
        get_top_features(adata, adata.obs.keys()[-1], topn=10).to_csv(markers, index=True, sep='\t')
    except:
        with open(markers, 'w') as mrs_file:
            mrs_file.write("This solver needs samples of at least 2 classes in the data, but the data contains only one class: '0'")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Define 10x version")
    parser.add_argument('--h5', type=str, required=True,
                        help='h5 path')
    parser.add_argument('--sample_id', type=str, required=True,
                        help='Number of threads for parallel-fastq-dump program')
    parser.add_argument('--plt_before', type=str, required=True,
                        help='Filename: plot with seurat labels')
    parser.add_argument('--plt_after', type=str, required=True,
                        help='Filename: plot with sccaf labels')
    parser.add_argument('--log', type=str, required=True,
                        help='Filename for log file')
    parser.add_argument('--markers', type=str, required=True,
                        help='Filename for markers file')
    args = parser.parse_args()
    pipeline(args.h5, args.sample_id, args.plt_before, args.plt_after, args.log, args.markers)
