import os
import sys
import csv
import logging

import pandas as pd
from sklearn.tree import DecisionTreeClassifier

from ml_utils import load_polya_seq_df, map_clvs, compare, KARBOR_FEATURE_COLS
from cluster import cluster_clv_sites
from kleat.misc.utils import backup_file


logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')


def load_df(sample_id):
    infile = os.path.abspath(
        './benchmark_transcriptome/{0}.ml_ready.pkl'.format(sample_id))
    logging.info('reading {0} ...'.format(infile))
    return pd.read_pickle(infile)


def map_df(df, sample_id, mapping_cutoff=50):
    df_ref = load_polya_seq_df(sample_id)

    df = map_clvs(df, df_ref)
    df['is_tp'] = df.abs_dist < mapping_cutoff
    return df


def train(df_tr_mapped, max_depth):
    Xs = df_tr_mapped[KARBOR_FEATURE_COLS]
    ys = df_tr_mapped.is_tp.values
    clf = DecisionTreeClassifier(max_depth=max_depth)
    logging.info('training in process ...')
    clf.fit(Xs, ys)
    return clf


def predict(test_sample_id, clf):
    df_te = load_df(test_sample_id)
    df_te_mapped = map_df(df_te, test_sample_id)

    Xs_te = df_te_mapped[KARBOR_FEATURE_COLS]
    df_te_mapped['predicted'] = clf.predict(Xs_te)

    clv_key_cols = ['seqname', 'strand', 'clv']
    out_df = df_te_mapped.query('predicted')[clv_key_cols].drop_duplicates()
    return out_df


def cluster_clv(df, cutoff=20):
    grped = df.groupby(['seqname', 'strand'])
    applied = grped.apply(cluster_clv_sites, cutoff)
    dedupped = applied[['seqname', 'strand', 'mclv']].drop_duplicates()
    out = dedupped.rename(columns={'mclv': 'clv'}).reset_index(drop=True)
    return out


MAP_CUTOFF = 50
SAMPLE_IDS = [
    'HBRC4',
    'HBRC6',
    'UHRC1',
    'UHRC2'
]

ref_sample_id = sys.argv[1]
df_tr = load_df(ref_sample_id)
df_tr_mapped = map_df(df_tr, ref_sample_id)

max_depth_list = range(2, 10, 1)
outfile = './out.csv'
backup_file(outfile)
with open(outfile, 'wt') as opf:
    csvwriter = csv.writer(opf)
    csvwriter.writerow(
        ['sample_id', 'precision', 'recall', 'f1', 'tree_max_depth']
    )
    for max_depth in max_depth_list:
        print(f'training on {ref_sample_id} with max tree depth: {max_depth}')
        clf = train(df_tr_mapped, max_depth)

        for test_sample_id in SAMPLE_IDS:
            print(f'testing on {test_sample_id}')
            df_predicted = predict(test_sample_id, clf)
            df_clustered = cluster_clv(df_predicted)
            df_ref = load_polya_seq_df(test_sample_id)

            reca, prec, f1 = compare(df_clustered, df_ref)
            csvwriter.writerow(
                [test_sample_id, prec, reca, f1, max_depth]
            )
