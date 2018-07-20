"""
check the maximum recall of predicted clvs without filtering, but clustered
"""

import sys
import logging

import pandas as pd

from kleat.post import cluster_clv_parallel
from train_arbor import map_to_ref, load_polya_seq_df


logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')


if __name__ == "__main__":
    sample_id = sys.argv[1]
    infile = f'./benchmark_transcriptome/tcga-run-3/{sample_id}.pkl'
    print(f'reading {infile}')
    adf = pd.read_pickle(infile)
    print(f'before clustering df.shape: {adf.shape}')

    print('clustering ...')
    adf_clustered = cluster_clv_parallel(adf)
    logging.info('dropping duplicates after clustering ...')

    bdf = adf_clustered[['seqname', 'strand', 'mode_clv']].drop_duplicates()
    bdf = bdf.rename(columns={'mode_clv': 'clv'}).reset_index(drop=True)
    print(f'after clustering df.shape: {bdf.shape}')

    print('mapping ...')
    map_cutoff = 50
    df_ref = load_polya_seq_df(sample_id)
    df_mapped = map_to_ref(bdf, df_ref, map_cutoff)
    df_mapped['is_tp'] = df_mapped.abs_dist < map_cutoff

    print('calculating metrics ...')
    recall = df_mapped.query('is_tp').shape[0] / df_ref.shape[0]
    precision = df_mapped.query('is_tp').shape[0] / df_mapped.shape[0]
    f1 = (2 * recall * precision) / (recall + precision)

    outfile = f'max_recall.{sample_id}'
    print(f'writing results to {outfile} ...')
    with open(outfile, 'wt') as opf:
        opf.write('sample_id\tprecision\trecall\tf1\n')
        opf.write(f'{sample_id}\t{precision}\t{recall}\t{f1}\n')
