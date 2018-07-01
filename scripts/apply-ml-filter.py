import os
import sys
import pickle

import numpy as np
import pandas as pd

import pysam

from utils.cluster import cluster_clv_sites
from utils import search_hexamer

kleat_version = sys.argv[1]
train_data_sample_id = sys.argv[2]
SAMPLE_ID = sys.argv[3]


if kleat_version == 'kleat2':
    from utils.utils import KLEAT2_FEATURE_COLS as feature_cols
elif kleat_version == 'kleat3':
    from utils.utils import KLEAT3_FEATURE_COLS as feature_cols
else:
    raise

print(f'working on {kleat_version} {SAMPLE_ID}')

depths = range(2, 20, 1)
for max_depth in depths:
    print(f'working on depth: {max_depth}')
    outdir = f'./{kleat_version}_ml/pkl/DT/{train_data_sample_id}/max_depth{max_depth}'
    clf_pkl = f'{outdir}/clf.pkl'

    fdf = pd.read_csv(
        f'./{kleat_version}_ml/{SAMPLE_ID}_ml_ready.csv', low_memory=False)

    with open(clf_pkl, 'rb') as inf:
        clf = pickle.load(inf)

    fdf['predicted'] = clf.predict(fdf[feature_cols])

    pre_gdf = fdf.query('predicted')[['seqname', 'strand', 'clv']]

    CUTOFF = 20
    gdf = pre_gdf.groupby(['seqname', 'strand']).apply(
        cluster_clv_sites, CUTOFF).reset_index(drop=True)

    hdf = gdf[['seqname', 'strand', 'mclv']].drop_duplicates()

    df_clv_sc = pd.read_csv('/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/data/reference_data/annotated-clv-sc-mapping.csv.gz')
    annot_clvs = df_clv_sc.groupby(['seqname', 'strand']).apply(lambda g: g.clv.sort_values().values)

    def calc_abs_dist_to_annot_clv(grp):
    #     print(grp.name)
        if grp.name in annot_clvs.index:
             # grp.name holds the group key
            aclvs = annot_clvs.loc[grp.name]
            aclvs = annot_clvs.loc[grp.name]  # grp.name holds the group key
            bcast = np.broadcast_to(grp.clv.values, (aclvs.shape[0], grp.shape[0])).T
            grp['abs_dist_to_aclv'] = np.min(np.abs(bcast - aclvs), axis=1)
            return grp
        else:
            # return None, this group will be gone, which is good as they won't
            # be clv of targeted genes for sure, increasing precision
            return 

    idf = hdf.rename(columns={'mclv': 'clv'}).groupby(['seqname', 'strand']).apply(
        calc_abs_dist_to_annot_clv).reset_index(drop=True)

    def search_hexamer_wrapper(refseq, chrm, clv, strand, window=50):
        res = search_hexamer.search(refseq, chrm, clv, strand, window)
        if res is None:
            res = ['NA', -1, -1]
        return pd.Series(res, index=['hexamer', 'hexamer_id', 'hexamer_loc'])

    REF_FA = '/gsc/www/bcgsc.ca/downloads/zxue/tasrkleat-static/off-cloud/reference_data/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa'
    refseq = pysam.FastaFile(REF_FA)

    hxm_df = idf.apply(
        lambda row: search_hexamer_wrapper(
            refseq, row.seqname, row.clv, row.strand), axis=1)
    idf_hxm = pd.concat([idf, hxm_df], axis=1)

    outdir = f'./{kleat_version}_ml/cv/trained_on_{train_data_sample_id}/max_depth{max_depth}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    out_csv = f'{outdir}/{SAMPLE_ID}.csv'
    idf_hxm.rename(columns={'mclv': 'clv'}).to_csv(out_csv, index=False)
