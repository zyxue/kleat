import sys
import pickle

import numpy as np
import pandas as pd

import pysam

from cluster import cluster_clv_sites
from utils import FEATURE_COLS

SAMPLE_ID = sys.argv[1]         # e.g. HBRC4
print(f'working on {SAMPLE_ID}')

for max_depth in list(range(2, 10)):
    print(f'working on depth: {max_depth}')
    classifier_pkl = f'./pkl/DT_HBRC4_max_depth{max_depth}.pkl'

    edf = pd.read_csv(f'./kleat3_{SAMPLE_ID}_ml_ready.csv', low_memory=False)

    with open(classifier_pkl, 'rb') as inf:
        clf = pickle.load(inf)

    fdf = edf

    fdf['predicted'] = clf.predict(fdf[FEATURE_COLS])

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
            aclvs = annot_clvs.loc[grp.name] # grp.name holds the group key
            bcast = np.broadcast_to(grp.clv.values, (aclvs.shape[0], grp.shape[0])).T
            grp['abs_dist_to_aclv'] = np.min(np.abs(bcast - aclvs), axis=1)
            return grp
        else:
            # return None, this group will be gone, which is good 
            # as they won't be clv of targeted genes for sure, increasing precision
            return 

    idf = hdf.rename(columns={'mclv': 'clv'}).groupby(['seqname', 'strand']).apply(
        calc_abs_dist_to_annot_clv).reset_index(drop=True)

    import search_hexamer

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

    out_csv = f'./cv/kleat3_{SAMPLE_ID}_ml_filtered_clustered_max_depth{max_depth}.csv'
    idf_hxm.rename(columns={'mclv': 'clv'}).to_csv(out_csv, index=False)
