import pickle

import pandas as pd
from sklearn import metrics
from sklearn.tree import DecisionTreeClassifier

from compare_to_polyA_Seq import load_polya_df, map_clvs


SAMPLE_ID = 'HBRC4'
edf = pd.read_csv(f'./kleat3_{SAMPLE_ID}_ml_ready.csv', low_memory=False)

FEATURE_COLS = [
    'num_contig_tail_reads', 'num_bridge_reads', 'num_link_reads',
    'num_blank_contigs',
    'contig_len',
    'contig_mapq', 'contig_tail_length', 'max_bridge_tail_length',
    'abs_dist_to_aclv'] + [
        'AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA',
        'CATAAA', 'GATAAA', 'AATATA', 'AATACA',
        'AATAGA', 'AAAAAG', 'ACTAAA', 'AAGAAA',
        'AATGAA', 'TTTAAA', 'AAAACA', 'GGGGCT', 'NA'
    ]

rootdir = '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk'
truth_csv = f'{rootdir}/{SAMPLE_ID[:3]}/{SAMPLE_ID[3:]}/polyA-Seq/polyA-Seq-truth-114-genes.csv'
df_ref = load_polya_df(truth_csv)
df_ref['seqname'] = df_ref.seqname.str.replace('chr', '').replace('M', 'MT')

fdf = edf.copy()
fdf_mapped = map_clvs(fdf, df_ref)
DIST_CUTOFF = 50
fdf_mapped['is_tp'] = fdf_mapped.abs_dist < DIST_CUTOFF

Xs = fdf[FEATURE_COLS]


for max_depth in list(range(2, 10)):
    print(f'working on depth: {max_depth}')
    clf = DecisionTreeClassifier(max_depth=max_depth)
    clf.fit(Xs, fdf_mapped.is_tp.values)

    fdf['predicted'] = clf.predict(Xs)
    with open(f'./pkl/DT_{SAMPLE_ID}_max_depth{max_depth}.pkl', 'wb') as opf:
        pickle.dump(clf, opf)
