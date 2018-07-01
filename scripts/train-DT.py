import os
import sys
import pickle

import pandas as pd
from sklearn.tree import DecisionTreeClassifier

from compare_to_polyA_Seq import load_polya_df, map_clvs

kleat_version = sys.argv[1]
SAMPLE_ID = sys.argv[2]

if kleat_version == 'kleat2':
    from utils.utils import KLEAT2_FEATURE_COLS as feature_cols
elif kleat_version == 'kleat3':
    from utils.utils import KLEAT3_FEATURE_COLS as feature_cols
else:
    raise

fdf = pd.read_csv(f'./{kleat_version}_ml/{SAMPLE_ID}_ml_ready.csv', low_memory=False)

rootdir = '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk'
truth_csv = f'{rootdir}/{SAMPLE_ID[:3]}/{SAMPLE_ID[3:]}/polyA-Seq/polyA-Seq-truth-114-genes.csv'
df_ref = load_polya_df(truth_csv)
df_ref['seqname'] = df_ref.seqname.str.replace('chr', '').replace('M', 'MT')

fdf_mapped = map_clvs(fdf, df_ref)
DIST_CUTOFF = 50
# the row orders are different from fdf
fdf_mapped['is_tp'] = fdf_mapped.abs_dist < DIST_CUTOFF

Xs = fdf_mapped[feature_cols]


depths = range(2, 20, 1)
for max_depth in depths:
    print(f'working on depth: {max_depth}')
    clf = DecisionTreeClassifier(max_depth=max_depth)
    clf.fit(Xs, fdf_mapped.is_tp.values)

    fdf['predicted'] = clf.predict(Xs)

    outdir = f'./{kleat_version}_ml/pkl/DT/{SAMPLE_ID}/max_depth{max_depth}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    out_pkl = f'{outdir}/clf.pkl'
    print(f'writing {out_pkl}')
    with open(out_pkl, 'wb') as opf:
        pickle.dump(clf, opf)
