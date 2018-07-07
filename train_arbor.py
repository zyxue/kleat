import os
import sys
import pickle

import pandas as pd
from sklearn.tree import DecisionTreeClassifier

from ml_utils import load_polya_seq_df, map_clvs, KARBOR_FEATURE_COLS


kleat_version = 'karbor'
map_cutoff = 50
sample_id = sys.argv[1]

fdf = pd.read_csv(f'./benchmark_114genes/{kleat_version}_ml/{sample_id}_ml_ready.csv', low_memory=False)

df_ref = load_polya_seq_df(sample_id)

fdf_mapped = map_clvs(fdf, df_ref)

# the row orders are different from fdf
fdf_mapped['is_tp'] = fdf_mapped.abs_dist < map_cutoff

Xs = fdf_mapped[KARBOR_FEATURE_COLS]

depths = range(2, 20, 1)
for max_depth in depths:
    print(f'working on depth: {max_depth}')
    clf = DecisionTreeClassifier(max_depth=max_depth)
    clf.fit(Xs, fdf_mapped.is_tp.values)

    fdf['predicted'] = clf.predict(Xs)

    outdir = f'./benchmark_114genes/{kleat_version}_ml/pkl/DT/{sample_id}/max_depth{max_depth}'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    out_pkl = f'{outdir}/clf.pkl'
    print(f'writing {out_pkl}')
    with open(out_pkl, 'wb') as opf:
        pickle.dump(clf, opf)
