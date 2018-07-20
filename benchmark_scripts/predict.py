import os
import argparse
import csv
import pickle
import logging
import multiprocessing
import subprocess


import pandas as pd
from sklearn.tree import DecisionTreeClassifier, export_graphviz

from kleat.misc.utils import backup_file
from ml_utils import load_polya_seq_df, map_clvs, compare, KARBOR_FEATURE_COLS
import train_arbor as TA
from cluster import cluster_clv_sites




logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')


def get_args():
    parser = argparse.ArgumentParser(
        description=('train & test a decision tree for different max_depths '
                     'for UHRC1/2 and HBRC4/6 samples.')
    )

    parser.add_argument(
        '--classifier-pkl', type=str,
        help='only necessary when --no-training is specified'
    )

    parser.add_argument(
        '-e', '--test-sample-ids', type=str, nargs='+', required=True,
        help='the sample IDs to test the tree on, e.g. UHRC1 UHRC2 HBRC4 HBRC6'
    )

    parser.add_argument(
        '-o', '--output', type=str, default='./out.csv',
        help='output for recording testing results in csv'
    )

    # arguments below are often good in default
    parser.add_argument(
        '--map-cutoff', type=int, default=50,
        help=('the cutoff below which the predicted clv and the reference clv '
              '(e.g. polyA-Seq) are considered the same')
    )

    parser.add_argument(
        '--num-cpus', type=int,
        default=min(multiprocessing.cpu_count(), 16),
        help='the number of CPUs to use for parallel training and testing'
    )
    return parser.parse_args()


def main():
    args = get_args()

    with open(args.classifier_pkl, 'rb') as inf:
        clf = pickle.load(inf)

    with open(args.output, 'wt') as opf:
        csvwriter = csv.writer(opf)
        csvwriter.writerow(
            ['sample_id', 'precision', 'recall', 'f1', 'tree_max_depth']
        )

        for test_sample_id in args.test_sample_ids:
            df_te = TA.load_df(test_sample_id)
            df_ref = load_polya_seq_df(test_sample_id)
            row_results = TA.run_test(test_sample_id, df_te, df_ref, clf, args.map_cutoff)

            csvwriter.writerow(row_results)


if __name__ == "__main__":
    main()
