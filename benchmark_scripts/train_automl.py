import os
import argparse
import csv
import pickle
import logging
import multiprocessing
import subprocess

import autosklearn.classification

import pandas as pd

from ml_utils import load_polya_seq_df, map_clvs, compare, KARBOR_FEATURE_COLS
from kleat.misc.cluster import cluster_clv_sites
from kleat.misc.utils import backup_file


logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')


def load_df(sample_id):
    # kind of hacky
    if sample_id.endswith('114genes'):
        sample_id = sample_id.replace('_114genes', '')
        infile = os.path.abspath(f'./benchmark_114genes/{sample_id}.ml_ready.pkl')
    else:
        infile = os.path.abspath(f'./benchmark_transcriptome/{sample_id}.ml_ready.pkl')
    logging.info(f'reading {infile} ...')
    df = pd.read_pickle(infile)
    return df


def map_to_ref(df, df_ref, map_cutoff):
    logging.info('mapping predicted clv to ground truth ...')
    df = map_clvs(df, df_ref)
    df['is_tp'] = df.abs_dist < map_cutoff
    return df


def prepare_args(df_tr_mapped, max_depth_list, clf_type):
    Xs = df_tr_mapped[KARBOR_FEATURE_COLS]
    ys = df_tr_mapped.is_tp.values
    res = []
    for d in max_depth_list:
        clf = autosklearn.classification.AutoSklearnClassifier(ml_memory_limit=30720)
        res.append([clf, Xs, ys])
    return res


def train_it(clf, Xs, ys):
    logging.info(f'training {clf}')
    clf.fit(Xs, ys)
    return clf


def train_it_wrapper(args):
    return train_it(*args)


def predict(df_te, clf):
    Xs_te = df_te[KARBOR_FEATURE_COLS]
    df_te['predicted'] = clf.predict(Xs_te)

    # dedup here may not be necessary
    clv_key_cols = ['seqname', 'strand', 'clv']
    out_df = df_te.query('predicted')[clv_key_cols].drop_duplicates()
    return out_df


def cluster_clv(df, cutoff=20):
    """
    Serial version, as tests are already paralleled, otherwise

        AssertionError: daemonic processes are not allowed to have children

    see max_recall_analyiss.cluster_clv_parallel for a parallel version

    return clustered clv in new dataframe with three columns:

    - seqname
    - strand
    - clv, i.e. the representative mode clv for each cluster
    """
    grped = df.groupby(['seqname', 'strand'])
    applied = grped.apply(cluster_clv_sites, cutoff)
    dedupped = applied[['seqname', 'strand', 'mode_clv']].drop_duplicates()
    out = dedupped.rename(columns={'mode_clv': 'clv'}).reset_index(drop=True)
    return out


def run_test(test_sample_id, df_te, df_ref, clf, map_cutoff):
    """args as list, otherwise this function can't be passed
    multiprocessing.Pool"""
    logging.info(f'testing on {test_sample_id} with {clf}')

    df_predicted = predict(df_te, clf)
    df_clustered = cluster_clv(df_predicted)

    recall, precision, f1 = compare(df_clustered, df_ref, map_cutoff)
    return test_sample_id, precision, recall, f1


def run_test_wrapper(args):
    return run_test(*args)


def run_test_in_parallel(clf_list, test_sample_ids, map_cutoff, num_cpus):
    all_results = []
    for test_sample_id in test_sample_ids:
        df_te = load_df(test_sample_id)
        df_ref = load_polya_seq_df(test_sample_id)

        args_list_for_test = []
        for clf in clf_list:
            args_list_for_test.append(
                (test_sample_id, df_te, df_ref, clf, map_cutoff)
            )

        with multiprocessing.Pool(num_cpus) as p:
            logging.info(f'start parallel testing with {num_cpus} CPUs ...')
            results = p.map(run_test_wrapper, args_list_for_test)
            all_results.extend(results)
    return all_results


def gen_clf_output(output, classifier_type, depth):
    """
    generate path to output clf pickl based on output_csv and depth of the tree
    """
    out_pkl = os.path.join(
        os.path.dirname(output),
        os.path.basename(output).replace('.csv', f'.{classifier_type}.tree_max_depth_{depth}.pkl')
    )
    return out_pkl


def gen_tree_vis_outputs(output, classifier_type, depth):
    out_pkl = gen_clf_output(output, classifier_type, depth)
    out_dot = out_pkl[:-3] + 'dot'
    out_png = out_pkl[:-3] + 'png'
    return out_dot, out_png


def pickle_clf(clf, out_pkl):
    with open(out_pkl, 'wb') as opf:
        pickle.dump(clf, opf)


def get_args():
    parser = argparse.ArgumentParser(
        description=('train & test a decision tree for different max_depths '
                     'for UHRC1/2 and HBRC4/6 samples.')
    )

    parser.add_argument(
        '-r', '--train-sample-id', type=str, required=True,
        help='the sample ID to train a decision tree on, e.g. HBRC4'
    )
    parser.add_argument(
        '-e', '--test-sample-ids', type=str, nargs='+', required=True,
        help='the sample IDs to test the tree on, e.g. UHRC1 UHRC2 HBRC4 HBRC6'
    )
    parser.add_argument(
        '-d', '--max-depths', type=int, nargs=3, default=[2, 20, 1],
        help='the range of max_depth values for training the tree, [beg, end, step]'
    )
    # just pickle them all
    # parser.add_argument(
    #     '-p', '--pickle-depths', type=int, nargs='+', default=[],
    #     help=('for the particular max_depth values, '
    #           'the tree would be pickled and a visualization would be produced ')
    # )
    parser.add_argument(
        '-o', '--output', type=str, default='./out.csv',
        help='output for recording testing results in csv'
    )
    parser.add_argument(
        '-t', '--classifier-type', type=str, default='DecisionTreeClassifier',
        choices=['DecisionTreeClassifier', 'GradientBoostingClassifier'],
        help=('the tree-based classifier type intended, the name follows '
              'the corresponding classifier class name in scikit-learn')
    )
    parser.add_argument(
        '--no-testing', action='store_true',
        help='if specified, only training is conducted, no testing would be run'
    )

    # arguments below are often good in default
    parser.add_argument(
        '--map-cutoff', type=int, default=50,
        help=('the cutoff below which the predicted clv and the reference clv '
              '(e.g. polyA-Seq) are considered the same')
    )
    parser.add_argument(
        '--vis-tree-max-depth', type=int, default=3,
        help=('specified max_depth when visualizing the tree, '
              'a high value would result in a big tree')
    )
    parser.add_argument(
        '--num-cpus', type=int,
        default=min(multiprocessing.cpu_count(), 16),
        help='the number of CPUs to use for parallel training and testing'
    )
    return parser.parse_args()


def main():
    args = get_args()
    train_sample_id = args.train_sample_id
    test_sample_ids = args.test_sample_ids
    clf_type = args.classifier_type

    df_tr = load_df(train_sample_id)
    df_tr_ref = load_polya_seq_df(train_sample_id)
    df_tr_mapped = map_to_ref(df_tr, df_tr_ref, args.map_cutoff)

    beg, end, step = args.max_depths
    max_depth_list = range(beg, end, step)
    logging.info(f'max_depth list trees: {max_depth_list}')
    train_args = prepare_args(df_tr_mapped, max_depth_list, clf_type)

    logging.info('prepare for TRAINing ...')
    with multiprocessing.Pool(args.num_cpus) as p:
        logging.info(f'start parallel training with {args.num_cpus} CPUs ...')
        clf_list = [train_it_wrapper(train_args[0])]

    clf_dd = dict(zip(max_depth_list, clf_list))
    for depth in clf_dd:
        clf_pkl = gen_clf_output(args.output, clf_type, depth)
        backup_file(clf_pkl)
        logging.info(f'pickling {clf_pkl}')
        pickle_clf(clf_dd[depth], clf_pkl)

    if args.no_testing:
        logging.info('--no-testing is specified, no testing will be run.')
    else:
        logging.info('prepare for TESTing ...')
        backup_file(args.output)
        with open(args.output, 'wt') as opf:
            csvwriter = csv.writer(opf)
            csvwriter.writerow(
                ['sample_id', 'precision', 'recall', 'f1']
            )

            test_results = run_test_in_parallel(
                clf_list, test_sample_ids, args.map_cutoff, args.num_cpus)

            for row in test_results:
                csvwriter.writerow(row)


if __name__ == "__main__":
    main()
