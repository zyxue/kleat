#!/usr/bin/env python

import os
import logging
import multiprocessing

import pandas as pd

from kleat import polya
from kleat.args import get_args
from kleat.post import (
    cluster_clv_parallel,
    aggregate_polya_evidence,
    add_annot_info,
    add_hex_dist,
    add_extra
)
from kleat.misc import utils as U
from kleat.misc import settings as S

logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')

logger = logging.getLogger(__name__)


def gen_output(args_output, output_format):
    format2extension_dd = {
        'csv': 'csv',
        'tsv': 'tsv',
        'pickle': 'pkl',
        'pkl': 'pkl',
    }
    if args_output is None:
        args_output = './output.{0}'.format(format2extension_dd[output_format])
    return os.path.abspath(args_output)


def dump_output_df(out_df, output, output_format):
    if output_format == 'csv':
        out_df.to_csv(output, index=False)
    elif output_format == 'tsv':
        out_df.to_csv(output, sep='\t', index=False)
    elif output_format in ['pkl', 'pickle']:
        out_df.to_pickle(output)
    else:
        raise ValueError('unknown output format: {0}'.format(output_format))


def main():
    args = get_args()
    c2g_bam_file = args.contigs_to_genome
    r2c_bam_file = args.reads_to_contigs
    ref_fa_file = args.reference_genome
    output = gen_output(args.output, args.output_format.lower())
    U.backup_file(output)

    args_list = polya.prepare_args_for_collect_polya_evidence(
        args.num_cpus, output, c2g_bam_file,
        r2c_bam_file, ref_fa_file, args.bridge_skip_check_size
    )

    logger.info('Processing contigs in parallel with {0} CPUs...'.format(args.num_cpus))
    with multiprocessing.Pool(args.num_cpus) as p:
        tmp_tsv_files = p.map(polya.collect_polya_evidence_wrapper, args_list)

    logger.info('Reading {0} files into a single pandas.DataFrame...'.format(len(tmp_tsv_files)))
    dfs = []
    for f in tmp_tsv_files:
        _df = pd.read_csv(f, keep_default_na=False, sep='\t')
        dfs.append(_df)
    df_clv = pd.concat(dfs)
    logger.info('df.shape: {0}'.format(df_clv.shape))

    # remove file in a separate loop in case exit none-zero and still have
    # partial results for debugging
    logger.info('removing {0} tmp files ...'.format(len(tmp_tsv_files)))
    for f in tmp_tsv_files:
        os.remove(f)

    if args.keep_pre_aggregation_tmp_file:
        tmp_output = polya.gen_tmp_output(output)
        U.backup_file(tmp_output)
        logger.info('Dumping raw results before aggregation to {0}'.format(tmp_output))
        df_clv.to_csv(tmp_output, sep='\t', index=False)

    if args.cluster_first_then_aggregate:
        logger.info('Clustering clv since --cluster-first-then-aggregate is specified ...')
        df_clustered = cluster_clv_parallel(df_clv, args.cluster_cutoff, args.num_cpus)
        df_clustered['clv'] = df_clustered['mode_clv']
        df_clustered.drop(['cluster_id', 'mode_clv'], axis=1, inplace=True)
        df_clv = df_clustered

    logger.info('Aggregating polya evidence for each (seqname, strand, clv)...')
    df_agg = aggregate_polya_evidence(df_clv, args.num_cpus)

    logger.info('Calculating closest annotated clv...')
    df_ant_dist = add_annot_info(df_agg, args.karbor_clv_annotation)

    logger.info('calculating distance between PAS hexamers and clvs ...')
    df_hex_dist = add_hex_dist(df_ant_dist)
    add_extra(df_hex_dist)

    logger.info('Writing to {0}...'.format(output))
    out_df = df_hex_dist.rename(columns=S.FORMAT_OUTPUT_HEADER_DD)

    out_df = out_df[S.OUTPUT_HEADER]
    out_df.sort_values(['seqname', 'strand', 'clv'], inplace=True)
    dump_output_df(out_df, output, args.output_format)
    logger.info('Completed writing to {0}...'.format(output))


if __name__ == "__main__":
    main()
