#!/usr/bin/env python

import os
import logging
import csv

import pysam
from tqdm import tqdm
import pandas as pd

from kleat.misc import apautils
from kleat.args import get_args
from kleat.proc import process_suffix, process_bridge_and_link, process_blank
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


def collect_polya_evidence(c2g_bam, r2c_bam, ref_fa, csvwriter, bridge_skip_check_size):
    """loop through each contig and collect polyA evidence"""
    gen_key = apautils.gen_clv_key_tuple_from_clv_record
    iters = tqdm(enumerate(c2g_bam), desc='processed', unit=' contigs')
    for k, contig in iters:
        if contig.is_unmapped:
            continue

        ascs = []           # already supported clvs
        rec = process_suffix(
            contig, r2c_bam, ref_fa, csvwriter)
        if rec is not None:
            ascs.append(gen_key(rec))

        for rec in process_bridge_and_link(
                contig, r2c_bam, ref_fa, csvwriter, bridge_skip_check_size):
            # TODO: with either bridge or link, they probably won't support
            # clv of the other strand
            ascs.append(gen_key(rec))

        if not apautils.has_tail(contig):
            process_blank(contig, ref_fa, csvwriter, ascs)


def gen_tmp_output(output):
    return os.path.join(
        os.path.dirname(output), '__tmp_{0}.tsv'.format(os.path.basename(output)))


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
    c2g_bam = pysam.AlignmentFile(args.contigs_to_genome)
    r2c_bam = pysam.AlignmentFile(args.reads_to_contigs)
    ref_fa = pysam.FastaFile(args.reference_genome)
    output = gen_output(args.output, args.output_format.lower())
    tmp_output = gen_tmp_output(output)
    for o in [output, tmp_output]:
        if os.path.exists(o):
            U.backup_file(o)

    with open(tmp_output, 'wt') as opf:
        csvwriter = csv.writer(opf, delimiter='\t')
        csvwriter.writerow(S.HEADER)
        collect_polya_evidence(
            c2g_bam, r2c_bam, ref_fa, csvwriter,
            args.bridge_skip_check_size
        )

    logger.info('Reading {0} into a pandas.DataFrame...'.format(tmp_output))
    df_clv = U.timeit(pd.read_csv)(tmp_output, keep_default_na=False, sep='\t')
    logger.info('df.shape: {0}'.format(df_clv.shape))

    logger.info('Clustering clv ...')
    df_clustered = cluster_clv_parallel(df_clv, args.cluster_cutoff, args.num_cpus)
    df_clustered['clv'] = df_clustered['mode_clv']
    df_clustered.drop(['cluster_id', 'mode_clv'], axis=1, inplace=True)

    logger.info('Aggregating polya evidence for each (seqname, strand, clv)...')
    df_agg = aggregate_polya_evidence(df_clustered, args.num_cpus)

    logger.info('Calculating closest annotated clv...')
    df_ant_dist = add_annot_info(df_agg, args.karbor_clv_annotation)

    logger.info('calculating distance between PAS hexamers and clvs ...')
    df_hex_dist = add_hex_dist(df_ant_dist)

    logger.info('Writing to {0}...'.format(output))
    out_df = df_hex_dist.rename(columns=S.FORMAT_OUTPUT_HEADER_DD)
    add_extra(out_df)
    out_df = out_df[S.OUTPUT_HEADER]
    dump_output_df(out_df, output, args.output_format)

    if not args.keep_pre_aggregation_tmp_file:
        os.remove(tmp_output)


if __name__ == "__main__":
    main()
