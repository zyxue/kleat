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
from kleat.post import aggregate_polya_evidence, add_annot_info
from kleat.misc import utils as U
from kleat.misc import settings as S


logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')

logger = logging.getLogger(__name__)


def collect_polya_evidence(c2g_bam, r2c_bam, ref_fa, csvwriter):
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

        for rec in process_bridge_and_link(contig, r2c_bam,
                                           ref_fa, csvwriter):
            # TODO: with either bridge or link, they probably won't support
            # clv of the other strand
            ascs.append(gen_key(rec))

        if not apautils.has_tail(contig):
            process_blank(contig, ref_fa, csvwriter, ascs)


def gen_tmp_output(output):
    return os.path.join(
        os.path.dirname(output), '__tmp_{0}'.format(os.path.basename(output)))


def main():
    args = get_args()
    c2g_bam = pysam.AlignmentFile(args.contig_to_genome)
    r2c_bam = pysam.AlignmentFile(args.read_to_contig)
    ref_fa = pysam.FastaFile(args.reference_genome)
    output = os.path.abspath(args.output)
    tmp_output = gen_tmp_output(output)
    for o in [output, tmp_output]:
        if os.path.exists(o):
            U.backup_file(o)

    with open(tmp_output, 'wt') as opf:
        csvwriter = csv.writer(opf, delimiter='\t')
        csvwriter.writerow(S.HEADER)
        collect_polya_evidence(c2g_bam, r2c_bam, ref_fa, csvwriter)

    logger.info('Reading {0} into a pandas.DataFrame...'.format(tmp_output))
    df_clv = U.timeit(pd.read_csv)(tmp_output, keep_default_na=False, sep='\t')
    logger.info('df.shape: {0}'.format(df_clv.shape))

    # TODO: slow step, maybe parallize it
    logger.info('Aggregating polya evidence for each (seqname, strand, clv)...')
    df_agg = aggregate_polya_evidence(df_clv, args.num_cpus)

    logger.info('Calculating closest annotated clv...')
    df_clv_with_adist = add_annot_info(df_agg, args.karbor_clv_annotation)

    logger.info('Writing to {0}...'.format(output))
    df_clv_with_adist.to_csv(output, sep='\t', index=False)

    if not args.keep_pre_aggregation_tmp_file:
        os.remove(tmp_output)


if __name__ == "__main__":
    main()
