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
from kleat.post import agg_polya_evidence, add_abs_dist_to_annot_clv
from kleat.misc import utils as U
from kleat.misc import settings as S


logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')

logger = logging.getLogger(__name__)


def main():
    args = get_args()
    c2g_bam = pysam.AlignmentFile(args.contig_to_genome)
    r2c_bam = pysam.AlignmentFile(args.read_to_contig)
    ref_fa = pysam.FastaFile(args.reference_genome)
    clv_sc_mapping = args.clv_sc_mapping
    output = args.output
    tmp_output = os.path.join(os.path.dirname(output),
                              'tmp_{0}'.format(os.path.basename(output)))
    if os.path.exists(output):
        U.backup_file(output)

    with open(tmp_output, 'wt') as opf:
        csvwriter = csv.writer(opf, delimiter='\t')
        csvwriter.writerow(S.HEADER)

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

    logger.info('Reading {0} into a pandas.DataFrame...'.format(tmp_output))
    df_clv = U.timeit(pd.read_csv)(tmp_output, keep_default_na=False, sep='\t')
    logger.info('df.shape: {0}'.format(df_clv.shape))

    # TODO: slow step, maybe parallize it
    logger.info('Aggregating evidence for the same (seqname, strand, clv)...')
    grped = df_clv.groupby(['seqname', 'strand', 'clv'])
    pre_df_clv_agg = U.timeit(grped.apply)(agg_polya_evidence)
    df_clv_agg = pre_df_clv_agg.reset_index()[S.HEADER]

    logger.info('Reading {0}'.format(clv_sc_mapping))
    df_mapping = pd.read_pickle(clv_sc_mapping)
    logger.info('df.shape: {0}'.format(df_mapping.shape))

    logger.info('Calculating closest annotated clv...')
    df_clv_with_adist = add_abs_dist_to_annot_clv(df_clv_agg, df_mapping)

    logger.info('Writing to {0}...'.format(output))
    df_clv_with_adist.to_csv(output, sep='\t', index=False)

    # TODO: remove tmp_output


if __name__ == "__main__":
    main()
