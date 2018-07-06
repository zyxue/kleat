#!/usr/bin/env python

import os
import logging
import csv
import argparse

import pysam
from tqdm import tqdm
import pandas as pd
import numpy as np

from kleat.evidence import suffix, bridge, link, blank
from kleat.misc import apautils
from kleat.misc import utils as U
from kleat.misc import settings as S


logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')

logger = logging.getLogger(__name__)


def get_args():
    parser = argparse.ArgumentParser(
        description='KLEAT: cleavage site detection via de novo assembly')
    parser.add_argument(
        '-c', '--contig-to-genome', type=str, required=True,
        help='input contig-to-genome alignment BAM file'
    )
    parser.add_argument(
        '-r', '--read-to-contig', type=str, required=True,
        help='input read-to-contig alignment BAM file'
    )
    parser.add_argument(
        '-f', '--reference-genome', type=str, default=None,
        help=('(optional) reference genome FASTA file, if provided, '
              'KLEAT will search polyadenylation signal (PAS) hexamer in '
              'both contig and reference genome, which is useful for '
              'checking mutations that may affect PAS hexmaer.  '
              'Note this fasta file needs to be consistent with the one '
              'used for generating the read-to-contig BAM alignments')
    )
    parser.add_argument(
        '-m', '--clv-sc-mapping', type=str, required=True,
        help=('the mapping pickle (TODO: support CSV) file of clv-to-stop '
              'codon extracted from annotation')
    )
    parser.add_argument(
        '-o', '--output', type=str, default='./output.tsv',
        help='output tsv file'
    )
    return parser.parse_args()


def process_suffix(contig, r2c_bam, ref_fa, csvwriter):
    tail_side = apautils.has_tail(contig)
    if tail_side is not None:
        clv_record = suffix.gen_clv_record(contig, r2c_bam, tail_side, ref_fa)
        apautils.write_row(clv_record, csvwriter)
        return clv_record


def process_bridge_and_link(contig, r2c_bam, ref_fa, csvwriter):
    # bridge & link
    aligned_reads = r2c_bam.fetch(contig.query_name)
    dd_bridge, dd_link = extract_bridge_and_link(contig, aligned_reads)
    bdg_clvs, lnk_clvs = [], []
    if len(dd_bridge['num_reads']) > 0:
        bdg_clvs.extend(
            bridge.write_evidence(dd_bridge, contig, ref_fa, csvwriter))
    if len(dd_link['num_reads']) > 0:
        lnk_clvs.extend(
            link.write_evidence(dd_link, contig, ref_fa, csvwriter))
    return bdg_clvs + lnk_clvs


def process_blank(contig, ref_fa, csvwriter, asp_clv_keys):
    """:param asp_clv_keys: list of already supported clv key tuples"""
    asp_clv_keys = set(asp_clv_keys)
    for clv_rec in blank.gen_two_clv_records(contig, ref_fa, asp_clv_keys):
        apautils.write_row(clv_rec, csvwriter)


def extract_bridge_and_link(contig, aligned_reads):
    """
    bridge and link are processed together by looping throught reads aligned to
    the contig

    :param contig: a pysam.libcalignedsegment.AlignedSegment instance
    :param aligned_reads: a pysam.libcalignmentfile.IteratorRowRegion instance
    """
    dd_bridge = bridge.init_evidence_holder()
    dd_link = link.init_evidence_holder()
    for read in aligned_reads:
        if bridge.is_a_bridge_read(read):
            bdg_evid = bridge.analyze_bridge(contig, read, dd_bridge)
            bridge.update_evidence(bdg_evid, dd_bridge)
        elif link.is_a_link_read(read):
            link_evid = link.analyze_link(contig, read)
            link.update_evidence(link_evid, dd_link)
    return dd_bridge, dd_link


def gen_ref_fa(ref_genome_file):
    if ref_genome_file is not None:
        return pysam.FastaFile(ref_genome_file)


def calc_abs_dist_to_annot_clv(grp, annot_clvs):
    aclvs = annot_clvs.loc[grp.name] # grp.name holds the group key
    bcast = np.broadcast_to(grp.clv.values, (aclvs.shape[0], grp.shape[0])).T

    sgn_dists = bcast - aclvs   # sgn: signed
    abs_dists = np.abs(sgn_dists)
    min_idxes = np.argmin(abs_dists, axis=1)

    nrows = abs_dists.shape[0]
    dists = sgn_dists[np.arange(nrows), min_idxes]

    grp['aclv'] = aclvs[min_idxes]
    grp['signed_dist_to_aclv'] = dists
    return grp


def add_abs_dist_to_annot_clv(df_clv, df_mapping):
    """
    add absolute distance to the closest annotated clv as an addition column

    :param df_mapping: clv-stop codon mapping in dataframe
    """
    # do some checking about which version of seqnames are used, use whatever
    # is used by df_clv as the reference
    mapping_seqname_is_ucsc = df_mapping.seqname.values[0] in S.UCSC_SEQNAMES
    cleavge_seqname_is_ucsc = df_clv.seqname.values[0] in S.UCSC_SEQNAMES

    if mapping_seqname_is_ucsc != cleavge_seqname_is_ucsc:
        if mapping_seqname_is_ucsc:
            df_mapping.seqname = df_mapping.seqname.replace(S.UCSC_TO_ENSEMBL_SEQNAME)
        else:
            df_mapping.seqname = df_mapping.seqname.replace(S.ENSEMBL_TO_UCSC_SEQNAME)

    annot_clvs = df_mapping.groupby(['seqname', 'strand']).apply(
        lambda g: g.clv.sort_values().values)

    # remove patch chromosomes
    if cleavge_seqname_is_ucsc:
        ndf_clv = df_clv.query('seqname in {0}'.format(S.UCSC_SEQNAMES))
    else:
        ndf_clv = df_clv.query('seqname in {0}'.format(S.ENSEMBL_SEQNAMES))

    logger.info('calculating absolute distances to annotated cleavage sites')
    timed = U.timeit(
        lambda _df: _df.groupby(['seqname', 'strand'])\
        .apply(calc_abs_dist_to_annot_clv, annot_clvs=annot_clvs)
    )
    out = timed(ndf_clv)
    return out


def main():
    args = get_args()
    c2g_bam = pysam.AlignmentFile(args.contig_to_genome)
    r2c_bam = pysam.AlignmentFile(args.read_to_contig)
    ref_fa = gen_ref_fa(args.reference_genome)
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

    logger.info('reading {0}'.format(clv_sc_mapping))
    df_mapping = pd.read_pickle(clv_sc_mapping)
    logger.info('df.shape: {0}'.format(df_mapping.shape))

    logger.info('reading {0} into a pandas.DataFrame'.format(tmp_output))
    df_clv = U.timeit(pd.read_csv)(tmp_output, sep='\t')
    logger.info('df.shape: {0}'.format(df_clv.shape))

    df_clv_with_adist = add_abs_dist_to_annot_clv(df_clv, df_mapping)
    df_clv_with_adist.to_csv(output, sep='\t', index=False)

    # TODO: remove tmp_output


if __name__ == "__main__":
    main()
