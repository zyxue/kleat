#!/usr/bin/env python

import os
import logging
import csv
import argparse
from collections import defaultdict

import pysam
from tqdm import tqdm

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
        '-o', '--output', type=str, default='./output.tsv',
        help='output tsv file'
    )
    return parser.parse_args()


def process_suffix(contig, r2c_bam, ref_fa, csvwriter):
    tail_side = apautils.has_tail(contig)
    if tail_side is not None:
        clv_record = suffix.gen_clv_record(contig, r2c_bam, tail_side, ref_fa)
        apautils.write_row(clv_record, csvwriter)


def process_bridge_and_link(contig, r2c_bam, ref_fa, csvwriter):
    # bridge & link
    aligned_reads = r2c_bam.fetch(contig.query_name)
    dd_bridge, dd_link = extract_bridge_and_link(contig, aligned_reads)
    if len(dd_bridge['num_reads']) > 0:
        bridge.write_evidence(dd_bridge, contig, ref_fa, csvwriter)
    if len(dd_link['num_reads']) > 0:
        link.write_evidence(dd_link, contig, ref_fa, csvwriter)


def process_blank(contig, csvwriter):
    for clv_rec in blank.gen_two_clv_records(contig):
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
            bdg_evid = bridge.analyze_bridge(contig, read)
            bridge.update_evidence(bdg_evid, dd_bridge)
        elif link.is_a_link_read(read):
            link_evid = link.analyze_link(contig, read)
            link.update_evidence(link_evid, dd_link)
    return dd_bridge, dd_link


def gen_ref_fa(ref_genome_file):
    if ref_genome_file is not None:
        return pysam.FastaFile(ref_genome_file)


def main():
    args = get_args()
    c2g_bam = pysam.AlignmentFile(args.contig_to_genome)
    r2c_bam = pysam.AlignmentFile(args.read_to_contig)
    ref_fa = gen_ref_fa(args.reference_genome)
    output = args.output
    if os.path.exists(output):
        U.backup_file(output)

    with open(output, 'wt') as opf:
        csvwriter = csv.writer(opf, delimiter='\t')
        csvwriter.writerow(S.HEADER)

        iters = tqdm(enumerate(c2g_bam), desc='processed', unit=' contigs')
        for k, contig in iters:
            if contig.is_unmapped:
                continue

            process_suffix(
                contig, r2c_bam, ref_fa, csvwriter)
            process_bridge_and_link(
                contig, r2c_bam, ref_fa, csvwriter)
            # process_blank(contig, csvwriter)

if __name__ == "__main__":
    main()
