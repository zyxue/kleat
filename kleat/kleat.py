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


def write_row(clv_record, csvwriter):
    csvwriter.writerow([getattr(clv_record, _) for _ in S.HEADER])


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
        '-o', '--output', type=str, default='./output.csv',
        help='output csv file'
    )
    return parser.parse_args()


def process_suffix(contig, r2c_bam, csvwriter):
    clv_record = suffix.gen_clv_record(contig, r2c_bam)
    write_row(clv_record, csvwriter)


def process_bridge_and_link(contig, r2c_bam, csvwriter):
    # bridge & link
    aligned_reads = r2c_bam.fetch(contig.query_name)
    dd_bridge, dd_link = extract_bridge_and_link(contig, aligned_reads)
    write_bridge(dd_bridge, contig, csvwriter)
    write_link(dd_link, contig, csvwriter)


def process_blank(contig, csvwriter):
    for clv_rec in blank.gen_two_clv_records(contig):
        write_row(clv_rec, csvwriter)


def init_bridge_evidence_holder():
    """
    initialize holders for bridge and link evidence of a given contig
    """
    return {
        'num_reads': defaultdict(int),
        'max_tail_len': defaultdict(int),
    }


def init_link_evidence_holder():
    return {
        'num_reads': defaultdict(int)
    }


def is_bridge_read(read):
    return not read.is_unmapped and apautils.has_tail(read)


def update_bridge(evid_tuple, evid_holder):
    """
    update the bridge evidence holder with extracted bridge evidence

    :param evid_tuple: evidence tuple
    :param evid_holder: A dict holder bridge evidence for a given contig
    """
    seqname, strand, ref_clv, tail_len = evid_tuple
    clv_key = apautils.gen_clv_key_tuple(seqname, strand, ref_clv)
    evid_holder['num_reads'][clv_key] += 1
    evid_holder['max_tail_len'][clv_key] = max(
        evid_holder['max_tail_len'][clv_key], tail_len)


def is_link_read(read):
    # Here we focused on the unmapped all A/T read, but in
    # principle, we could also check from the perspecitve and
    # the mate of a link read, but it would be harder to verify
    # the sequence composition of the link read
    return (not read.mate_is_unmapped
            # assume reference_id comparison is faster than
            # reference_name
            and read.reference_id == read.next_reference_id
            and set(read.query_sequence) in [{'A'}, {'T'}])


def update_link(evid_tuple, evid_holder):
    seqname, strand, ref_clv = evid_tuple
    clv_key = apautils.gen_clv_key_tuple(seqname, strand, ref_clv)
    evid_holder['num_reads'][clv_key] += 1


def extract_bridge_and_link(contig, aligned_reads):
    """
    bridge and link are processed together by looping throught reads aligned to
    the contig

    :param contig: a pysam.libcalignedsegment.AlignedSegment instance
    :param aligned_reads: a pysam.libcalignmentfile.IteratorRowRegion instance
    """
    dd_bridge = init_bridge_evidence_holder()
    dd_link = init_link_evidence_holder()
    for read in aligned_reads:
        if is_bridge_read(read):
            bdg_evid = bridge.analyze_bridge_read(contig, read)
            update_bridge(bdg_evid, dd_bridge)
        elif is_link_read(read):
            link_evid = link.analyze_link(contig, read)
            update_link(link_evid, dd_link)
    return dd_bridge, dd_link


def write_bridge(dd_bridge, contig, csvwriter):
    if len(dd_bridge['num_reads']) == 0:
        return
    for clv_key in dd_bridge['num_reads']:
        clv_record = bridge.gen_clv_record(
            contig, clv_key,
            dd_bridge['num_reads'][clv_key],
            dd_bridge['max_tail_len'][clv_key]
        )
        write_row(clv_record, csvwriter)


def write_link(dd_link, contig, csvwriter):
    if len(dd_link['num_reads']) == 0:
        return
    for clv_key in dd_link['num_reads']:
        clv_record = link.gen_clv_record(
            contig, clv_key, dd_link['num_reads'][clv_key])
        write_row(clv_record, csvwriter)


def main():
    args = get_args()
    c2g_bam = pysam.AlignmentFile(args.contig_to_genome)
    r2c_bam = pysam.AlignmentFile(args.read_to_contig)
    output = args.output
    if os.path.exists(output):
        U.backup_file(output)

    with open(output, 'wt') as opf:
        csvwriter = csv.writer(opf)
        csvwriter.writerow(S.HEADER)

        for k, contig in tqdm(enumerate(c2g_bam)):
            if contig.is_unmapped:
                continue

            if apautils.has_tail(contig):
                process_suffix(contig, r2c_bam, csvwriter)
            process_bridge_and_link(contig, r2c_bam, csvwriter)
            process_blank(contig, csvwriter)

if __name__ == "__main__":
    main()
