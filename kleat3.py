#!/usr/bin/env python

import logging
import csv
import argparse
from collections import defaultdict

import pysam
from tqdm import tqdm

from evidence import suffix, bridge, link, blank
from settings import HEADER
import utils as U


logging.basicConfig(
    level=logging.DEBUG, format='%(asctime)s|%(levelname)s|%(message)s')

logger = logging.getLogger(__name__)


def write_row(clv_record, csvwriter):
    csvwriter.writerow([getattr(clv_record, _) for _ in HEADER])


def get_args():
    parser = argparse.ArgumentParser(
        description='KLEAT3: cleavage site detection via de novo assembly')
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


if __name__ == "__main__":
    args = get_args()

    c2g_bam = pysam.AlignmentFile(args.contig_to_genome)
    r2c_bam = pysam.AlignmentFile(args.read_to_contig)
    output = args.output
    U.backup_file(output)

    with open(output, 'wt') as opf:
        csvwriter = csv.writer(opf)
        csvwriter.writerow(HEADER)

        for k, contig in tqdm(enumerate(c2g_bam)):
            if contig.is_unmapped:
                continue

            # suffix evidence
            if U.has_tail(contig):
                clv_record = suffix.gen_clv_record(contig, r2c_bam)
                write_row(clv_record, csvwriter)

            dd_bridge = {
                'num_reads': defaultdict(int),
                'max_tail_len': defaultdict(int),
            }

            dd_link = {
                'num_reads': defaultdict(int)
            }

            for read in r2c_bam.fetch(contig.query_name):
                # bridge evidence
                if not read.is_unmapped:
                    if U.has_tail(read):
                        seqname, strand, ref_clv, tail_len = \
                            bridge.analyze_bridge_read(contig, read)

                        clv_key = U.gen_clv_key_tuple(seqname, strand, ref_clv)
                        dd_bridge['num_reads'][clv_key] += 1
                        dd_bridge['max_tail_len'][clv_key] = max(
                            dd_bridge['max_tail_len'][clv_key], tail_len)
                # link evidence
                else:
                    # Here we focused on the unmapped all A/T read, but in
                    # principle, we could also check from the perspecitve and
                    # the mate of a link read, but it would be harder to verify
                    # the sequence composition of the link read
                    if (not read.mate_is_unmapped
                        # assume reference_id comparison is faster than
                        # reference_name
                        and read.reference_id == read.next_reference_id
                        and set(read.query_sequence) in [{'A'}, {'T'}]):
                        seqname, strand, ref_clv = \
                            link.analyze_link(contig, read)

                        clv_key = U.gen_clv_key_tuple(seqname, strand, ref_clv)
                        dd_link['num_reads'][clv_key] += 1

            if len(dd_bridge['num_reads']) > 0:
                for clv_key in dd_bridge['num_reads']:
                    clv_record = bridge.gen_clv_record(
                        contig, clv_key,
                        dd_bridge['num_reads'][clv_key],
                        dd_bridge['max_tail_len'][clv_key]
                    )
                    write_row(clv_record, csvwriter)

            if len(dd_link['num_reads']) > 0:
                for ref_clv in dd_link['num_reads']:
                    clv_record = link.gen_clv_record(
                        contig, clv_key, dd_link['num_reads'][clv_key])
                    write_row(clv_record, csvwriter)

            # blank evidence
            for clv_rec in blank.gen_two_clv_records(contig):
                write_row(clv_rec, csvwriter)
