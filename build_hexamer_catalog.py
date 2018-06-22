import argparse
import csv
import sys

from Bio import Seq, SeqIO


"""
This script builds a comprehensive catalog of polyadenylation signal for a
given reference genome
"""


# numbers indicate their strengths, 1 is the strongest
CANDIDATE_HEXAMERS = [
    ('AATAAA', 1),
    ('ATTAAA', 2),
    ('AGTAAA', 3),
    ('TATAAA', 4),
    ('CATAAA', 5),
    ('GATAAA', 6),
    ('AATATA', 7),
    ('AATACA', 8),
    ('AATAGA', 9),
    ('AAAAAG', 10),
    ('ACTAAA', 11),
    ('AAGAAA', 12),
    ('AATGAA', 13),
    ('TTTAAA', 14),
    ('AAAACA', 15),
    ('GGGGCT', 16)
]


CANDIDATE_HEXAMERS_REVERSE_COMPLEMENTED = [
    (str(Seq.Seq(i).reverse_complement()), j) for (i, j) in CANDIDATE_HEXAMERS
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input-fa', type=str,
        help='e.g. a fasta file of a reference genome'
    )
    parser.add_argument(
        '-o', '--output-csv', type=str,
        help='output csv file'
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    io = SeqIO.parse(args.input_fa, 'fasta')
    with open(args.output_csv, 'wt') as opf:
        csvwriter = csv.writer(opf)
        csvwriter.writerow([
            'seqname', 'strand', 'position', 'hexamer', 'hexamer_id'
        ])
        for rec in io:
            seq = rec.seq.upper()
            seqname = rec.id
            print(f'working on {seqname}')
            for strand in ['+', '-']:
                if strand == '+':
                    hexamer_list = CANDIDATE_HEXAMERS
                else:
                    hexamer_list = CANDIDATE_HEXAMERS_REVERSE_COMPLEMENTED

                for (hexamer, hexamer_id) in hexamer_list:
                    index = 0
                    while True:
                        index = seq.find(hexamer, start=index)
                        if index == -1:
                            break
                        csvwriter.writerow([
                            seqname, strand, index, hexamer, hexamer_id
                        ])
                        index += 1
