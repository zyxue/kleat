import glob
import os
import argparse
import csv
import tempfile
import multiprocessing

import pandas as pd
from Bio import Seq, SeqIO


"""
This script builds a comprehensive catalog of polyadenylation signal for a
given reference genome
"""


def rev_comp(s):
    return str(Seq.Seq(s).reverse_complement())


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
    (rev_comp(i), j) for (i, j) in CANDIDATE_HEXAMERS
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--input-fa', type=str,
        help='e.g. a fasta file of a reference genome'
    )
    parser.add_argument(
        '-o', '--output-pkl', type=str,
        help='output pickle file'
    )
    return parser.parse_args()


def gen_output_per_rec(rec):
    seqname = rec.id
    return f'./tmp/{seqname}.csv'


def build_catalog(rec):
    """build catalog for one chromosome"""
    seq = rec.seq.upper()
    seqname = rec.id

    output = gen_output_per_rec(rec)
    with open(output, 'wt') as opf:
        csvwriter = csv.writer(opf)
        csvwriter.writerow([
            'seqname', 'strand', 'position', 'hexamer', 'hexamer_id'
        ])

        print(f'working on {output}')
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

                    # reverse back
                    hmr = rev_comp(hexamer) if strand == '-' else hexamer
                    csvwriter.writerow([
                        seqname, strand, index, hmr, hexamer_id
                    ])
                    index += 1
    print(f'{output} finished.')


if __name__ == "__main__":
    args = parse_args()

    io = SeqIO.parse(args.input_fa, 'fasta')

    if not os.path.exists('./tmp'):
        os.mkdir('./tmp')

    recs = []
    num_recs = 0
    for rec in io:
        num_recs += 1

        output = gen_output_per_rec(rec)
        if os.path.exists(output):
            print(f'{output} exists.')
            continue

        print(f'adding {rec.id} to the list')
        recs.append(rec)

    if len(recs) > 0:
        num_procs = min(len(recs), multiprocessing.cpu_count())
        pool = multiprocessing.Pool(num_procs)
        pool.map(build_catalog, recs)

    if num_recs > 1:
        print(f'prepare to concatenate all chromosome files into one...')
        dfs = []
        for f in glob.glob('./tmp/*.csv'):
            print(f'reading {f}')
            _df = pd.read_csv(f)
            dfs.append(_df)
        print(f'concatenating...')
    output_df = pd.concat(dfs)
    output_df.to_pickle(args.output_pkl)
    print(f'concatenated into f{output_pkl}')
