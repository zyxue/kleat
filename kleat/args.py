import argparse


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
        '-f', '--reference-genome', type=str, required=True,
        help=('reference genome FASTA file, if provided, '
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
    parser.add_argument(
        '-p', '--num-cpus', type=int, default=1,
        help=('parallize the step of aggregating polya evidence for each '
              '(seqname, strand, clv)')
    )
    return parser.parse_args()
