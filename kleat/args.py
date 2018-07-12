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
        '-a', '--karbor-clv-annotation', type=str, required=True,
        help=('the annotated clv pickle formatted for karbor with '
              '(seqname, strand, clv, gene_ids, gene_names) columns '
              'this file is processed from GTF annotation file')
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
    parser.add_argument(
        '--keep-pre-aggregation-tmp-file', action='store_true',
        help=('specify this if you would like to keep the tmp file before '
              'aggregating polyA evidence per cleavage site, mostly for '
              'debugging purpose')
    )
    return parser.parse_args()
