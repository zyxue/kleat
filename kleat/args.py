import argparse


def get_args():
    parser = argparse.ArgumentParser(
        description='KLEAT: cleavage site detection via de novo assembly')
    parser.add_argument(
        '-c', '--contigs-to-genome', type=str, required=True,
        help='input contig-to-genome alignment BAM file'
    )
    parser.add_argument(
        '-r', '--reads-to-contigs', type=str, required=True,
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
        '-o', '--output', type=str, default=None,
        help=('output tsv file, if not specified, it will use prefix output, '
              'and the extension depends on the value of --output-format. '
              'e.g. output.csv, output.pickle, etc.')
    )

    parser.add_argument(
        '-m', '--output-format', type=str, default='csv',
        help='also support tsv, pickle (python)'
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

    parser.add_argument(
        '--bridge-skip-check-size', type=int, default=3,
        help=('the size beyond which the clv is predicted be on the next '
              'matched region. Otherwise, clv is predicted to be at the edge '
              'of the prevision match. The argument is added because '
              'inconsistent break points are observed between read '
              '(softclip as in r2c alignment) and contig '
              '(boundry between BAM_CMATCH and BAM_CREF_SKIP)')
    )

    return parser.parse_args()
