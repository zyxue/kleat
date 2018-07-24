import logging

import numpy as np
import pandas as pd

from kleat.misc.settings import CANDIDATE_HEXAMERS_WITH_NA


CTG_HEX_DUMMY_COLS = ['ctg_{0}'.format(_[0]) for _ in CANDIDATE_HEXAMERS_WITH_NA]
REF_HEX_DUMMY_COLS = ['ref_{0}'.format(_[0]) for _ in CANDIDATE_HEXAMERS_WITH_NA]

KARBOR_FEATURE_COLS = (
    [
        'signed_dist_to_aclv',
        # 'abs_dist_to_aclv',

        'max_contig_len',
        'max_contig_mapq',
        'any_contig_is_hardclipped',

        'num_suffix_contigs',
        'num_bridge_contigs',
        'num_link_contigs',
        'num_blank_contigs',
        'num_total_contigs',

        'num_suffix_reads',
        'max_suffix_contig_tail_len',

        'num_bridge_reads',
        'max_bridge_read_tail_len',

        'num_link_reads',

        'ctg_hex_dist',
        # 'ref_hex_dist',
    ]
    + CTG_HEX_DUMMY_COLS        # PAS hexamer shown on contig
    # + REF_HEX_DUMMY_COLS        # PAS hexamer shown on reference genome
)


def map_clvs(df_pred, df_ref):
    ref_dd = df_ref.groupby(['seqname', 'strand']).apply(
        lambda g: g.clv.values).to_dict()

    _dfs = []
    for k, g in df_pred.groupby(['seqname', 'strand']):
        _df = g.copy()

        ref_clvs = ref_dd.get(k)
        if ref_clvs is None:
            _df['mapped_ref_clv'] = np.nan
        else:
            _df['mapped_ref_clv'] = g.clv.apply(
                lambda v: ref_clvs[np.argmin(np.abs(v - ref_clvs))])
        _dfs.append(_df)

    df_mapped = pd.concat(_dfs)

    df_mapped['dist'] = df_mapped.mapped_ref_clv - df_mapped.clv
    df_mapped['abs_dist'] = np.abs(df_mapped['dist'])
    return df_mapped


def compare(df_pred, df_ref, map_cutoff):
    """
    :param dist_cutoff: below which the predicted clv and the annotated clv are
    considered the same
    """
    df_mapped = map_clvs(df_pred, df_ref)
    df_mapped['is_tp'] = df_mapped.abs_dist < map_cutoff

    sensitivity = df_mapped.query('is_tp').shape[0] / df_ref.shape[0]
    precision = df_mapped.query('is_tp').shape[0] / df_mapped.shape[0]
    f1 = (2 * sensitivity * precision) / (sensitivity + precision)

    return sensitivity, precision, f1


def load_bed(input_bed):
    """"mostly for whole transcriptome"""
    return pd.read_csv(
        input_bed, header=None, sep='\t',
        names=['seqname', 'clv', 'clv1', 'hexamer_str', 'unknown', 'strand'])


def load_polya_seq_df(sample_id):
    dd = {
        'UHRC1': '/projects/cheny_prj/KLEAT_benchmarking/polyA_seq/UHR1.bed',
        'UHRC2': '/projects/cheny_prj/KLEAT_benchmarking/polyA_seq/UHR2.bed',
        'HBRC4': '/projects/cheny_prj/KLEAT_benchmarking/polyA_seq/Brain1.bed',
        'HBRC6': '/projects/cheny_prj/KLEAT_benchmarking/polyA_seq/Brain2.bed',

        'UHRC1_114genes': '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/UHR/C1/polyA-Seq/polyA-Seq-truth-114-genes.csv',
        'UHRC2_114genes': '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/UHR/C2/polyA-Seq/polyA-Seq-truth-114-genes.csv',
        'HBRC4_114genes': '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/HBR/C4/polyA-Seq/polyA-Seq-truth-114-genes.csv',
        'HBRC6_114genes': '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/HBR/C6/polyA-Seq/polyA-Seq-truth-114-genes.csv'
    }

    infile = dd[sample_id]

    logging.info(f'reading {infile} ...')
    if infile.endswith('csv'):
        df = pd.read_csv(infile)
        logging.info('reading done'.format(infile))
    elif infile.endswith('bed'):
        df = load_bed(infile)
    else:
        raise ValueError('unknown file extension: {0}'.format(infile))
    return df
