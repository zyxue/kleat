import numpy as np
import pandas as pd


KARBOR_FEATURE_COLS = [
    'contig_len',
    'contig_mapq',
    'num_suffix_reads',
    'max_suffix_contig_tail_len',
    'num_suffix_contigs',

    'num_bridge_reads',
    'max_bridge_read_tail_len',
    'max_bridge_read_tail_len',

    'num_link_reads',
    'num_link_contigs',

    'num_blank_contigs',

    'signed_dist_to_aclv',
    # 'abs_dist_to_aclv',

    # PAS hexamer shown on contig
    'ctg_AATAAA', 'ctg_ATTAAA', 'ctg_AGTAAA', 'ctg_TATAAA',
    'ctg_CATAAA', 'ctg_GATAAA', 'ctg_AATATA', 'ctg_AATACA',
    'ctg_AATAGA', 'ctg_AAAAAG', 'ctg_ACTAAA', 'ctg_AAGAAA',
    'ctg_AATGAA', 'ctg_TTTAAA', 'ctg_AAAACA', 'ctg_GGGGCT', 'ctg_NA',

    # # PAS hexamer shown on reference genome
    # 'ref_AATAAA', 'ref_ATTAAA', 'ref_AGTAAA', 'ref_TATAAA',
    # 'ref_CATAAA', 'ref_GATAAA', 'ref_AATATA', 'ref_AATACA',
    # 'ref_AATAGA', 'ref_AAAAAG', 'ref_ACTAAA', 'ref_AAGAAA',
    # 'ref_AATGAA', 'ref_TTTAAA', 'ref_AAAACA', 'ref_GGGGCT', 'ref_NA',
]


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


def compare(df_pred, df_ref, dist_cutoff=50):
    df_mapped = map_clvs(df_pred, df_ref)
    df_mapped['is_tp'] = df_mapped.abs_dist < dist_cutoff

    sensitivity = df_mapped.query('is_tp').shape[0] / df_ref.shape[0]
    precision = df_mapped.query('is_tp').shape[0] / df_mapped.shape[0]
    f1 = (2 * sensitivity * precision) / (sensitivity + precision)

    return sensitivity, precision, f1


def load_polya_seq_df_114genes(sample_id):
    dd = {
        'UHRC1': '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/UHR/C1/polyA-Seq/polyA-Seq-truth-114-genes.csv',
        'UHRC2': '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/UHR/C2/polyA-Seq/polyA-Seq-truth-114-genes.csv',
        'HBRC4': '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/HBR/C4/polyA-Seq/polyA-Seq-truth-114-genes.csv',
        'HBRC6': '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/HBR/C6/polyA-Seq/polyA-Seq-truth-114-genes.csv'
    }
    return pd.read_csv(dd[sample_id])


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
        'HBRC6': '/projects/cheny_prj/KLEAT_benchmarking/polyA_seq/Brain2.bed'
    }
    return load_bed(dd[sample_id])
