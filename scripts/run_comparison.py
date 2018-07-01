import sys
import os
import glob


from compare_to_polyA_Seq import (
    # after ml filter, load_kleat3_df applies to kleat2 results, too
    load_kleat3_df, load_polya_df, load_bed, replace_seqname, compare
)

kleat_version = sys.argv[1]
train_sample = sys.argv[2]

glob_pattern = f'./{kleat_version}_ml/cv/trained_on_{train_sample}/max_depth*/*.csv'
# print(glob_pattern)
pred_files = glob.glob(glob_pattern)
for pred_file in pred_files:
    bn = os.path.basename(pred_file)
    sample_id = bn.split('.')[0]
    tree_max_depth = int(os.path.basename(os.path.dirname(pred_file)).split('h')[-1])  # dirty hack!

    rootdir = '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk'
    if sample_id in ['UHRC1', 'UHRC2', 'HBRC4', 'HBRC6']:
        truth_csv = f'{rootdir}/{sample_id[:3]}/{sample_id[3:]}/polyA-Seq/polyA-Seq-truth-114-genes.csv'
        df_ref = load_polya_df(truth_csv)
    elif sample_id == 'Brain-C4_S3_RNABloom':
        bed = '/projects/cheny_prj/KLEAT_benchmarking/polyA_seq/Brain1.bed'
        df_ref = load_bed(bed)
    elif sample_id == 'Brain-C6_S4_RNABloom':
        bed = '/projects/cheny_prj/KLEAT_benchmarking/polyA_seq/Brain2.bed'
        df_ref = load_bed(bed)
    elif sample_id == 'UHRR-C1_S1_RNABloom':
        bed = '/projects/cheny_prj/KLEAT_benchmarking/polyA_seq/UHR1.bed'
        df_ref = load_bed(bed)
    elif sample_id == 'UHRR-C2_S2_RNABloom':
        bed = '/projects/cheny_prj/KLEAT_benchmarking/polyA_seq/UHR2.bed'
        df_ref = load_bed(bed)
    else:
        raise

    df_pred = load_kleat3_df(pred_file)

    replace_seqname(df_ref)
    replace_seqname(df_pred)

    se, pr, f1 = compare(df_pred, df_ref)
    print(f'{sample_id}\t{tree_max_depth}\t{pr}\t{se}\t{f1}')
