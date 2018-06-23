import os
import glob


from compare_to_polyA_Seq import (
    load_kleat3_df, load_polya_df, replace_seqname, compare
)


for pred_file in glob.glob('./cv/kleat3_*.csv'):
    bn = os.path.basename(pred_file)
    sample_id = pred_file.split('_')[1]
    tree_max_depth = int(pred_file.split('.')[-2].split('h')[-1])  # dirty hack!
    # print(f'working on {sample_id}, max_depth: {max_depth}')

    './pkl/kleat3_HBRC4_ml_filtered_clustered_max_depth2.csv'

    rootdir = '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk'
    if sample_id == 'UHRC1':
        truth_csv = f'{rootdir}/UHR/C1/polyA-Seq/polyA-Seq-truth-114-genes.csv'
    elif sample_id == 'UHRC2':
        truth_csv = f'{rootdir}/UHR/C2/polyA-Seq/polyA-Seq-truth-114-genes.csv'
    elif sample_id == 'HBRC4':
        truth_csv = f'{rootdir}/HBR/C4/polyA-Seq/polyA-Seq-truth-114-genes.csv'
    elif sample_id == 'HBRC6':
        truth_csv = f'{rootdir}/HBR/C6/polyA-Seq/polyA-Seq-truth-114-genes.csv'
    else:
        raise

    df_ref = load_polya_df(truth_csv)
    df_pred = load_kleat3_df(pred_file)

    replace_seqname(df_ref)
    replace_seqname(df_pred)

    se, pr, f1 = compare(df_pred, df_ref)
    print(f'{sample_id}\t{tree_max_depth}\t{pr}\t{se}\t{f1}')
