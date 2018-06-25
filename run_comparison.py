import sys
import os
import glob


from compare_to_polyA_Seq import (
    # after ml filter, load_kleat3_df applies to kleat2 results, too
    load_kleat3_df, load_polya_df, replace_seqname, compare
)

kleat_version = sys.argv[1]

pred_files = glob.glob(f'./{kleat_version}_ml/cv/max_depth*/*.csv')
for pred_file in pred_files:
    bn = os.path.basename(pred_file)
    sample_id = bn.split('.')[0]
    tree_max_depth = int(os.path.basename(os.path.dirname(pred_file)).split('h')[-1])  # dirty hack!

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
        raise ValueError(f'unknown smaple_id: {sample_id}')

    df_ref = load_polya_df(truth_csv)
    df_pred = load_kleat3_df(pred_file)

    replace_seqname(df_ref)
    replace_seqname(df_pred)

    se, pr, f1 = compare(df_pred, df_ref)
    print(f'{sample_id}\t{tree_max_depth}\t{pr}\t{se}\t{f1}')
