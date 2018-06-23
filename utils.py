import pysam

REF_FA = '/gsc/www/bcgsc.ca/downloads/zxue/tasrkleat-static/off-cloud/reference_data/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa'
REF_SEQ = pysam.FastaFile(REF_FA)


# used when merging polyA evidence
COLS_TO_SUM = [
    'num_contig_tail_reads', 
    'num_bridge_reads', 
    'num_link_reads', 'num_blank_contigs'
]

COLS_TO_MAX = ['contig_len', 'contig_mapq', 'contig_tail_length', 'max_bridge_tail_length']

FEATURE_COLS = [
    'num_contig_tail_reads', 'num_bridge_reads', 'num_link_reads',
    'num_blank_contigs',
    'contig_len',
    'contig_mapq', 'contig_tail_length', 'max_bridge_tail_length',
    'abs_dist_to_aclv'] + [
        'AATAAA', 'ATTAAA', 'AGTAAA', 'TATAAA',
        'CATAAA', 'GATAAA', 'AATATA', 'AATACA',
        'AATAGA', 'AAAAAG', 'ACTAAA', 'AAGAAA',
        'AATGAA', 'TTTAAA', 'AAAACA', 'GGGGCT', 'NA'
    ]




# CLV_SC_MAPPING_CSV = '/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/data/reference_data/annotated-clv-sc-mapping.csv.gz'
# _df = pd.read_csv(CLV_SC_MAPPING_CSV)
# _sr = df_clv_sc.groupby(['seqname', 'strand']).apply(lambda g: g.clv.sort_values().values)

# def calc_abs_dist_to_annot_clv(grp):
# #     print(grp.name)
#     if grp.name in annot_clvs.index:
#          # grp.name holds the group key
#         aclvs = annot_clvs.loc[grp.name]
#         aclvs = annot_clvs.loc[grp.name] # grp.name holds the group key
#         bcast = np.broadcast_to(grp.clv.values, (aclvs.shape[0], grp.shape[0])).T
#         grp['abs_dist_to_aclv'] = np.min(np.abs(bcast - aclvs), axis=1)
#         return grp
#     else:
#         # return None, this group will be gone, which is good 
#         # as they won't be clv of targeted genes for sure, increasing precision
#         return 
