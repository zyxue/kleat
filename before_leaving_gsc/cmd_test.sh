set -o xtrace

DATA_DIR='/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/UHR/C1/tasrkleat-results'
kleat \
    --contigs-to-genome ${DATA_DIR}/align_contigs2genome/cba.sorted.bam \
    --reads-to-contigs ${DATA_DIR}/align_reads2contigs/cba.sorted.bam \
    --reference-genome /gsc/www/bcgsc.ca/downloads/tasrkleat-static/on-cloud/hg19.fa \
    --karbor-clv-annotation ./karbor-annot-clv-hg19.pkl \
    --keep-pre-aggregation-tmp-file \
    --output-format tsv \
    --output ./test-on-114genes-UHRC1.tsv \
    --num-cpus 30

# DATA_DIR=./benchmark_transcriptome/
# sid=UHRC1
# kleat \
#     --contigs-to-genome ${DATA_DIR}/${sid}/c2g.sorted.bam \
#     --reads-to-contigs ${DATA_DIR}/${sid}/r2c_sorted.bam \
#     --reference-genome /gsc/www/bcgsc.ca/downloads/tasrkleat-static/on-cloud/hg19.fa \
#     --karbor-clv-annotation ./karbor-annot-clv-hg19.pkl \
#     --keep-pre-aggregation-tmp-file \
#     --output-format tsv \
#     --output ./test-on-UHRC1.tsv \
#     --num-cpus 10


# # DATA_DIR=/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-PAS-dynamics-analysis/../tcga/BRCA/RNA/RNA-Seq/UNC-LCCC/ILLUMINA/UNCID_2185669.176b2422-3599-41c5-be62-ae6e28d60b90.121205_UNC14-SN744_0276_AC19W1ACXX_6_CAGATC.tar.gz/tasrkleat-results/
# DATA_DIR=/projects/btl/zxue/tasrkleat-TCGA-results/tcga/BRCA/RNA/RNA-Seq/UNC-LCCC/ILLUMINA/UNCID_2169005.7cd7e1cf-4390-4b0e-884c-fc0a26c7e1bf.130606_UNC9-SN296_0372_AC277JACXX_5_TAGCTT.tar.gz/tasrkleat-results/
# kleat \
#         --contigs-to-genome ${DATA_DIR}/align_contigs2genome/cba.sorted.bam \
#         --reads-to-contigs ${DATA_DIR}/align_reads2contigs/cba.sorted.bam \
#         --reference-genome /gsc/www/bcgsc.ca/downloads/tasrkleat-static/on-cloud/hg19.fa \
#         --karbor-clv-annotation ./karbor-annot-clv-hg19.pkl \
#         --output-format tsv \
#         --output ./test-on-tcga-BRCA.tsv \
#         --num-cpus 10


# # DATA_DIR=/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-PAS-dynamics-analysis/../tcga/BRCA/RNA/RNA-Seq/UNC-LCCC/ILLUMINA/UNCID_2185669.176b2422-3599-41c5-be62-ae6e28d60b90.121205_UNC14-SN744_0276_AC19W1ACXX_6_CAGATC.tar.gz/tasrkleat-results/
# # DATA_DIR=/projects/btl/zxue/tasrkleat-TCGA-results/tcga/BRCA/RNA/RNA-Seq/UNC-LCCC/ILLUMINA/UNCID_2169005.7cd7e1cf-4390-4b0e-884c-fc0a26c7e1bf.130606_UNC9-SN296_0372_AC277JACXX_5_TAGCTT.tar.gz/tasrkleat-results/
# DATA_DIR=/projects/btl/zxue/tasrkleat-TCGA-results/tcga/BRCA/RNA/RNA-Seq/UNC-LCCC/ILLUMINA/UNCID_2185124.49a0f1eb-6ece-48fc-970f-3922263fd20b.130124_UNC9-SN296_0328_AD188HACXX_3_TAGCTT.tar.gz/tasrkleat-results/
# kleat \
#     --contigs-to-genome ${DATA_DIR}/align_contigs2genome/cba.sorted.bam \
#     --reads-to-contigs ${DATA_DIR}/align_reads2contigs/cba.sorted.bam \
#     --reference-genome /gsc/www/bcgsc.ca/downloads/tasrkleat-static/on-cloud/hg19.fa \
#     --karbor-clv-annotation ./karbor-annot-clv-hg19.pkl \
#     --output-format tsv \
#     --output ./test-on-tcga-BRCA.tsv \
#     --num-cpus 10
