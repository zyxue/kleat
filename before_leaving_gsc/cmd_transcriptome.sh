# set -o xtrace

# # transcriptome

# SAMPLE_IDS="UHRC1 UHRC2 HBRC4 HBRC6"
# # SAMPLE_IDS="UHRC2 HBRC4"
# for sid in ${SAMPLE_IDS}; do
#     DATA_DIR=./benchmark_transcriptome

#     cmd="kleat \
#             --contigs-to-genome ${DATA_DIR}/${sid}/c2g.sorted.bam \
#             --reads-to-contigs ${DATA_DIR}/${sid}/r2c_sorted.bam \
#             --reference-genome /gsc/www/bcgsc.ca/downloads/tasrkleat-static/on-cloud/hg19.fa \
#             --karbor-clv-annotation ./karbor-annot-clv-hg19.pkl \
#             --output ${DATA_DIR}/ctg-hex-comparison/${sid}.pkl \
#             --keep-pre-aggregation-tmp-file \
#             --output-format pickle \
#             --num-cpus 30"
#     echo ${cmd}
# done


# 114 genes

KEYS="UHR/C1 UHR/C2 HBR/C4 HBR/C6"
for key in ${KEYS}; do
    sid=${key/\/}
    DATA_DIR=/projects/btl/zxue/tasrkleat-TCGA-results/tasrkleat-TCGA-analysis-scripts/benchmark-kleat.bk/${key}/tasrkleat-results

    cmd="kleat \
        --contigs-to-genome ${DATA_DIR}/align_contigs2genome/cba.sorted.bam \
        --reads-to-contigs ${DATA_DIR}/align_reads2contigs/cba.sorted.bam \
        --reference-genome /gsc/www/bcgsc.ca/downloads/tasrkleat-static/on-cloud/hg19.fa \
        --karbor-clv-annotation ./karbor-annot-clv-hg19.pkl \
        --keep-pre-aggregation-tmp-file \
        --output-format pickle \
        --num-cpus 20 \
        --output ./benchmark_114genes/tcga-run-4/${sid}.pkl"
    echo ${cmd}
done
