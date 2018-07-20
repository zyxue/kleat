Only tested with Python-3.6

```
for i in benchmark_transcriptome/{UHRC[12],HBRC[46]}.pkl ; do python ./benchmark_scripts/prepare_for_ml.py $i; done

time python benchmark_scripts/train_arbor.py --train-sample-id HBRC4 --test-sample-ids HBRC4 HBRC6 UHRC1 UHRC2 --max-depths 1 40 1 --pickle-depth 6 8 --output bm_output.csv --num-cpus 20

for i in UHRC1 UHRC2 HBRC4 HBRC6 ; do echo "python ./benchmark_scripts/max_recall_analysis.py $i"; done  | parallel --ungroup
```

python benchmark_scripts/predict.py --test-sample-ids HBRC4_114genes HBRC6_114genes UHRC1_114genes UHRC2_114genes --classifier-pkl benchmark_transcriptome/tcga-run-3/benchmark.tree_max_depth_7.pkl --output benchmark_114genes/tcga-run-3/benchmark.csv
