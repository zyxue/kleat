Only tested with Python-3.6

```
for i in benchmark_transcriptome/{UHRC[12],HBRC[46]}.tsv ; do python ./benchmark_scripts/prepare_for_ml.py $i; done

python train_arbor.py HBRC4

for i in UHRC1 UHRC2 HBRC4 HBRC6 ; do echo "python ./benchmark_scripts/max_recall_analysis.py $i"; done  | parallel --ungroup
```

