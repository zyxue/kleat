Only tested with Python-3.6

```
for i in benchmark_transcriptome/{UHRC[12],HBRC[46]}.tsv ; do python ./benchmark_scripts/prepare_for_ml.py $i; done

python train_arbor.py HBRC4
```
