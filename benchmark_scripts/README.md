```
for i in benchmark_transcriptome/{UHRC[12],HBRC[46]}.tsv ; do python prepare_for_ml.py $i; done

python train_arbor.py HBRC4
```
