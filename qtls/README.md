# Compute predictor weights, eQTLs, and Bonferroni correction

*For each chromosome, scan over every gene and compute the weights, eQTLs, and Bonferroni correction*

This can be parallelized on a compute cluster

```
for c in `seq 1 22`; do
bash assoc_chr.sh $c
done
```

Internally, each job runs `compute_bonferroni.R` to do the weighted Bonferroni correction.

Finally, we merge all of the chromosome files into one:

```
cat chr*.bonf.tab > all_genes_bonferroni.tab
```
