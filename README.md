# Genomic Deep Learning Benchmarks

Working off the analyses in [Huang et al.](https://www.nature.com/articles/s41588-023-01574-w) to further compare genomic DL models to predixcan/QTL and naive promoter SNP benchmarks. This analysis uses the public GEUVADIS / 1000 Genomes data.

Currently there are two analyses:
* A simple re-analysis of the Huang et al. data to compare prediction R2 with that of a promoter SNS. This analysis is simple and can be easily reproduced and extended with a single R script if a new predicted expression matrix is available.
* eQTL discovery followed by weighted FDR to identify eGenes. This analysis is a bit more intensive to extend with new prediction weights because SNP weights have to be recomputed for every gene.

## Prediction R2

This is a basic re-analysis of the Huang et al. data to also include a promoter SNP as baseline. It can be reproduced by calling:

```
R --slave < analysis_rsq.R
```

Which will load in the predictions from all methods (including our bespoke promoter SNP annotation), compute the adjusted R2 for all genes, and generate the figures `rsq_scatterplot` and `rsq_barplot`. Additional predictive models can be added to the `method_files` and `method_names` variables.

## eQTL Re-weighting

This analysis estimates eQTLs, computes weights for the genomically (or Predixcan) predicted expression, and performs weighted Bonferroni correction. Details on the analysis (which is split up by chromosome) are in the `qtls` README. The complete analysis produces the `qtls/all_genes_bonferroni.tab` file.

Then the following command will generate the corresponding barplots ( `qvalue` package assumed to be installed in R):

```
R --slave < analysis_qtl.R
```

Which will generate the `plot_eqtl.pdf` figure.