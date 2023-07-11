# MTClass
Using machine learning to classify gene-SNP pairs

## Classification of gene-SNP pairs
We adopted an ensemble learning framework to classify an individual's genotype based on the vector of expression levels from multiple tissues or multiple exons. In doing so, we hope to uncover eQTLs that have effects on gene expression across multiple phenotypes.

The general syntax of the script is as follows:
```
python3 classify.py exp gt out_dir
```
where:
  * exp = path to gene expression file (in .txt or .txt.gz format)
  * gt = path to genotypes file (in .txt or .txt.gz format)
  * out_dir = output directory to store results

### Gene expression file
The gene expression file must be structured in the following way, with the first two columns being ```gene``` and ```donor```. The rest of the columns are the feature names:
| gene | donor | pheno1 | pheno2 | pheno3 | ... |
| --- | --- | --- | --- | --- | --- |
| HBB | Sample1 | 8.19 | 7.12 | 11.47 | ... |
| HBB | Sample2 | 5.01 | 12.70 | 2.15 | ... |

### Genotypes file
The genotypes file must be structured in the following way, with the first two columns being ```gene``` and ```ID```. The rest of the columns are the sample names. The sample names must match those in the ```donor``` column of the gene expression file.

**Note**: The formatting of this genotype file is very similar to that of a typical VCF, except with the addition of the ```gene``` column and the binarization of the genotypes.
| gene | ID | Sample1 | Sample2 | ... |
| --- | --- | --- | --- | --- |
| HBB | chr11_12581527_A_T_b38 | 0 | 1 | ... |
| HBB | chr11_12592567_G_C_b38 | 1 | 1 | ... |

### Example usage
We provide an example of using the MTClass classification script with GTEx (https://gtexportal.org) expression levels of the HBB gene from 9 tissues. To preserve subject confidentiality, however, the true genotypes of the GTEx donors were randomized, so classification performance on this example dataset is poor for all eQTLs.
``` 
python3 classify.py ExampleData/test_exp.txt ExampleData/test_gt.txt ./
```

## GWAS Colocalization Analysis
One metric of assessing the functionality of our top SNPs is by calculating the colocalization with known GWAS signals. We did this by downloading the entire GWAS Catalog (https://www.ebi.ac.uk/gwas/docs/file-downloads) and searching 10kb upstream and downstream of a given SNP for known GWAS associations. We compared our machine learning framework to state-of-the-art linear multivariate approaches, namely MultiPhen and MANOVA, which both output a nominal p-value. We have created a script and sampled 100,000 gene-SNP pairs from our results to serve as an example of the utility of this method.
