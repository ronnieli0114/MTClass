# MTClass
Using machine learning to classify gene-SNP pairs based on multiple phenotypes. A new way to think about multivariate genome-wide association.

## Classification of gene-SNP pairs
We adopted an ensemble learning framework to classify an individual's genotype based on the vector of expression levels from multiple tissues or multiple exons. In doing so, we hope to uncover eQTLs that have effects on gene expression across multiple phenotypes.

The MTClass script utilizes the following dependencies:
* numpy >= 1.25.0
* pandas >= 2.0.3
* scikit-learn >= 1.3.0

The general syntax of the main MTClass script is as follows:
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
One metric for assessing the functionality of our top SNPs is by calculating the colocalization with known GWAS signals. We did this by downloading the entire GWAS Catalog (https://www.ebi.ac.uk/gwas/docs/file-downloads) and searching 10kb upstream and downstream of a given SNP for known GWAS associations. We compared our machine learning framework to state-of-the-art linear multivariate approaches, namely MultiPhen and MANOVA, which both output a nominal p-value. We have created a script and sampled 100,000 gene-SNP pairs from our results to serve as an example of the utility of this method.

The GWAS colocalization analysis script utilizes the following dependencies:
* matplotlib >= 3.7.2
* numpy >= 1.25.0
* pandas >= 2.0.3
* requests >= 2.31.0

The general syntax of the GWAS colocalization analysis script is as follows:
```
python3 gwas_colocalization.py mtclass mtclass_metric multiphen manova out_dir
```

This script requires that you __EITHER__ specify where the downloaded GWAS Catalog is with ```--gwas_catalog ${PATH_TO_GWAS_CATALOG}```, __OR__ you can choose to download the GWAS Catalog into the current working directory with ```--download```. The NHGRI-EBI GWAS Catalog is about 270 MB in size.

The MTClass results must be formatted in the following way:
| gene | variant | Metric1 | Metric2 | ... |
| --- | --- | --- | --- | --- |
| HBB | chr11_12581527_A_T_b38 | 0.751 | 0.732 | ... |

The MultiPhen and MANOVA results are formatted similarly, with ```gene``` and ```variant``` columns, and with a ```pval``` column afterwards instead of classification metrics.

An example dataset using a random subset of 500,000 gene-SNP pairs has been provided in the code. To run this analysis, simply use the following command, which will download the GWAS Catalog, and save a plot and a table in the current working directory:
```
python3 gwas_colocalization.py ExampleResults/mtclass.txt.gz f1_macro_median ExampleResults/multiphen.txt.gz ExampleResults/manova.txt.gz ./ --download
```
