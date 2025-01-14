# MTClass
Using machine learning to classify eGene-eQTL pairs based on multiple phenotypes. The MTClass algorithm provides a new framework to think about and conduct multivariate genome-wide association.

## Classification of gene-SNP pairs
We adopted an ensemble machine learning a to classify an individual's genotype based on the vector of expression levels from multiple phenotypes (tissues, exons, isoforms, cell types, etc.). In doing so, we hope to uncover eQTLs that have broad effects on gene expression across multiple phenotypes.

## GWAS Colocalization Analysis
One metric for assessing the functionality of our top SNPs is by calculating the colocalization with known GWAS signals. We did this by downloading the entire GWAS Catalog (https://www.ebi.ac.uk/gwas/docs/file-downloads) and searching 10kb upstream and downstream of a given eQTL for known trait associations. We compared our machine learning framework to state-of-the-art linear multivariate approaches, namely MultiPhen and MANOVA, which both output a nominal p-value.