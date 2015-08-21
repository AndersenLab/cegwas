# cegwas

## Installation

**NOTE**

This package requires the installation of the `biomaRt` package from bioconducter::

```r
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
```

Install package:

```r
devtools::install_github("AndersenLab/cegwas")
```

A set of functions to process phenotype data, perform GWAS, and perform post-mapping data processing for C. elegans.

##### The pipeline is split into three steps
1. Process phenotypes
2. Perform GWAS mappings
3. Process mapping data frame

##### Process Phenotypes
Input data frame for this step contains properly formatted phenotype data. The first column should be named *trait* all additional columns should be strains. One row corresponding to one trait for all strains.

_**Example Usage**_

```r
processed_phenotypes <- process_pheno(data)
``` 

This function outputs a list object. Outputs a list. The first element of the list is an ordered vector of traits. The second element of the list is a dataframe containing one column for each strain, with values corresponding to traits in element 1 for rows.

##### GWAS Mappings
Input data for this step is the output from the `process_pheno` function. GWAS mappings are performed using the `GWAS` function from the `rrBLUP` package with a 5% minor allele frequency cutoff for SNPs. Additional input data for this function are built into the package (SNP set & kinship matrix)

_**Example Usage**_

```r
mapping_df <- gwas_mappings(processed_phenotypes, cores = 4, only_sig = TRUE)
``` 

The output for function is a data frame that contains SNP information, trait information, and log transformed p-values.

##### Process Mappings
The input data sets for this step are:

1. The output from the mapping step
2. The processed phenotype data set
3. SNP set

_**Example Usage**_

```r
processed_mapping_df <- process_mappings(mapping_df, snp_df = snps, processed_phenotypes, CI_size = 50, snp_grouping = 200)
```

The resulting dataframe contains all information output from the `gwas_mappings` function as well as 

1. Variance explained calculations for SNPs above the bonferroni corrected p-value
2. Confidence intervals information for all identified peaks

**NOTE** the process_mappings function is also broken up into three separate functions, which have their own documentation:

1. `calculate_VE`
2. `find_peaks`
3. `identify_CI`