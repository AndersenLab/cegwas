# data files

- 124strain_5perc-MAF_imputed_SNPs.txt : SNPs identified by RADseq (Andersen et al. 2012) data set lifted to WS245 by Dan Cook. SNPs were imputed and only those with greater than 5% minor allele frequency were kept. These are used for mappings. 

- whole_genome_kinship.Rdata : 124 x 124 strain kinship matrix generated using the a.mate function in the rrBLUP package. Whole-genome SNP data was used to generate the relatedness of the strains in this file. 

- fill_error.Rda : empty data frame that resembles emma mapping output that is used for error handling in mapping function.


