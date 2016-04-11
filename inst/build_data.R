# This script processes data.

# gene Ids
library("readr")
library("dplyr")
library("devtools")

# Download
# url <- "ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/c_elegans.PRJNA13758.current.geneIDs.txt.gz"
# tmp <- tempfile()
# download.file(url,tmp)
# gene_info <- readr::read_csv(gzfile(tmp), col_names = c("taxon","WBID","name","protein","live")) %>%
#   dplyr::select(-taxon)
# 
# gene_ids <- as.list(gene_info$name)
# names(gene_ids) <- gene_info$WBID
# 
# save(gene_info, file = "data/gene_info.rda")
# save(gene_ids, file = "data/gene_ids.rda")

load("R/sysdata.rda")

# Strain Isotype File
strain_isotype <- rio::import("https://docs.google.com/spreadsheets/d/1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI/pub?output=tsv") %>%
    dplyr::arrange(strain, isotype)

vcf_version <- 20160408

wb_build <- 245
vcf_path <- paste0("~/Dropbox/Andersenlab/Reagents/WormReagents/Variation/Andersen_VCF/WI.", vcf_version, ".impute.vcf.gz")
kinship <- generate_kinship(vcf_path)
snps <- generate_mapping(vcf_path)


# Save Datasets
devtools::use_data(kinship, snps, strain_isotype, mapping_snps, gene_ids, wb_build, vcf_version, internal = T, overwrite = T)
