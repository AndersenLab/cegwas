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

# Strain Isotype File
strain_isotype <- readr::read_tsv("http://storage.googleapis.com/andersen_lab_strains/processed/strain_isotype.tsv") %>%
  dplyr::arrange(strain, isotype)

load("data/gene_ids.rda")
load("data/kinship.rda")
load("data/snps.rda")

vcf_date <- "20160106"
vcf_path <- paste0("~/Dropbox/AndersenLab/Reagents/WormReagents/Variation/Andersen_VCF/WI.", vcf_date, ".vcf.gz")
kinship <- generate_kinship(vcf_path)
snps <- generate_mapping(vcf_path)

# Save Datasets
devtools::use_data(kinship, snps, strain_isotype, gene_ids, wb_build, internal = FALSE, overwrite = T)

# generate  isotype_location.rda file 
# isotype_location <- data.table::fread("~/AndersenLab/Github_Repos/Andersen-Lab-Strains/processed/strain_info_join.tsv",
#                                       header = T) %>%
#   filter(!is.na(longitude),sequenced==TRUE) %>%
#   distinct(isotype)
