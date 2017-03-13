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
strain_isotype <- data.table::fread("https://docs.google.com/spreadsheets/d/1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI/pub?output=tsv") %>%
    dplyr::arrange(strain, isotype)

# Create Strain --> Isotype Table
strain_isotype_mapping <- dplyr::bind_rows(strain_isotype %>% 
  dplyr::select(previous_names, isotype) %>%
  dplyr::filter(complete.cases(.), isotype != "") %>%
  dplyr::mutate(previous_names = stringr::str_split(previous_names, "\\|")) %>%
  tidyr::unnest(previous_names) %>%
  dplyr::rename(strain = previous_names), 
strain_isotype %>%
  dplyr::select(strain, isotype) %>%
  dplyr::filter(complete.cases(.), isotype != ""))  %>%
  dplyr::distinct() %>%
  dplyr::select(strain, isotype)

vcf_version <- 20170312

wb_build <- 245
vcf_path <- paste0("~/Dropbox/Andersenlab/Reagents/WormReagents/Variation/Andersen_VCF/WI.", vcf_version, ".impute.vcf.gz")
kinship <- generate_kinship(vcf_path)
snps <- generate_mapping(vcf_path)


# Save Datasets
devtools::use_data(kinship, snps, strain_isotype, mapping_snps, gene_ids, wb_build, vcf_version, strain_isotype_mapping, internal = F, overwrite = T)
