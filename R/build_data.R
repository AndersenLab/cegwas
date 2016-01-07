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
strain_isotype <- readr::read_tsv("inst/data-raw/mapping_strain_isotype.tsv") %>%
  dplyr::arrange(strain, isotype)

load("data/gene_ids.rda")
load("data/kinship.rda")
load("data/snps.rda")

# cegwas db; if db file not exist or outdated then download.
wb_build <- 245
file_path <- paste0("~/.WS", wb_build, ".elegans_gff.db")
if (file.info(file_path)$size ==0 | is.na(file.info(file_path)$size == 0)) {
  message(paste0("Downloading Gene Database to ", file_path))
  url <- paste0("http://storage.googleapis.com/cegwas/WS", wb_build, ".celegans_gff.db")
  download.file(url, file_path)
}


# Save Datasets
devtools::use_data(kinship, snps, strain_isotype, gene_ids, wb_build, internal = FALSE, overwrite = T)

