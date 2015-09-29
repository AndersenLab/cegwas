
# gene Ids
library("readr")

# Download
url <- "ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/c_elegans.PRJNA13758.current.geneIDs.txt.gz"
tmp <- tempfile()
download.file(url,tmp)
gene_info <- readr::read_csv(gzfile(tmp), col_names = c("taxon","WBID","name","protein","live")) %>%
  select(-taxon)

gene_ids <- as.list(gene_info$name)
names(gene_ids) <- gene_info$WBID

save(gene_ids, file = "data/gene_ids.Rda")

# Isotypes

strain_isotype <- read_tsv("data-raw/mapping_strain_isotype.tsv")
save(strain_isotype, file = "data/strain_isotype.Rda")
