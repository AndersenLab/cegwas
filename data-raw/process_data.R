
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

# Generate kinship and mapping snp sets
snps <- generate_mapping("~/Dropbox/Andersenlab/Reagents/WormReagents/Variation/Andersen_VCF/20150731_WI_PASS.impute.gt.vcf.gz")
kinship <- generate_kinship("~/Dropbox/Andersenlab/Reagents/WormReagents/Variation/Andersen_VCF/20150731_WI_PASS.impute.gt.vcf.gz")

# Isotypes
strain_warnings <- read_tsv("data-raw/mapping_strain_isotype.tsv")
strain_isotype <- strain_warnings$isotype
names(strain_isotype) <- strain_warnings$strain
strain_warnings <- filter(strain_warnings, warning_msg != "") %>%
                   select(strain, warning_msg)
sw <- strain_warnings$warning_msg
names(sw) <- strain_warnings$strain
strain_warnings <- sw
devtools::use_data(kinship, snps, gene_ids, strain_isotype, strain_warnings, internal = TRUE, overwrite = T)

