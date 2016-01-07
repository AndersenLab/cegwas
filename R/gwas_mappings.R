#' GWAS Mappings
#'
#' \code{gwas_mappings} uses the rrBLUP package function \code{GWAS} to perform association mapping
#' for C. elegans. Uses 5\% MAF SNPs from RADseq dataset from Andersen et al. 2012 and a 
#' kinship matrix generated from whole-genome sequence data. Can use user developed kinship matrix as well. 
#'
#' This is the detail section if you want to fill out in the future
#'
#' @param data two element list. element 1 : traits. element 2: trait values with strains in columns
#' with each row corresponding to trait in element 1
#' @param cores number of cores on computer that you want to allocate for mapping. Default value is 4
#' @param kin_matrix is a strainXstrain matrix. default kinship matrix is described above.
#' @param min.MAF Allele frequency at which to remove SNPs from mapping snpset (dependent on available strains).
#' @param mapping_snp_set is a logical statement that, if TRUE, uses a pruned SNP set that was identified from simulation studies.
#' @return Outputs a data frame with the following columns : marker, CHROM, POS, trait, log10p, where marker is CHROM_POS.
#' @importFrom foreach %dopar%
#' @export

gwas_mappings <- function(data, cores = parallel::detectCores(), kin_matrix = kinship, snpset = snps, min.MAF = 0.05, mapping_snp_set = TRUE){
  # phenotype prep
  x <- data.frame(trait = data[[1]], data[[2]]) %>%
    tidyr::gather(strain, value, -trait) %>%
    tidyr::spread(trait, value) # extract phenotypes from phenotype object
  
  # add marker column to snp set
  y <- data.frame(marker = paste(snpset$CHROM,snpset$POS,sep="_"),
                  snpset)

  # option to only keep snp set identified using simulations
  if(mapping_snp_set == TRUE){
    y <- y %>% 
      dplyr::filter(marker %in% data.frame(mapping_snps)$CHROM_POS)
  }
  
  kin <- as.matrix(kin_matrix)
  
  strains <- data.frame(strain = x[,1])
  
  x <- data.frame(x[,2:ncol(x)])
  
  colnames(x) <- data[[1]]
  
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  # run mappings
  system.time(
    maps <- foreach::foreach(i = iterators::iter(x, by = 'col')) %dopar% {
      rrBLUP::GWAS(pheno = data.frame(strains, i),
                   geno = y,
                   K = kin,
                   min.MAF = min.MAF,
                   n.core = 1,
                   P3D = FALSE,
                   plot = FALSE)
    })
  
  parallel::stopCluster(cl)   
  
  for(i in 1:ncol(x)){
    colnames(maps[[i]]) <- c("marker", "CHROM",  "POS",   "log10p")
    maps[[i]]$trait <-  colnames(x)[i]
  }
  
  mappings <- dplyr::bind_rows(maps)
  
  return(mappings)
}


#' cegwas_map
#'
#' \code{cegwas_map} is a convenience function takes trait data for a set of strains and performs
#' a mapping - returning a processed mapping data frame. \code{map} wraps \code{\link{process_pheno}}, \code{\link{gwas_mappings}}, 
#' and \code{\link{process_mappings}} into a single function. 
#'
#' @param data two element list. element 1 : traits. element 2: trait values with strains in columns
#' with each row corresponding to trait in element 1
#' @param cores number of cores on computer that you want to allocate for mapping. Default value is 4
#' @param BF defines a custom bonferroni correction.
#' @param remove_strains Remove strains with no known isotype. Default is FALSE.
#' @param duplicate_method Method for dealing with the presence of multiple strains falling into the same isotype. Either \code{"average"} to average phenotypes or \code{"first"} to take the first observation.
#' @param kin_matrix is a strainXstrain matrix. default kinship matrix is described above.
#' @param snps is a set of mapping snps.
#' @param mapping_snp_set Use simulation based snps when TRUE. Use 5\% cut when FALSE.
#' @return Outputs a two element list that contains two dataframes. 
#' The first data frame is a processed mappings dataframe that contains the same columns
#' as the output of \code{\link{gwas_mappings}} with two additional columns. One that contains
#' the bonferroni corrected p-value (BF) and another that contains an identifier 1,0 if 
#' the indicated SNP has a higher -log10 value than the bonferroni cut off or not, respectively
#' The second data frame contains the variance explained data as well as all of the information from the first element.
#' @export

cegwas_map <- function(trait_data, 
                       cores = parallel::detectCores(),
                       remove_strains = TRUE, 
                       kin_matrix = kinship,
                       snpset = snps,
                       duplicate_method = "first",
                       BF = NA,
                       mapping_snp_set = TRUE) {
  processed_phenotypes <- process_pheno(trait_data, remove_strains = remove_strains, duplicate_method = "first")
  mapping_df <- gwas_mappings(processed_phenotypes, kin_matrix = kin_matrix, snpset = snpset, cores = cores, mapping_snp_set = mapping_snp_set)
  processed_mapping_df <- process_mappings(mapping_df, phenotype_df = processed_phenotypes, CI_size = 50, snp_grouping = 200, BF = BF)
}




