#' Process QTL Intervals
#'
#' \code{variant_correlation} Returns all highly correlated variants in a QTL confidence interval
#'
#' This function losely identifies unique intervals for QTL in a data set and browses a whole-genome variant set for variants that are highly correlated with the phenotype that gave rise to the QTL. Spearman rank correlation is used. 
#' Heterozygotes are ignored in the analysis. Correlations are only calculated for variants that are present in at least 5% of the assayed strains. Additionally, variant information needs to have been
#' acquired for at least 80% of the phenotyped strains, this removes the possibility of discrepency between correlated variants and the correlation that led to the QTL.
#'
#' @param df is a dataframe that is output from the \code{process_mappings} function
#' @param quantile_cutoff is a quantile cutoff that determines what variants to keep, default is to keep all variants with correlation coefficients greater than the 90th quantile
#' @param variant_severity what variants to look at from snpeff output
#' @param gene_types what gene types to look at from snpeff output
#' @param kin is a strain by strain relatedness matrix you want to correct you trait with
#' @return Outputs a list. Each list contains two data frames, the first contains mapping information (e.g. log10p, confidence interval start and stop), phenotype information, and gene ids. 
#' The second element of the list contains more detailed gene information
#' @importFrom dplyr %>%
#' @export

variant_correlation <- function(df, 
                                quantile_cutoff = 0.9, 
                                variant_severity = c("MODERATE", "SEVERE"),
                                gene_types = "ALL",
                                kin = kinship,
                                condition_trait) {
  
  # source("~/Dropbox/Andersenlab/WormReagents/Variation/Andersen_VCF/read_vcf.R") # Get snpeff function
  
  
  # loosely identify unique peaks
  intervals <- df %>% na.omit() %>% 
    dplyr::distinct(CHROM, startPOS, endPOS, trait, .keep_all = TRUE) %>%
    dplyr::distinct(CHROM, startPOS, trait, .keep_all = TRUE) %>% 
    dplyr::distinct(CHROM, endPOS, trait, .keep_all = TRUE) %>% 
    dplyr::arrange(CHROM,  startPOS)
  
  if(condition_trait == T){
    intervals <- df %>% na.omit() %>% 
      tidyr::separate(trait, into = c("condition", "trait"), sep = "_")%>%
      dplyr::distinct(CHROM, startPOS, endPOS, condition, .keep_all = TRUE) %>%
      dplyr::distinct(CHROM, startPOS, condition, .keep_all = TRUE) %>% 
      dplyr::distinct(CHROM, endPOS, condition, .keep_all = TRUE) %>% 
      dplyr::arrange(CHROM,  startPOS)%>%
      tidyr::unite(trait, condition, trait, sep = "_")
  }
  
  strains <- as.character(na.omit(unique(df$strain)))
  intervalGENES <- list()
  
  for (i in 1:nrow(intervals)) {
    print(paste(100 * signif(i/nrow(intervals), 3), "%", 
                sep = ""))
    nstrains <- data.frame(df) %>% na.omit() %>% dplyr::filter(trait == 
                                                                 as.character(intervals[i, "trait"]))
    nstrains <- length(unique(nstrains$strain))
    chr <- as.character(intervals[i, ]$CHROM)
    left <- intervals[i, ]$startPOS
    right <- intervals[i, ]$endPOS
    region_of_interest <- paste0(chr, ":", left, "-", right)
    snpeff_output <- snpeff(region = region_of_interest, severity = variant_severity, elements = gene_types)
    pruned_snpeff_output <- snpeff_output %>% 
      dplyr::filter(strain %in% strains) %>% 
      dplyr::filter(!is.na(impact)) %>% 
      dplyr::distinct(CHROM, POS, strain, effect, gene_id, .keep_all = TRUE) %>% 
      dplyr::arrange(effect) %>% 
      dplyr::select(CHROM, POS, REF, ALT, GT, effect, nt_change, 
                    aa_change, gene_name, gene_id, transcript_biotype, 
                    feature_type, strain)  %>% 
      dplyr::group_by(CHROM, POS, effect) %>%
      dplyr::filter(!is.na(GT), GT != "HET") %>% 
      dplyr::mutate(num_allele = ifelse(GT == "REF", 0, ifelse(GT == "ALT", 1, NA))) %>% 
      dplyr::mutate(num_alt_allele = sum(num_allele,  na.rm = T), num_strains = n()) %>% 
      dplyr::filter(num_alt_allele/num_strains > 0.05) %>% 
      dplyr::filter(num_strains > nstrains * 0.8) %>% 
      dplyr::ungroup()
    if (nrow(pruned_snpeff_output) > 0) {
      interval_df <- df %>% 
        dplyr::filter(CHROM == chr, startPOS == left, endPOS == right) %>% 
        dplyr::group_by(trait, CHROM, startPOS, endPOS) %>% 
        dplyr::filter(log10p ==  max(log10p)) %>% 
        dplyr::distinct(trait, startPOS, endPOS, peakPOS, strain, .keep_all = TRUE) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(trait, startPOS, endPOS, peakPOS, 
                      strain, log10p, pheno_value = value) %>%
        na.omit()
      
      
      
      
      correct_it <- list()
      for(j in 1:length(unique(interval_df$trait))){
        
        temp_pheno <- dplyr::filter(interval_df, trait == unique(interval_df$trait)[j])%>%
          dplyr::select(trait, strain, value = pheno_value)
        
        correct_it[[j]] <- kinship_correction(temp_pheno) %>%
          dplyr::mutate(trait = unique(interval_df$trait)[j])
      }
      
      correct_df <- dplyr::bind_rows(correct_it)
      
      
      pheno_snpeff_df <- pruned_snpeff_output %>% 
        dplyr::left_join(., interval_df, by = "strain", copy = TRUE) %>% 
        dplyr::distinct(strain, trait, pheno_value, gene_id,  CHROM, POS, aa_change, .keep_all = TRUE) %>% 
        dplyr::group_by(trait, CHROM, POS, effect, feature_type) %>% 
        # na.omit()%>%
        dplyr::left_join(., correct_df, by=c("strain","trait")) %>%
        dplyr::filter(!is.na(trait))%>%
        dplyr::mutate(corrected_spearman_cor_p = cor.test(corrected_pheno, num_allele, method = "spearman", use = "pairwise.complete.obs",exact = F)$p.value,
                      spearman_cor_p = cor.test(pheno_value, num_allele, method = "spearman", use = "pairwise.complete.obs",exact = F)$p.value) %>% 
        dplyr::ungroup() %>% 
        dplyr::filter(corrected_spearman_cor_p < quantile(corrected_spearman_cor_p, 
                                                  probs = quantile_cutoff, na.rm = T)) %>% 
        dplyr::ungroup() %>% 
        dplyr::arrange(corrected_spearman_cor_p)
      
      intervalGENES[[i]] <- pheno_snpeff_df
    }
    else {
      intervalGENES[[i]] <- NA
    }
  }
  return(intervalGENES)
}

#' Combine Variant Correlation Data
#'
#' \code{process_correlations} Combines all variant correlation information with phenotype data from GWAS mappings
#'
#' Function to combine phenotype, mapping, and gene information from \code{variant_correlation} function output. 
#'
#' @param df is a list object that is output from the \code{variant_correlation} function
#' @param gene_information is a data.frame that contains gene information, the column gene_id (WBGene....) must be present to join to snpeff output. Default is a data frame that contains gene_id, concise_description, provisional_description, and gene_class_description
#' @return Outputs a data frame that contains phenotype data, mapping data, and gene information for highly correlated variants in a particular QTL confidence interval.
#' @export

process_correlations <- function(df, gene_information = gene_functions){
  
  cors <- list()
  for (i in 1:length(df)) {
    if (length(unique(df[[i]])) > 1) {
      cors[[i]] <- df[[i]] %>% dplyr::filter(!grepl("green", 
                                                    trait)) %>% dplyr::filter(trait != "")
    }
  }
  variant_pheno <- dplyr::bind_rows(cors) %>% dplyr::select(CHROM,POS, REF, ALT, nt_change, aa_change, gene_name, transcript_biotype, 
                                                            gene_id, effect, num_alt_allele, num_strains, strain, 
                                                            GT, trait, pheno_value, startPOS, endPOS, log10p, spearman_cor_p, corrected_spearman_cor_p) %>%
    dplyr::arrange(corrected_spearman_cor_p, desc(pheno_value)) %>% 
    dplyr::distinct(CHROM, POS, REF,  ALT, strain, gene_id, trait, .keep_all = TRUE) %>% 
    dplyr::arrange(corrected_spearman_cor_p,  CHROM, POS, desc(pheno_value)) %>% 
    dplyr::left_join(., gene_information, by = "gene_id")
  return(variant_pheno)
}

#' Interval Summary
#'
#' \code{interval_summary} Summarizes variation within an interval.
#' Can also be used to summarize regions by wormbase identifier or locus name.
#' 
#' @param query A query such as \code{pot-2}, \code{I:1-1000}, \code{WBGene00010195}.
#' @param filter_variants Filter out poor quality variants. Default is TRUE. Not compatible with imputed data.
#' @param impute Use imputed data. Default is FALSE.
#' @return data frame of summarized information.
#' @examples interval_summary("pot-1")
#' @export

interval_summary <- function(query, filter_variants = T, impute = F) {
  elegans_gff <- get_db()
  if (!grepl("(I|II|III|IV|V|X|MtDNA).*", query)) {
    interval <- (dplyr::collect(elegans_gff %>%
                                  dplyr::filter(locus == query | gene_id == query | sequence_name == query) %>%
                                  dplyr::select(chrom, start, end) %>% 
                                  dplyr::summarize(chrom = chrom, start = min(start), end = max(end))) %>%
                   dplyr::mutate(interval = paste0(chrom, ":", start, "-", end)))$interval
  } else {
    interval <- query
  }
  
  # Parse interval
  q_interval <- stringr::str_split(interval, "[\\.:-]+")[[1]]
  qchrom <- q_interval[[1]]
  qstart <- as.integer(q_interval[[2]])
  qend   <- as.integer(q_interval[[3]])
  
  region_elements <- dplyr::collect(elegans_gff %>%
                                      dplyr::filter(chrom == qchrom, start >= qstart, end <= qend) %>% 
                                      dplyr::group_by(biotype) %>%
                                      dplyr::select(biotype, locus, gene_id) %>%
                                      dplyr::distinct(.keep_all = TRUE)) %>%
    dplyr::mutate(locus = ifelse(is.na(locus), gene_id, locus)) %>%
    dplyr::select(biotype, locus) %>%
    dplyr::group_by(biotype)
  
  variants <- snpeff(interval, severity = "ALL", elements = "ALL", impute = impute) %>%
    dplyr::filter(GT != "REF") 
  
  if (filter_variants == TRUE & impute == FALSE) {
    variants <- dplyr::filter(variants, FILTER == "PASS", FT == "PASS")
  }
  
  variant_gene_effects <- variants %>% dplyr::select(gene_id, gene_name, effect, impact) %>% 
    dplyr::mutate(gene_name = ifelse(gene_name == "", gene_id, gene_name)) %>%
    dplyr::group_by(gene_id, effect) %>%
    dplyr::filter(!grepl("-", gene_id)) %>% # Remove intergene regions
    dplyr::mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
    dplyr::distinct(.keep_all = TRUE) 
  
  variant_gene_summary <- variant_gene_effects %>%
    dplyr::group_by(gene_name, effect, impact) %>%
    dplyr::summarize(n = n())
  
  # Calculate max severity variants
  mvariants <- variants %>% dplyr::group_by(CHROM, POS, effect) %>%
    dplyr::select(CHROM, POS, impact, effect) %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    tidyr::spread(effect, impact) %>%
    dplyr::mutate_each(dplyr::funs(sev), matches("_")) %>%
    dplyr::ungroup()
  
  mvariants$max_severity <- rseverity[apply(mvariants %>% dplyr::select(-CHROM, -POS) %>% t(), 2, max)]
  
  # Take all annotations and filter for the maximal impact at each position (so only maximal variant effects are counted)
  max_severity <- mvariants %>% dplyr::group_by(max_severity) %>% dplyr::summarize(n = n()) 
  
  # Calculate raw number of annotations
  variant_summary <- as.data.frame(t(mvariants %>% dplyr::select(-CHROM, -POS, -max_severity) %>% dplyr::summarize_each(dplyr::funs(sum(. > 0))))) %>%
    dplyr::add_rownames() %>%
    dplyr::rename(effect = rowname, n_variants = V1) %>%
    dplyr::left_join(variant_gene_summary, .)
  
  # Count total number of snps
  total_snps <- nrow(variants %>% dplyr::select(CHROM, POS) %>% dplyr::distinct(.keep_all = TRUE))
  
  # Genes with variants:
  total_genes <- nrow(region_elements)
  genes_w_variants <- length(unique(unlist((variant_gene_summary %>% dplyr::filter(impact != "MODIFIER"))$gene_name)))
  genes_w_MOD_HIGH <- length(unique(unlist((variant_gene_summary %>% dplyr::filter(impact %in% c("MODERATE","HIGH")))$gene_name)))
  genes_w_HIGH <-  length(unique(unlist((variant_gene_summary %>% dplyr::filter(impact %in% c("HIGH")))$gene_name)))
  
  results <-  list("query" = query,
                   "region" = interval,
                   "chrom" = qchrom,
                   "start" = qstart,
                   "end" = qend,
                   "snps" = total_snps,
                   "genes" = total_genes,
                   "genes_w_variants" = genes_w_variants,
                   "genes_w_MOD_HIGH" = genes_w_MOD_HIGH,
                   "genes_w_HIGH" = genes_w_HIGH,
                   "variant_gene_effects" = variant_gene_effects, 
                   "variant_summary" = variant_summary,
                   "max_severity_variants" = max_severity,
                   "region_genes" = region_elements
  )
  results
  
}

tf <- function(x) { ifelse(!is.na(x), T, F)}
severity <- list("MODIFIER" = 1, "LOW" = 2, "MODERATE" = 3, "HIGH" = 4)
rseverity <- names(severity)
names(rseverity) <- 1:4
sev <- function(col) { sapply(col, function(x) { 
  if (is.na(x)) { 
    0              
  } else {
    severity[[x]] 
  }
})
}


# # y is a data.frame that contains a trait, strain, and value column
# # function is adapted from the kinship.on.the.fly function from the cape package
kinship_correction <- function(y, kin = kinship) {
  
  K <- kin[colnames(kin) %in% y$strain, rownames(kin) %in% y$strain]
  y <- y %>% dplyr::filter(strain %in% colnames(K))
  
  model = regress::regress(as.vector(y$value)~1,~K, pos = c(TRUE, TRUE))	
  
  #This err.cov is the same as err.cov in Dan's code using estVC
  err.cov = regress::summary.regress(model)$sigma[1]*K+regress::summary.regress(model)$sigma[2]*diag(nrow(K))
  
  eW = eigen(err.cov, symmetric = TRUE)
  
  if(min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)){
  }else{
    eW$values[eW$values <= 0] = Inf
  } 
  
  err.cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)
  new.pheno <- err.cov %*% y$value
  
  new.pheno <- data.frame(strain = y$strain, corrected_pheno = new.pheno)
  
  return(new.pheno)
}


