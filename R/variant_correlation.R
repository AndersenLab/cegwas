#' Process QTL Intervals
#'
#' \code{variant_correlation} Returns all highly correlated variants in a QTL confidence interval
#'
#' This function losely identifies unique intervals for QTL in a data set and browses a whole-genome variant set for variants that are highly correlated with the phenotype that gave rise to the QTL. Spearman rank correlation is used. 
#' Heterozygotes are ignored in the analysis. Correlations are only calculated for variants that are present in at least 5% of the assayed strains. Additionally, variant information needs to have been
#' acquired for at least 80% of the phenotyped strains, this removes the possibility of discrepency between correlated variants and the correlation that led to the QTL.
#'
#' @param df is a dataframe that is output from the \code{process_mappings} function
#' @param quantile_cutoff_high is a quantile cutoff that determines what variants to keep, default is to keep all variants with correlation coefficients greater than the 90th quantile
#' @param quantile_cutoff_low is a quantile cutoff that determines what variants to keep, default is to keep all variants with correlation coefficients less than the 10th quantile
#' @param variant_severity what variants to look at from snpeff output
#' @param gene_types what gene types to look at from snpeff output
#' @param kin is a strain by strain relatedness matrix you want to correct you trait with
#' @return Outputs a list. Each list contains two data frames, the first contains mapping information (e.g. log10p, confidence interval start and stop), phenotype information, and gene ids. 
#' The second element of the list contains more detailed gene information
#' @importFrom dplyr %>%
#' @export

variant_correlation <- function(df, 
                                quantile_cutoff_high = .9, 
                                quantile_cutoff_low = .1,
                                variant_severity = c("MODERATE", "SEVERE"),
                                gene_types = "ALL",
                                kin = cegwas::kinship){
  
  # source("~/Dropbox/Andersenlab/WormReagents/Variation/Andersen_VCF/read_vcf.R") # Get snpeff function
  
  
  # loosely identify unique peaks
  intervals <- df %>% na.omit() %>% 
    dplyr::distinct(CHROM, startPOS, endPOS) %>%
    dplyr::distinct(CHROM, startPOS) %>% 
    dplyr::distinct(CHROM, endPOS) %>% 
    dplyr::arrange(CHROM,  startPOS)
  
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
      dplyr::distinct(CHROM, POS, strain, effect, gene_id) %>% 
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
        dplyr::distinct(trait, startPOS, endPOS, peakPOS, strain) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(trait, startPOS, endPOS, peakPOS, 
                      strain, log10p, pheno_value = value) %>%
        na.omit()
      
      correct_it <- list()
      for(j in 1:length(unique(interval_df$trait))){
        
        temp_pheno <- filter(interval_df, trait == unique(interval_df$trait)[j])
        
        k <- kin[row.names(kin)%in%temp_pheno$strain,colnames(kin)%in%temp_pheno$strain]
        
        correct_it[[j]] <- data.frame(t(temp_pheno$pheno_value) %*% k) %>%
          tidyr::gather(strain, corrected_pheno)  %>%
          dplyr::mutate(trait = unique(interval_df$trait)[j])
      }
      
      correct_df <- dplyr::rbind_all(correct_it)
      
      
      pheno_snpeff_df <- pruned_snpeff_output %>% 
        dplyr::left_join(., interval_df, by = "strain", copy = TRUE) %>% 
        dplyr::distinct(strain, trait, pheno_value, gene_id,  CHROM, POS, aa_change) %>% 
        dplyr::group_by(trait, CHROM, POS, effect, feature_type) %>% 
        na.omit()%>%
        dplyr::left_join(., correct_df, by=c("strain","trait")) %>%
        dplyr::mutate(corrected_spearman_cor = cor(corrected_pheno, num_allele, method = "spearman", use = "pairwise.complete.obs"),
                      spearman_cor = cor( pheno_value, num_allele, method = "spearman", use = "pairwise.complete.obs")) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(abs_spearman_cor = abs(corrected_spearman_cor)) %>% 
        dplyr::filter(abs_spearman_cor > quantile(abs_spearman_cor, 
                                                  probs = quantile_cutoff_high, na.rm = T)) %>% 
        dplyr::ungroup() %>% 
        dplyr::arrange(desc(abs_spearman_cor))
      
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
  variant_pheno <- dplyr::rbind_all(cors) %>% dplyr::select(CHROM,POS, REF, ALT, nt_change, aa_change, gene_name, transcript_biotype, 
                                                            gene_id, effect, num_alt_allele, num_strains, strain, 
                                                            GT, trait, pheno_value, startPOS, endPOS, log10p, spearman_cor, corrected_spearman_cor,
                                                            abs_spearman_cor) %>%
    dplyr::arrange(desc(abs_spearman_cor), desc(pheno_value)) %>% 
    dplyr::distinct(CHROM, POS, REF,  ALT, strain, gene_id, trait) %>% dplyr::arrange(desc(abs_spearman_cor),  CHROM, POS, desc(pheno_value)) %>% 
    dplyr::left_join(., gene_information, by = "gene_id")
  return(variant_pheno)
}


#' Browse Variant Info
#'
#' \code{snpeff} enables you to query variants called and annotated by the \href{http://www.andersenlab.org}{Andersen Lab}. 
#' 
#' @param ... Gene names, regions, or wormbase identifiers to query.
#' @param severity A vector with variants of given severities (LOW, MODERATE, HIGH, MODIFIER). Default takes moderate and high. Use "ALL" to return all variants.
#' @param elements A vector containing gene structural elements (CDS, five_prime_UTR, exon, intron, three_prime_UTR). Use "ALL" to return all variants.
#' @param long Return dataset in long or wide format. Default is to return in long format.
#' @param remote Use remote data. Checks for local data if possible. False by default.
#' @param impute Use imputed data. Default is FALSE.
#' @return Outputs a data frame that contains phenotype data, mapping data, and gene information for highly correlated variants in a particular QTL confidence interval.
#' @examples snpeff("pot-2","II:1-10000","WBGene00010785")
#' @export

snpeff <- function(...,
                   severity = c("HIGH","MODERATE"),
                   elements = c("exon"),
                   long = TRUE,
                   remote = FALSE,
                   impute = FALSE) {
  
  regions <- unlist(list(...))
  
  # Allow user to specify 'ALL'
  if ("ALL" %in% severity) {
    severity <-  c("LOW", "MODERATE", "HIGH", 'MODIFIER')
  }
  if ("ALL" %in% elements) {
    elements <- c("CDS", "five_prime_UTR", "exon", "intron", "three_prime_UTR")
  }
  
  # Ensure that bcftools is available:
  bcftools_version <- as.double(stringr::str_extract(readLines(pipe("bcftools --version"))[1], "[0-9]+\\.[0-9]+"))
  if(is.na(bcftools_version) | bcftools_version < 1.2) {
    stop("bcftools 1.2+ required for this function")
  }
  
  results <- suppressWarnings(lapply(regions, function(query) {
  # Save region as query

  # Fix region specifications
  query <- gsub("\\.\\.", "-", query)
  query <- gsub(",", "", query)

  # Resolve region names
  if (!grepl("(I|II|III|IV|V|X|MtDNA).*", query)) {
    elegans_gff <- get_db()
    # Pull out regions by element type.
    region <- paste((dplyr::bind_rows(lapply(elements, function(e) {
      dplyr::collect(dplyr::filter(elegans_gff, locus == query | gene_id == query | sequence_name == query, type_of == e) %>%
                       dplyr::select(chrom, start, end, gene_id, biotype, type_of, locus, sequence_name) %>%
                       dplyr::distinct())
      })) %>%
      dplyr::summarize(chrom = chrom[1], start = min(start), end = max(end)) %>%
      dplyr::mutate(region_format = paste0(chrom, ":", start, "-", end)) %>%
        dplyr::select(region_format) %>%
        dplyr::distinct())$region_format, collapse = ",")
      if (stringr::str_length(regions[[1]]) == 0) {
        message(paste0(query, " not found."))
        region <- NA
      }
  } else {
    region <- query
  }

  vcf_path <- get_vcf(remote = remote, impute = impute)
  
  sample_names <- readr::read_lines(suppressWarnings(pipe(paste("bcftools","query","-l",vcf_path))))
  base_header <- c("CHROM", "POS", "REF","ALT","FILTER")
  ANN_header = c("allele", "effect", "impact",
                "gene_name", "gene_id", "feature_type", 
                "feature_id", "transcript_biotype","exon_intron_rank",
                "nt_change", "aa_change", "cDNA_position/cDNA_len", 
                "protein_position", "distance_to_feature", "error", "extra")
  # If using long format provide additional information.
  format <- "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%ANN[\\t%TGT]\\n'"
  if (long == T) {
    format <- "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%ANN[\\t%TGT!%FT!%DP!%DP4!%SP!%HP]\\n'"
  }
  command <- paste("bcftools","query","--regions", region, "-f", format ,vcf_path)
  if (!is.na(region)) {
    message(paste0("Query: ", query, "; region - ", region, "; "))
    result <- try(dplyr::tbl_df(data.table::fread(command, col.names = c(base_header, "ANN", sample_names ), sep = "\t")), silent = TRUE)
    if(!grepl("^Error.*", result[[1]][1])) {
    tsv <- result %>%
           dplyr::mutate(REF = ifelse(REF==TRUE, "T", REF), # T nucleotides are converted to 'true'
                  ALT = ifelse(ALT==TRUE, "T", ALT))
    } else {
      tsv <- as.data.frame(NULL)
    }
  # If no results are returned, stop.
    if (typeof(tsv) == "character" | nrow(tsv) == 0) {
      warning("No Variants")
      NA
    } else {
      tsv <-  dplyr::mutate(tsv, ANN=strsplit(ANN,",")) %>%
              tidyr::unnest(ANN) %>%
              tidyr::separate(ANN, into = ANN_header, sep = "\\|") %>%
              dplyr::select(one_of(c(base_header, ANN_header)), everything(), -extra) %>%
              dplyr::mutate(gene_name = as.character(gene_ids[gene_name])) %>%
              dplyr::mutate(query = query, region = region) %>%
              dplyr::select(CHROM, POS, query, region, everything())
      
      tsv <-  dplyr::filter(tsv, impact %in% severity) 
      if (nrow(tsv) == 0) {
        message(paste("No Results for", region, "after filtering"))
      }
      if (long == FALSE) {
        tsv
      } else {
        tsv <- tidyr::gather_(tsv, "strain", "GT", names(tsv)[23:length(tsv)])  %>%
          tidyr::separate(GT, into=c("a1","a2", "FT", "DP", "DP4", "SP", "HP"), sep="/|\\||\\!", remove=T) %>%
          dplyr::mutate(a1=ifelse(a1 == ".", NA, a1)) %>%
          dplyr::mutate(a2=ifelse(a2 == ".", NA, a2)) %>%
          dplyr::mutate(GT = NA) %>%
          dplyr::mutate(GT = ifelse(a1 == REF & a2 == REF & !is.na(a1), "REF",GT)) %>%
          dplyr::mutate(GT = ifelse(a1 != a2 & !is.na(a1), "HET",GT)) %>%
          dplyr::mutate(GT = ifelse(a1 == a2 & a1 != REF & !is.na(a1), "ALT",GT)) %>%
          dplyr::select(CHROM, POS, strain, REF, ALT, a1, a2, GT, FT, FILTER, DP, DP4, SP, HP, everything()) %>%
          dplyr::arrange(CHROM, POS) 
      }
        tsv
    }
  }
  }
  ))
  results <- do.call(rbind, results)
  if (!"CHROM" %in% names(results)){
    stop("No Results")
  }
  results <- results %>% dplyr::filter(!is.na(CHROM))
  results
}

#' Generate directory path to VCF file
#'
#' \code{get_vcf} Generate directory path to VCF file
#' 
#' @param remote logical, to use remote VCF stored on google or local in Andersen lab dropbox
#' @param impute logical, to use imputed VCF file or non imputed
#' @return character value corresponding to VCF location
#' @export

get_vcf <- function(remote = F, impute = T) {
  
  # Set vcf path; determine whether local or remote
  if (impute == F) {
    vcf_name = paste0("WI.", vcf_version, ".snpeff.vcf.gz")
  } else {
    vcf_name = paste0("WI.", vcf_version, ".impute.snpeff.vcf.gz")
  }
  
  
  vcf_path <- paste0("~/Dropbox/Andersenlab/Reagents/WormReagents/Variation/Andersen_VCF/", vcf_name)
  # Use remote if not available.
  local_or_remote <- "locally"
  if (!file.exists(vcf_path) | remote == T) {
    vcf_path <- paste0("http://storage.googleapis.com/andersen_dist/vcf/all/", vcf_version, "/", vcf_name)
    message("Using remote vcf")
  }
  if (local_or_remote == "locally") {
    system(paste0("touch ", vcf_path,".csi"))
    message("Using local vcf")
  }
  vcf_path
}


#' Get Database
#'
#' Return database; If not available, download.
#'
#' @export


get_db <- function() {
  file_path <- paste0("~/.WS", wb_build, ".elegans_gff.db")
  if (file.info(file_path)$size < 128 | is.na(file.info(file_path)$size == 0)) {
    message(paste0("Downloading Gene Database to ", file_path))
    url <- paste0("http://storage.googleapis.com/cegwas/WS", wb_build , ".celegans_gff.db")
    download.file(url, file_path)
  }
  dplyr::tbl(dplyr::src_sqlite(paste0("~/.WS", wb_build, ".elegans_gff.db")), "feature_set")
}



#' Fetch Variant Type
#'
#' \code{fetch_id_type} Fetches wormbase identifiers of a certain class. These can be piped into \code{snpeff}
#' 
#' @param id_type Type of genetic element. Omit for list of options.
#' @return Vector of identifiers.
#' @examples fetch_id_type("lincRNA")
#' @export

fetch_id_type <- function(id_type = NA) {
    elegans_gff <- get_db()
    valid_id_types <- dplyr::collect(elegans_gff %>%
                                       dplyr::select(biotype) %>%
                                       dplyr::distinct())$biotype
    if (!(id_type %in% valid_id_types)) {
        message("Available ID Types:")
        cat(valid_id_types, sep = "\n")
    } else {
        (dplyr::collect(elegans_gff %>%
        dplyr::filter(biotype == id_type, type_of == "exon") %>%
        dplyr::select(gene_id)))$gene_id
    }
  
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
    dplyr::distinct()) %>%
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
    dplyr::distinct() 

  variant_gene_summary <- variant_gene_effects %>%
      dplyr::group_by(gene_name, effect, impact) %>%
      dplyr::summarize(n = n())

  # Calculate max severity variants
  mvariants <- variants %>% dplyr::group_by(CHROM, POS, effect) %>%
              dplyr::select(CHROM, POS, impact, effect) %>%
              dplyr::distinct() %>%
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
  total_snps <- nrow(variants %>% dplyr::select(CHROM, POS) %>% dplyr::distinct())

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

