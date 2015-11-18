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
#' @return Outputs a list. Each list contains two data frames, the first contains mapping information (e.g. log10p, confidence interval start and stop), phenotype information, and gene ids. 
#' The second element of the list contains more detailed gene information
#' @importFrom dplyr %>%
#' @export

variant_correlation <- function(df, 
                                quantile_cutoff_high = .9, 
                                quantile_cutoff_low = .1){
  
  # source("~/Dropbox/Andersenlab/WormReagents/Variation/Andersen_VCF/read_vcf.R") # Get snpeff function
  
  
  # loosely identify unique peaks
    intervals <- df %>%
      na.omit() %>%
      dplyr::distinct(CHROM, startPOS, endPOS) %>%
      dplyr::distinct(CHROM, startPOS ) %>%
      dplyr::distinct(CHROM, endPOS ) %>%
      dplyr::arrange(CHROM, startPOS) 
  
  
  # unique strains to filter snpeff output for GWAS data - doesnt matter for genomic traits
  strains <- as.character(na.omit(unique(df$strain)))
  # set up the database to search for gene annotations using the biomart package
  ensembl = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="celegans_gene_ensembl")
  
  # initialize a list to store gene annotations for genes most highly correlated with phenotype
  intervalGENES <- list()
  
  # loop through all unique intervals
  for( i in 1:nrow(intervals)){
    
    print(paste(100*signif(i/nrow(intervals),3), "%",sep=""))
    
    
    nstrains <- data.frame(df) %>%
      na.omit() %>%
      dplyr::filter(trait == as.character(intervals[i,"trait"]))
    
    nstrains <- length(unique(nstrains$strain))
    
    # define chromosome and left and right bound for intervals
    chr <- as.character(intervals[i,]$CHROM)
    left <- intervals[i,]$startPOS
    right <- intervals[i,]$endPOS
    
    # define region of interest for Dan's snpeff function input
    region_of_interest <- paste0(chr,":",left,"-",right)
    
    # run variant effect prediction function
    snpeff_output <- snpeff(region = region_of_interest) 
    
    # prune snpeff outputs
    pruned_snpeff_output <- snpeff_output %>%
      dplyr::filter(strain %in% strains) %>% # only keep strains used in mappings
      dplyr::filter(!is.na(impact)) %>% # remove rows with NA in impact (looked to be mostly splice variants)
      dplyr::distinct(CHROM, POS, strain, effect, gene_id) %>% # remove duplicates
      dplyr::arrange(effect) %>% 
      # pull out columns of interest
      dplyr::select(CHROM, POS, REF, ALT, GT, effect, nt_change, aa_change, gene_name, gene_id, feature_type, strain) %>%
      dplyr::group_by(CHROM, POS, effect) %>% # group for individual genes and effects
      # make numeric allele column, this makes HETs and NAs in the GT column NAs - these are excluded from the correlation analysis
      # need to elimante hets and NAs from GT
      dplyr::filter(!is.na(GT), GT != "HET")%>%
      # make numeric
      dplyr::mutate(num_allele = ifelse(GT == "REF", 0, 
                                        ifelse(GT == "ALT", 1, NA)))%>%
      # determine if any alleles are present in less that 5% of the population
      dplyr::mutate(num_alt_allele = sum(num_allele, na.rm=T),
                    num_strains = n())%>%
      # if they are, eliminate
      dplyr::filter(num_alt_allele / num_strains > .05) %>%
      dplyr::filter(num_strains > nstrains*.8) %>%
      dplyr::ungroup()
    
    if( nrow(pruned_snpeff_output) > 0 ){
      # pull unique interval from processed mapping DF to recover, phenotypes, strains, log10p, phenotype value
      # this is useful to pull out all intervals with the same confidence interval that were pruned above.
      interval_df <- df %>%
        dplyr::filter(CHROM == chr, startPOS == left, endPOS == right)%>% # filter for confidence interval of interest
        dplyr::group_by(trait, CHROM, startPOS,endPOS) %>% # group by unique phenotype and interval
        dplyr::filter(log10p == max(log10p)) %>% # pull out most significant snp to minimize redundancy
        dplyr::distinct(trait, startPOS, endPOS, peakPOS, strain) %>% 
        dplyr::ungroup() %>%
        dplyr::select(trait, startPOS, endPOS, peakPOS, strain, log10p, pheno_value = value) 

      
      # calculate the correlation between interval variants and the phenotype 
      # pull out only the most correlated genes 
      
      pheno_snpeff_df <- pruned_snpeff_output %>%
        dplyr::left_join(., interval_df, by = "strain", copy = TRUE) %>% # join snpeff variant df to phenotype df for a particular interval
        dplyr::distinct(strain, trait, pheno_value, gene_id,CHROM,POS, aa_change) %>% # remove redundancy
        dplyr::group_by(trait, CHROM, POS, effect, feature_type) %>% # group_by unique variant and phenotype
        dplyr::mutate(spearman_cor = cor(pheno_value, num_allele, method = "spearman", use = "pairwise.complete.obs"))%>% # calculate correlation
        dplyr::ungroup()%>% # ungroup to calculate quantiles of correlations
        # dplyr::mutate(q90 = quantile(spearman_cor, probs = .9, na.rm = T) )%>%
        # we want to keep high positively correlated and high negatively correlated variants
        dplyr::mutate(abs_spearman_cor = abs(spearman_cor))%>%
        dplyr::filter(abs_spearman_cor > quantile(abs_spearman_cor, probs = quantile_cutoff_high, na.rm = T) )%>%
        dplyr::ungroup()%>%
        # organize DF by correlation
        dplyr::arrange(desc(abs_spearman_cor))
      
      
      # get gene annotations usining biomart package
      # attributes are the columns you want to return
      # filters are the columns you want to filter by, in this case we want to filter by wormbase_gene - e.g. WBGene00012953
      # values are the values you want to be present in the filter column, i.e the genes you want information from
      # this is pulled from the highly correlated variant DF above
      # mart is defined above as the annoted c.elegans genome
      gene_annotations <- getBM(attributes=c('entrezgene','go_id',"external_gene_name",
                                             "external_transcript_name","gene_biotype",
                                             "transcript_biotype","description", "family_description",
                                             "name_1006","wormbase_gene"), 
                                filters = "wormbase_gene",
                                values = unique(pheno_snpeff_df$gene_id),
                                mart=ensembl) %>%
        dplyr::distinct(entrezgene, go_id) %>% # pu;; distinct genes
        dplyr::rename(gene_id = wormbase_gene) # change column for joining
      
      # attach the correlation coefficient to gene annotation data frame to minimize looking at multiple data frames
      gene_cors <- pheno_snpeff_df %>%
        dplyr::select(gene_id, spearman_cor)%>%
        dplyr::distinct(gene_id, spearman_cor) %>%
        dplyr::left_join(gene_annotations, ., by = "gene_id") %>%
        dplyr::arrange(desc(spearman_cor))
      
      # append phenotype-snpeff-correlation DF and gene annotation DF to list for every unique interval
      intervalGENES[[i]] <- list(pheno_snpeff_df, gene_cors)
    } 
    else
    {
      intervalGENES[[i]] <- list(NA, NA)
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
#' @return Outputs a data frame that contains phenotype data, mapping data, and gene information for highly correlated variants in a particular QTL confidence interval.
#' @export

process_correlations <- function(df){
  
  # initialize list
  cors <-list()
  genes <- list()
  
  # separate list objects and remove green traits
  for(i in 1:length(df)){
    
    if(length(unique(df[[i]])) > 1){
      
      cors[[i]] <- df[[i]][[1]] %>%
        dplyr::filter(!grepl("green",trait)) %>%
        dplyr::filter(trait != "")
      
      genes[[i]] <- df[[i]][[2]] 
      
    }
  }
  
  # bind phenotype data
  variant_pheno <- dplyr::rbind_all(cors)%>%
    dplyr::select(CHROM, POS, REF, ALT, aa_change, gene_name, gene_id, num_alt_allele, num_strains, strain, GT, trait, pheno_value, startPOS, endPOS, log10p, spearman_cor, abs_spearman_cor) %>%
    dplyr::arrange(desc(abs_spearman_cor),desc(pheno_value))
  
  
  # bind gene data and join to phenotype data
  max_cor <- dplyr::rbind_all(genes) %>%
    dplyr::select(gene_id, external_gene_name, external_transcript_name, description, family_description, name_1006)%>%
    dplyr::left_join(variant_pheno, ., by = "gene_id") %>%
    dplyr::distinct(CHROM, POS, REF, ALT, strain, gene_id, trait, external_transcript_name) %>%
    dplyr::arrange(desc(abs_spearman_cor), CHROM, POS, desc(pheno_value))
 
  return(max_cor) 
}

#' Browse Variant Info
#'
#' \code{snpeff} enables you to query variants called and annotated by the \href{http://www.andersenlab.org}{Andersen Lab}. 
#' 
#' @param regions query A gene name, region, or wormbase identifier to query.
#' @param severity A vector with variants of given severities (LOW, MODERATE, HIGH, MODIFIER). Default takes moderate and high.
#' @param long Return dataset in long or wide format. Default is to return in long format.
#' @param remote Use remote data. Checks for local data if possible. False by default.
#' @return Outputs a data frame that contains phenotype data, mapping data, and gene information for highly correlated variants in a particular QTL confidence interval.
#' @examples snpeff(c("pot-2","II:1-10000","WBGene00010785"))
#' @export

snpeff <- function(regions,
                   severity = c("HIGH","MODERATE"),
                   long = TRUE,
                   remote = FALSE) {
  
  results <- lapply(regions, function(region) {
  # Fix region to allow wb type spec.
  
  # Set vcf path; determine whether local or remote
  vcf_name = "WI.20151118.snpeff.vcf.gz"
  
  # Resolve region names
  if (!grepl("(I|II|III|IV|V|X|MtDNA).*", region)) {
    # Resolve WB name if necessary
    gene_id <- names(which(gene_ids == region))
    if (length(gene_id) == 1) {
      wb_id <- gene_id
      region <- wb_id
    } else {
      wb_id <- region
    }
    
    if (gene_ids[wb_id] != "") {
      gene_message <- paste0("(", gene_ids[wb_id][[1]], ")")
    } else {
      gene_message <- ""
    }
    message(paste("Looking up", region, gene_message, "position"))
    wb_url <- paste0("http://api.wormbase.org/rest/field/gene/",wb_id, "/location/")
    wb_ret <- httr::GET(wb_url, httr::add_headers("Content-Type"="application/json"))
    region <- httr::content(wb_ret)$location$genomic_position$data[[1]]$label
  } 
  
  # Fix region specifications
  region <- gsub("\\.\\.", "-", region)
  region <- gsub(",", "", region)
  
  vcf_path <- paste0("~/Dropbox/Andersenlab/Reagents/WormReagents/Variation/Andersen_VCF/", vcf_name)
  # Use remote if not available.
  local_or_remote <- "locally"
  if (!file.exists(vcf_path) | remote == T) {
    vcf_path <- paste0("http://storage.googleapis.com/andersen/", vcf_name)
    local_or_remote <- "remotely"
  }
  
  sample_names <- readr::read_lines(pipe(paste("bcftools","query","-l",vcf_path)))
  base_header <- c("CHROM", "POS", "REF","ALT","FILTER")
  ANN_header = c("allele", "effect", "impact",
                "gene_name", "gene_id", "feature_type", 
                "feature_id", "transcript_biotype","exon_intron_rank",
                "nt_change", "aa_change", "cDNA_position/cDNA_len", 
                "protein_position", "distance_to_feature", "error", "extra")
  command <- paste("bcftools","query","--regions", region, "-f", "'%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%ANN[\t%TGT]\n'",vcf_path)
  
  message(paste("Querying", region, local_or_remote))
  tsv <- try(readr::read_tsv(pipe(command), col_names = c(base_header, "ANN", sample_names)), silent = T) %>%
         dplyr::mutate(REF = ifelse(REF==TRUE, "T", REF), # T nucleotides are converted to 'true'
                ALT = ifelse(ALT==TRUE, "T", ALT))
  # If no results are returned, stop.
  if (typeof(tsv) == "character") {
    stop("No Variants")
  } 
  tsv <-  dplyr::mutate(tsv, ANN=strsplit(ANN,",")) %>%
          tidyr::unnest(ANN) %>%
          tidyr::separate(ANN, into = ANN_header, sep = "\\|") %>%
          dplyr::select(one_of(c(base_header, ANN_header)), everything(), -extra) %>%
          dplyr::mutate(gene_name = as.character(gene_ids[gene_name]))
  
  tsv <-  dplyr::filter(tsv, impact %in% severity) 
  if (nrow(tsv) == 0) {
    message(paste("No Results for", region, "after filtering"))
  }
  if (long == FALSE) {
    tsv
  } else {
    tsv <- tidyr::gather_(tsv, "strain", "GT", names(tsv)[21:length(tsv)])  %>%
      tidyr::separate(GT, into=c("a1","a2"), sep="/|\\|", remove=T) %>%
      dplyr::mutate(a1=ifelse(a1 == ".", NA, a1)) %>%
      dplyr::mutate(a2=ifelse(a2 == ".", NA, a2)) %>%
      dplyr::mutate(GT = NA) %>%
      dplyr::mutate(GT = ifelse(a1 == REF & a2 == REF & !is.na(a1), "REF",GT)) %>%
      dplyr::mutate(GT = ifelse(a1 != a2 & !is.na(a1), "HET",GT)) %>%
      dplyr::mutate(GT = ifelse(a1 == a2 & a1 != REF & !is.na(a1), "ALT",GT)) %>%
      dplyr::select(CHROM, POS, strain, REF, ALT, a1, a2, GT, everything()) %>%
      dplyr::arrange(CHROM, POS) 
  }
    tsv
  })
  results <- do.call(rbind, results)
  results
}

