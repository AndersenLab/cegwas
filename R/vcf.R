
#' Browse Variant Info
#'
#' \code{snpeff} enables you to query variants called and annotated by the \href{http://www.andersenlab.org}{Andersen Lab}. 
#' 
#' @param ... Gene names, regions, or wormbase identifiers to query.
#' @param severity A vector with variants of given severities (LOW, MODERATE, HIGH, MODIFIER). Default takes moderate and high. Use "ALL" to return all variants.
#' @param elements A vector containing gene structural elements (CDS, five_prime_UTR, exon, intron, three_prime_UTR). Use "ALL" to return all variants.
#' @param long Return dataset in long or wide format. Default is to return in long format.
#' @param remote Use remote data. Checks for local data if possible. False by default.
#' @param use custom vcf file.
#' @return Outputs a data frame that contains phenotype data, mapping data, and gene information for highly correlated variants in a particular QTL confidence interval.
#' @examples snpeff("pot-2","II:1-10000","WBGene00010785")
#' @export

snpeff <- function(...,
                   severity = c("HIGH","MODERATE"),
                   elements = c("exon"),
                   long = TRUE,
                   remote = FALSE,
                   vcf = NA) {
  
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
                         dplyr::distinct( .keep_all = TRUE))
      })) %>%
        dplyr::summarize(chrom = chrom[1], start = min(start), end = max(end)) %>%
        dplyr::mutate(region_format = paste0(chrom, ":", start, "-", end)) %>%
        dplyr::select(region_format) %>%
        dplyr::distinct(.keep_all = TRUE))$region_format, collapse = ",")
      
      query_type <- "locus"
      
      if (stringr::str_length(regions[[1]]) == 0) {
        message(paste0(query, " not found."))
        region <- NA
      }
    } else {
      region <- query
      query_type <- "region"
    }
    
    if(is.na(vcf)) {
      vcf_path <- get_vcf(remote = remote)
    } else {
      vcf_path <- vcf
      gene_ids <- NA
    }
    
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
          dplyr::select(dplyr::one_of(c(base_header, ANN_header)), dplyr::everything(), -extra) %>%
          #dplyr::mutate(gene_name = as.character(gene_ids[gene_name])) %>%
          dplyr::mutate(query = query, region = region) %>%
          dplyr::select(CHROM, POS, query, region, dplyr::everything())
        
        # For locus queries, filter out non-matching genes.
        if (query_type == 'locus') {
          tsv <- dplyr::filter(tsv, (query == gene_id) | (query == gene_name) | (query == feature_id) )
        }
        
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
            dplyr::select(CHROM, POS, strain, REF, ALT, a1, a2, GT, FT, FILTER, DP, DP4, SP, HP, dplyr::everything()) %>%
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

get_vcf <- function(remote = F, version = vcf_version) {
  
  # Set vcf path; determine whether local or remote
  vcf_path <- paste0("~/Dropbox/Andersenlab/Reagents/WormReagents/_SEQ/WI/WI-", vcf_version, "/vcf/WI.", vcf_version, ".snpeff.vcf.gz")
  # Use remote if not available.
  local_or_remote <- "locally"
  if (!file.exists(vcf_path) | remote == T) {
    vcf_path <- paste0("http://storage.googleapis.com/elegansvariation.org/releases/", vcf_version, "/WI.", vcf_version, ".soft-filtered.vcf.gz")
    message("Using remote vcf")
  } else {
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
                                     dplyr::distinct(.keep_all = TRUE))$biotype
  if (!(id_type %in% valid_id_types)) {
    message("Available ID Types:")
    cat(valid_id_types, sep = "\n")
  } else {
    (dplyr::collect(elegans_gff %>%
                      dplyr::filter(biotype == id_type, type_of == "exon") %>%
                      dplyr::select(gene_id)))$gene_id
  }
  
}


#' VCF To Matrix
#'
#' \code{vcf_to_matrix} converts variant calls from a bcf, vcf.gz or vcf file into an R dataframe. 
#' Uses only biallelic variants. Heterozygous calls are ignored. Requires bcftools. 
#'
#' @param vcf a bcf, vcf, or vcf.gz file
#' @param allele_freq allele frequency to filter on. Default is 0
#' @param tag_snips A set of variants to filter on. By default, all variants are taken.
#' @param region A region to subset on.
#' @return Matrix of genotype calls
#' @seealso \link{generate_kinship} \link{generate_mapping}
#' @export

vcf_to_matrix <- function(vcf, allele_freq = 0.0, tag_snps = NA, region = NA) {
  gt_file <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = "")
  if (!is.na(tag_snps)) {
    tag_snps <- paste0("-T ", tag_snps)
  } else {
    tag_snps <- ""
  }
  if (!is.na(region)) {
    region <- paste0("--region", region)
  } else {
    region <- "" 
  }
  command <- paste0("bcftools view ", tag_snps, " -m2 -M2 --min-af ", allele_freq, " ", vcf," | ",
                    region,
                    "bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' | ",
                    "sed 's/[[# 0-9]*\\]//g' | ", 
                    "sed 's/:GT//g' | ",         
                    "sed 's/0|0/-1/g'   | ",
                    "sed 's/1|1/1/g'   | ",
                    "sed 's/0|1/NA/g'  | ",
                    "sed 's/1|0/NA/g'  | ",
                    "sed 's/.|./NA/g'  | ",
                    "sed 's/0\\/0/-1/g'  | ",
                    "sed 's/1\\/1/1/g'  | ",
                    "sed 's/0\\/1/NA/g' | ",
                    "sed 's/1\\/0/NA/g' | ",
                    "sed 's/.\\/./NA/g'  > ", gt_file)
  cat(command)
  # Generate simplified representation of genotypes
  system(command, intern = T)
  df <- dplyr::tbl_df(data.table::fread(gt_file))
  
  df
}


#' Generate Kinship
#'
#' \code{generate_kinship} uses \code{\link{vcf_to_matrix}} and the \code{\link[rrBLUP]{A.mat}} function from \code{\link{rrBLUP}}.
#' This function is used to generate a kinship matrix for mapping in conjunction with \code{\link{gwas_mappings}}.
#'
#' @param a vcf file.
#' @return A kinship matrix. 
#' @seealso \link{vcf_to_matrix} \link{generate_kinship}
#' @export

generate_kinship <- function(vcf, region=NA) {
  df <- vcf_to_matrix(vcf, region) %>%
    dplyr::select(-CHROM,-POS, -REF, -ALT)
  rrBLUP::A.mat(t(data.matrix(df)))
}

#' Generate Mapping variant set
#'
#' \code{generate_mapping} generates a dataframe suitable for mapping. This function uses bcftools to filter out variants with low allele frequencies <5%. This can be customized by setting
#' the af parameter. Additionally, this function uses a subset of variants from the \emph{C. elegans} genome that tag haplotype blocks. 
#' Tag SNPs were generated in \emph{\href{https://dx.doi.org/10.1038/ng.1050}{Andersen 2012 et al}}. Tag snps used can
#' be adjusted using a text file that specifies \emph{CHROM}    \emph{POS}. A custom set of variants can be used with the \code{variants} parameter.
#'
#' @return Mapping Snpset
#' @param vcf a bcf, vcf.gz, or vcf file
#' @param allele_freq minimum allele frequency. Default is >= 5\%
#' @param tag_snps A set of variants to filter on. By default, a set of tag snps is used 
#' @seealso \link{vcf_to_matrix} \link{generate_mapping}
#' @export

generate_mapping <- function(vcf, allele_freq = 0.00, tag_snps = paste0(path.package("cegwas"),"/41188.WS245.txt.gz")) {
  vcf_to_matrix(vcf, allele_freq = allele_freq, tag_snps = tag_snps) 
}
