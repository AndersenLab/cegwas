

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
    vcf_path <- paste0("https://storage.googleapis.com/elegansvariation.org/releases/", vcf_version, "/variation/WI.", vcf_version, ".soft-filter.vcf.gz")
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
    url <- paste0("https://storage.googleapis.com/cegwas/WS", wb_build , ".celegans_gff.db")
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
    region <- paste0("--regions ", region)
  } else {
    region <- ""
  }
  command <- paste0("bcftools view ", tag_snps, " ", region, " -m2 -M2 --min-af ", allele_freq, " ", vcf," | ",
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
  df <- vcf_to_matrix(vcf, region=region) %>%
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
