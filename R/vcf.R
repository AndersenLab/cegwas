#' VCF To Matrix
#'
#' \code{vcf_to_matrix} converts variant calls from a bcf, vcf.gz or vcf file into an R dataframe. 
#' Uses only biallelic variants. Heterozygous calls are ignored. Requires bcftools. 
#'
#' @param vcf a bcf, vcf, or vcf.gz file
#' @param allele_freq allele frequency to filter on. Default is 0
#' @param variants A set of variants to filter on. By default, all variants are taken.
#' @return Matrix of genotype calls
#' @seealso \link{generate_kinship} \link{generate_mapping}
#' @export

vcf_to_matrix <- function(vcf, allele_freq = 0.0, tag_snps = NA) {
  gt_file <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = "")
  if (!is.na(tag_snps)) {
    tag_snps <- paste0("-T ", tag_snps)
  } else {
    tag_snps <- ""
  }
  command <- paste0("bcftools view ", tag_snps, " -m2 -M2 --min-af ", allele_freq, " ", vcf," | ",
                    "bcftools query --print-header -f '%CHROM\\t%POS[\\t%GT]\\n' | ",
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

generate_kinship <- function(vcf) {
  df <- vcf_to_matrix(vcf) %>%
    dplyr::select(-CHROM,-POS)
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