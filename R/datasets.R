#' Strain Isotype (Dataset)
#'
#' Strain and isotype information.
#' 
#' @format tbl_df
#' @export
"strain_isotype"


#' SNP set (Dataset)
#'
#' Genotype information for C. elegans wild isolates. SNPs used from RADseq data set but 
#' includes all 124 wild isolates from whole-genome sequencing project
#'
#' @format A data frame with snps as rows and columns as strains:
#' \describe{
#'   \item{CHROM}{chromosome in roman numerals}
#'   \item{POS}{physical position in genome}
#'   \item{AB1}{first strain in collection, followed by all other strains}
#'   ...
#' }
#' @export
"snps"


#' Relatedness matrix (Dataset)
#'
#' Kinship matrix generated using the \code{a.mat} function in the rrBLUP package. 
#' Whole-genome SNP data was used to generate the relatedness of the strains in this file.
#' @format A Strain x Strain data frame for 152 strains
#' @export
"kinship"


#' Gene IDs (Dataset)
#'
#' A list of wormbase to standard gene identifier (e.g. _pot-2_) identifiers.
#' @format list
#' @export
"gene_ids"

#' VCF Version (Dataset)
#'
#' Current VCF Version.
#' @format vector
#' @export
"vcf_version"


#' Genome Build (Dataset)
#'
#' The wormbase build currently in use
#' @format A single element vector
#' @export
"wb_build"


#' Mapping SNPs (Dataset)
#'
#' Reduced SNP set
#' @format A data frame with snps as rows and columns as strains:
#' \describe{
#'   \item{CHROM}{chromosome in roman numerals}
#'   \item{POS}{physical position in genome}
#'   \item{AB1}{first strain in collection, followed by all other strains}
#'   ...
#' }
#' @export
"mapping_snps"

