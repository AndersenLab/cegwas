#' SNP set.
#'
#' Genotype information for C. elegans wild isolates. SNPs used from RADseq data set but 
#' includes all 124 wild isolates from whole-genome sequencing project
#' @name snps
#' @format A data frame with snps as rows and columns as strains:
#' \describe{
#'   \item{CHROM}{chromosome in roman numerals}
#'   \item{POS}{physical position in genome}
#'   \item{AB1}{first strain in collection, followed by all other strains}
#'   ...
#' }
NULL

#' Relatedness matrix
#'
#' Kinship matrix generated using the \code{a.mat} function in the rrBLUP package. 
#' Whole-genome SNP data was used to generate the relatedness of the strains in this file.
#' @name kinship
#' @format A Strain x Strain data frame for 124 strains
NULL

#' Gene IDs
#'
#' Mapping between wormbase ids and common names.
#' @name gene_ids
#' @format list
NULL

