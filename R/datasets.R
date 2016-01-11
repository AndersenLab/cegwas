#' SNP set (Dataset)
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

#' Relatedness matrix (Dataset)
#'
#' Kinship matrix generated using the \code{a.mat} function in the rrBLUP package. 
#' Whole-genome SNP data was used to generate the relatedness of the strains in this file.
#' @name kinship
#' @format A Strain x Strain data frame for 152 strains
NULL


#' Strain Isotype (Dataset)
#'
#' Table of strains and their respective isotypes in addition to strains with potential problems.
#' @name strain_isotype
#' @format tbl_df
NULL


#' Gene IDs (Dataset)
#'
#' A list of wormbase to standard gene identifier (e.g. _pot-2_) identifiers.
#' @name gene_ids
#' @format tbl_df
NULL

#' Genome Build (Dataset)
#'
#' The wormbase build currently in use
#' @name wb_build
#' @format A single element vector
NULL
