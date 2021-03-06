% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwas_mappings.R
\name{cegwas_map}
\alias{cegwas_map}
\title{cegwas_map}
\usage{
cegwas_map(trait_data, cores = parallel::detectCores(),
  remove_strains = TRUE, kin_matrix = kinship, snpset = snps,
  duplicate_method = "first", BF = NA, mapping_snp_set = TRUE,
  snp_grouping = 200)
}
\arguments{
\item{cores}{number of cores on computer that you want to allocate for mapping. Default value is 4}

\item{remove_strains}{Remove strains with no known isotype. Default is FALSE.}

\item{kin_matrix}{is a strainXstrain matrix. default kinship matrix is described above.}

\item{duplicate_method}{Method for dealing with the presence of multiple strains falling into the same isotype. Either \code{"average"} to average phenotypes or \code{"first"} to take the first observation.}

\item{BF}{defines a custom bonferroni correction.}

\item{mapping_snp_set}{Use simulation based snps when TRUE. Use 5\% cut when FALSE.}

\item{snp_grouping}{defines grouping of peaks. Defined further below, default is 200.}

\item{data}{two element list. element 1 : traits. element 2: trait values with strains in columns
with each row corresponding to trait in element 1}

\item{snps}{is a set of mapping snps.}
}
\value{
Outputs a two element list that contains two dataframes. 
The first data frame is a processed mappings dataframe that contains the same columns
as the output of \code{\link{gwas_mappings}} with two additional columns. One that contains
the bonferroni corrected p-value (BF) and another that contains an identifier 1,0 if 
the indicated SNP has a higher -log10 value than the bonferroni cut off or not, respectively
The second data frame contains the variance explained data as well as all of the information from the first element.
}
\description{
\code{cegwas_map} is a convenience function takes trait data for a set of strains and performs
a mapping - returning a processed mapping data frame. \code{map} wraps \code{\link{process_pheno}}, \code{\link{gwas_mappings}}, 
and \code{\link{process_mappings}} into a single function.
}
