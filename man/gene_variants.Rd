% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{gene_variants}
\alias{gene_variants}
\title{Plot variants for gene}
\usage{
gene_variants(gene)
}
\arguments{
\item{gene}{is the gene of interest in gene name (e.g. "top-2", "pot-2") or wormbase gene ID (e.g. "WBGene00010785") format.}
}
\value{
Ouput is a list of ggplot objects with strains on the Y axis and variants for gene of interest on X axis. Tiles are colored by variant or reference call.
}
\description{
\code{gene_variants} generates a plot to visualize presence of variants for a particular gene of interest
}
\examples{
 gene_variants(gene = c("top-2","pot-2"))
 test[[1]]
test[[2]]
}
