% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant_correlation.R
\name{interval_summary}
\alias{interval_summary}
\title{Interval Summary}
\usage{
interval_summary(query, filter_variants = T, impute = F)
}
\arguments{
\item{query}{A query such as \code{pot-2}, \code{I:1-1000}, \code{WBGene00010195}.}

\item{filter_variants}{Filter out poor quality variants. Default is TRUE. Not compatible with imputed data.}

\item{impute}{Use imputed data. Default is FALSE.}
}
\value{
data frame of summarized information.
}
\description{
\code{interval_summary} Summarizes variation within an interval.
Can also be used to summarize regions by wormbase identifier or locus name.
}
\examples{
interval_summary("pot-1")
}
