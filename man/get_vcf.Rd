% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf.R
\name{get_vcf}
\alias{get_vcf}
\title{Generate directory path to VCF file}
\usage{
get_vcf(remote = F, version = vcf_version)
}
\arguments{
\item{remote}{logical, to use remote VCF stored on google or local in Andersen lab dropbox}

\item{impute}{logical, to use imputed VCF file or non imputed}
}
\value{
character value corresponding to VCF location
}
\description{
\code{get_vcf} Generate directory path to VCF file
}
