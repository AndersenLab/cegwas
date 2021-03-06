% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_mappings.R
\name{find_peaks}
\alias{find_peaks}
\title{Find Peaks from GWAS Peaks}
\usage{
find_peaks(processed_mapping_df, CI_size = 50, snp_grouping = 200)
}
\arguments{
\item{processed_mapping_df}{The first element of the output list from the \code{calculate_VE} function.}

\item{CI_size}{defines the size (in # SNPs) of confidence intervals. Default is 50 and is defined in more detail below.}

\item{snp_grouping}{defines grouping of peaks. Defined further below, default is 200.}
}
\value{
Outputs a two element list that contains. 
First element - data frame containing all identified intervals
Second element - list containing one element for each interval
}
\description{
\code{find_peaks} Identifies QTL from GWAS mapping data set.
}
\details{
This function identifies QTL by looking at SNPs above the bonferroni corrected p-value.
If only one SNP passed the significance cutoff, then the confidence interval is defined
as +/- \code{CI_size} (number of SNPs; default 50) away from that SNP. If multiple SNPs 
are above the cutoff, the function asks if SNPs are within an arbitrary number of SNPs away
\code{snp_grouping} - default 200. If the significant SNPs are within this range, they are grouped into the same peak.
If they are greater than this distance, then the peaks are considered unique.
}
