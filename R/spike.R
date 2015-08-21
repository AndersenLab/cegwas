#' Spike in an effect for GWAS mapping
#'
#' This function generates fake data for testing GWAS mappings, with the option
#' to generate QTL. 
#'
#' @param meadiff The abolute difference in means for the two allele populations
#' @param sd The standard deviation for both populations
#' @param snps The snp set being used
#' @param snpindex The index of the snp upon which to spike the effect
#' @return A phenotype data frame with the fake spiked data, ready for process_pheno
#' @export


spike <- function(meandiff, sd, snps, snpindex) {
    set.seed(17)
    
    strains <- colnames(snps)[-c(1, 2)]
    trait_values <- seq_along(strains)
    snps <- data.frame(snps)
    for(i in seq_along(strains)) {
        if(snps[snpindex, i + 2] == 0) {
            trait_values[i] <- rnorm(1, mean = meandiff, sd = sd)
        } else {
            trait_values[i] <- rnorm(1, mean = 0, sd = sd)
        }
    }
    
    trait_values <- rbind(trait_values)
    
    pheno <- data.frame("test", trait_values)
    pheno[,2:ncol(pheno)] <- lapply(pheno[,2:ncol(pheno)], as.numeric)
    
    colnames(pheno) <- c("trait", strains)
    return(pheno)
}