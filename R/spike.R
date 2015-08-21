#' Spike in an effect for GWAS mapping
#'
#' This function generates fake data for testing GWAS mappings, with the option
#' to generate QTL. 
#'
#' @param snps The snp set being used
#' @param snpindex The index of the snp upon which to spike the effect
#' @return A phenotype data frame with the fake spiked data, ready for process_pheno
#' @export


spike <- function(snps, snpindex) {
    if (length(snpindex) == 1) {
        set.seed(17)
        
        strains <- colnames(snps)[-c(1, 2)]
        trait_values <- seq_along(strains)
        snps <- data.frame(snps)
        for(i in seq_along(strains)) {
            if(snps[snpindex, i + 2] == 0) {
                trait_values[i] <- rnorm(1, mean = 20, sd = 1)
            } else {
                trait_values[i] <- rnorm(1, mean = 0, sd = 1)
            }
        }
        
        trait_values <- rbind(trait_values)
        
        pheno <- data.frame("test", trait_values)
        pheno[,2:ncol(pheno)] <- lapply(pheno[,2:ncol(pheno)], as.numeric)
        
        colnames(pheno) <- c("trait", strains)
        return(pheno)
    } else {
        set.seed(17)
        strains <- colnames(snps)[-c(1, 2)]
        trait_values <- rep(0, times = length(strains))
        snps <- data.frame(snps)
        
        for (snp in snpindex) {
            meandiff <- rnorm(1, 30, sd = 5)
            sd <- rnorm(1, mean = 3, sd = 1)
            for(i in seq_along(strains)) {
                if(is.na(snps[snp, i + 2])){
                    trait_values[i] <- trait_values[i] + 0
                } else {
                    if(snps[snp, i + 2] == 0) {
                        trait_values[i] <- trait_values[i] + rnorm(1, mean = meandiff, sd = sd)
                    } else {
                        trait_values[i] <- trait_values[i] - rnorm(1, mean = meandiff, sd = sd)
                    }
                }
            }
        }
        
        trait_values <- rbind(trait_values)
        
        pheno <- data.frame("test", trait_values)
        pheno[,2:ncol(pheno)] <- lapply(pheno[,2:ncol(pheno)], as.numeric)
        
        colnames(pheno) <- c("trait", strains)
        return(pheno)
    }
}