#' Calculate Variance Explained for Significant SNPs
#'
#' \code{calculate_VE} calculates the variance explained (VE) for significant SNPs by using the 
#' the spearman rank correlation coefficient.
#'
#' This function requires three inputs, two of which are provided by the user and the other is loaded by the package.
#'
#' @param mapping_df the output from the \code{gwas_mappings} function. User input
#' @param phenotype_df two element list. element 1 : traits. element 2: trait values with strains in columns
#' with each row corresponding to trait in element 1
#' @return Outputs a two element list that contains two dataframes. 
#' The first data frame is a processed mappings dataframe that contains the same columns
#' as the output of \code{gwas_mappings} with two additional columns. One that contains
#' the bonferroni corrected p-value (BF) and another that contains an identifier 1,0 if 
#' the indicated SNP has a higher -log10 value than the bonferroni cut off or not, respectively
#' The second data frame contains the variance explained data as well as all of the information from the first element.
#' @export

calculate_VE <- function( mapping_df,
                          phenotype_df ) {
  
  # mapping_df <- sig_maps
  # snp_df <- snps
  # phenotype_df <- test
  
    # format phenotypes
    pheno <- phenotype_df[[2]]
    row.names(pheno) <- phenotype_df[[1]]
    pheno$trait <- phenotype_df[[1]]
  
  
  Processed <- mapping_df %>%
    dplyr::group_by( trait ) %>%
    dplyr::mutate( BF = -log10(.05/n()) ) %>% #  add BF threshold
    dplyr::group_by( trait ) %>%
    dplyr::mutate( aboveBF = ifelse(log10p >= BF, 1, 0) ) %>% #  label SNPs as significant
    dplyr::filter(sum(aboveBF) > 0) %>% # keep only significant mappings
    dplyr::ungroup()
  
  ## Select SNPs above BF
  snpsForVE <- Processed %>%
    dplyr::filter( aboveBF == 1 ) %>%
    dplyr::select( marker, trait )
  
  snpsForVE$trait <- as.character(snpsForVE$trait)
  
  ## Trim raw data and snp set to contain phenotypes and snps from BF mappings
  
  row.names(pheno) <- gsub("-", "\\.", row.names(pheno))
  pheno$trait <- gsub("-", "\\.", pheno$trait) 
  
  # Trim phenotypes and join to significant snps identified in mapping
  rawTr <- pheno[ row.names(pheno) %in% snpsForVE$trait,] %>%
    tidyr::gather( strain, value, -trait ) %>% # make long format
    dplyr::left_join( ., snpsForVE, 
                      by = "trait" ) # join to significant SNPs from mapping dataframe
  
  # make columns factored by dplyr into characters to minimize warnings
  rawTr$marker <- as.character(rawTr$marker)
  rawTr$strain <- as.character(rawTr$strain)
  
  # Trim snps to only contain those that are significant from mappings
  
  snp_df <- data.frame(snps)
  
  gINFO <- snp_df %>%
    dplyr::mutate( marker = paste(CHROM, POS, sep = "_")) %>%
    dplyr::filter( marker %in% snpsForVE$marker ) %>%
    tidyr::gather( strain, allele, -marker, -CHROM, -POS) 
    
  # make columns factored by dplyr into characters to minimize warnings
  gINFO$marker <- as.character(gINFO$marker)
  
  # combine genotype data, phenotype data, and significant snps from mappnings
  gINFO <- data.frame(gINFO) %>%
    dplyr::left_join( ., snpsForVE, by= "marker") %>% # join significant snps with genotypes
    dplyr::left_join( rawTr, ., by=c( "trait", "strain", "marker") ) # join to phenotypes
  
  # calculate variance explained using spearman correlation
  cors <- gINFO %>%
    # each significant snp contains genotype and phenotype information for all strains
    # so group by both to calculate variance explained for each significant snp
    dplyr::group_by( trait, marker ) %>% 
    # calculate correlation
    dplyr::mutate( var.exp = cor(value, allele, use = "pairwise.complete.obs", method="spearman")^2 )  
  
  # bring it all together, that :
  # # # genotypes
  # # # phenotypes
  # # # correlations
  # # # # # # # # # # # FOR ALL SIGNIFICANT SNPS
  
  CORmaps <- Processed %>%
    dplyr::left_join( ., cors, by=c("trait","marker","CHROM","POS"), copy=TRUE )
  
  return( list(Processed, CORmaps) )
}

#' Find Peaks from GWAS Peaks
#'
#' \code{find_peaks} Identifies QTL from GWAS mapping data set.
#'
#' This function identifies QTL by looking at SNPs above the bonferroni corrected p-value.
#' If only one SNP passed the significance cutoff, then the confidence interval is defined
#' as +/- \code{CI_size} (number of SNPs; default 50) away from that SNP. If multiple SNPs 
#' are above the cutoff, the function asks if SNPs are within an arbitrary number of SNPs away
#' \code{snp_grouping} - default 200. If the significant SNPs are within this range, they are grouped into the same peak.
#' If they are greater than this distance, then the peaks are considered unique. 
#' 
#' @param processed_mapping_df The first element of the output list from the \code{calculate_VE} function.
#' @param CI_size defines the size (in # SNPs) of confidence intervals. Default is 50 and is defined in more detail below.
#' @param snp_grouping defines grouping of peaks. Defined further below, default is 200.
#' @return Outputs a two element list that contains. 
#' First element - data frame containing all identified intervals
#' Second element - list containing one element for each interval
#' @export

find_peaks <- function( processed_mapping_df, 
                        CI_size = 50,
                        snp_grouping = 200 ) {
  
  # processed_mapping_df <- Processed
  # snp_grouping <- 200
  # CI_size <- 50
  
  # # # IDENTIFY PEAKS
  
  # PHENOTYPES THAT HAVE SIGNIFICANT MAPPINGS
  phenotypes <- as.character(unique(processed_mapping_df$trait))
  
  # INITIALIZE A LIST TO PUT INDIVIDUAL PHENOTYPE INTERVAL INFORMATION
  intervals <- list()
  
  # LOOP THROUGH ALL UNIQUE PHENOTYPES
  for( i in 1:length(phenotypes) ){
    
    print(paste(100*signif(i/length(phenotypes),3), "%",sep=""))
    
    # PREP DATA FRAME FOR PEAK IDENTIFICATION
    
    PeakDF <- processed_mapping_df %>%
      dplyr::filter( trait == phenotypes[i] ) %>%
      dplyr::group_by( CHROM, trait ) %>%
      dplyr::mutate( index = 1:n() ) %>% # SNP INDEX
      dplyr::mutate( peaks = cumsum(aboveBF) ) %>% # IDENTIFY PEAKS
      dplyr::filter( aboveBF == 1 )%>% # KEEP SNPS ABOVE BONFERRONI
      dplyr::group_by( CHROM, trait) %>%
      dplyr::mutate( nBF = n() ) %>% # COUNT NUMBER OF SNPS ABOVE BONFERRONI PER PHENOTYPE PER CHROMOSOME
      dplyr::group_by( CHROM, trait ) %>% 
      dplyr::arrange( CHROM, POS ) # ARRANGE DATA BY CHROMOSOME AND POSITION
    
    # generate a SNP index for SNPs on each chromosome
    SNPindex <- processed_mapping_df %>%
      dplyr::filter( trait == phenotypes[i] ) %>%
      dplyr::group_by( CHROM, trait ) %>%
      dplyr::mutate( index = 1:n() )%>%
      dplyr::distinct( CHROM, POS )%>%
      dplyr::select( CHROM, POS, index )%>%
      dplyr::filter( POS == min(POS) | POS == max(POS) )
    
    # FILTER COMPLETE DATA SET TO JUST LOOK AT ONE PHENOTYPE AT A TIME
    findPks <- PeakDF %>%
      dplyr::filter( trait == phenotypes[i] ) %>% 
      dplyr::group_by( CHROM ) %>% 
      dplyr::arrange( CHROM, POS )
    
    # IF ONLY ONE SNP PASSED SIGNIFICANCE THRESHOLD LABEL PEAK ID AS 1
    if ( findPks$nBF == 1 & length(unique(findPks$CHROM) ) == 1 ){
      
      findPks$pID <- 1
      
      # PLUS / MINUS 50 SNPS FROM PEAK SNP DEFINES CONFIDENCE INTERVAL
      findPks <- findPks %>%
        dplyr::group_by( CHROM, pID, trait ) %>%
        dplyr::mutate( start = min(index) - CI_size,
                       end = max(index) + CI_size )
      
      for( k in 1:nrow(findPks) ){
        
        tSNPs <- SNPindex %>%
          dplyr::filter( CHROM == findPks$CHROM[k] )
        
        if( findPks$start[k] < min(tSNPs$index) ){
          
          findPks$start[k] <- min(tSNPs$index)
          
        } else if( findPks$end[k] > max(tSNPs$index) ) {
          
          findPks$end[k] <- max(tSNPs$index)
          
        }
      }
      
      # APPEND TO LIST
      intervals[[i]] <- findPks %>%
        dplyr::ungroup()
      
    } 
    else 
    {
      # INITIALIZE PEAK ID COLUMN WITH 1'S :: GIVES YOU A STARTING POINT
      findPks$pID <- 1
      
      # LOOP THROUGH ROWS FOR EACH PHENOTYPE CORRESPONDING TO SNPS ABOVE BONFERRONI CORRECTION
      # START AT ROW 2 BECAUSE THERE WILL ALWAYS BE AT LEAST 1 UNIQUE PEAK
      for( j in 2:nrow(findPks) ){
        
        # IF 
        # SNP INDEX IS WITHIN A CERTAIN RANGE (snp_grouping) OF SNP FROM PREVIOUS ROW
        # AND
        # ON THE SAME CHROMOSOME AS SNP FROM PREVIOUS ROW
        # # # # CONSIDER THEM TO BE THE SAME PEAK
        
        # IF THE ABOVE CONDITIONS ARE NOT MET
        # ADD 1 TO THE PEAK ID (i.e. IDENTIFY AS A NEW PEAK)
        findPks$pID[j] <- ifelse( abs(findPks$index[j] - findPks$index[j-1]) < snp_grouping &
                                    findPks$CHROM[j] == findPks$CHROM[j-1],
                                  findPks$pID[j-1],
                                  findPks$pID[j-1]+1)
      }
      
      # PLUS / MINUS 50 SNPS FROM PEAK SNP DEFINES CONFIDENCE INTERVAL
      findPks <- findPks %>%
        dplyr::group_by( CHROM , pID, trait) %>%
        dplyr::mutate(start = min(index) - CI_size,
                      end = max(index) + CI_size)
      
      
      for( k in 1:nrow(findPks) ){
        
        tSNPs <- SNPindex %>%
          dplyr::filter( CHROM == findPks$CHROM[k] )
        
        if( findPks$start[k] < min(tSNPs$index) ){
          
          findPks$start[k] <- min(tSNPs$index)
          
        } else if( findPks$end[k] > max(tSNPs$index) ) {
          
          findPks$end[k] <- max(tSNPs$index)
          
        }
      } 
      
      
    }
    
    # APPEND TO LIST
    intervals[[i]] <- findPks %>% 
      dplyr::ungroup()
    
  }
  # BIND GENERATED LIST TOGETHER 
  intervalDF <- data.table::rbindlist(intervals)
  
  return( list(intervalDF, intervals) )
}

#' Identify confidence intervals associated with QTL. 
#'
#' \code{identify_CI} Identifies confidence intervals for identified QTL
#'
#' Function to combine all of the previously generated data into one data frame. Converts SNP index confidence intervals
#' into genomic position confidence intervals.
#' 
#' @param processed_mapping_df The first element of the \code{calculate_VE} function output
#' @param peak_df The first element of the \code{find_peaks} function output
#' @param peak_list The second element of the \code{find_peaks} function output 
#' @param correlation_df The second element of the \code{calculate_VE} function output
#' @return Outputs processed mapping dataframe that contains original mapping dataframe with appended information for significant SNPs only, including:
#' variance explained, confidence interval information, genotype information
#' @export

identify_CI <- function( processed_mapping_df, 
                         peak_df, 
                         peak_list, 
                         correlation_df ) {
  
  # processed_mapping_df <- Processed
  # peak_df = intervalDF 
  # peak_list = intervals
  # correlation_df = CORmaps 
  
  # FILTER COMPLETE MAPPING SET TO ONLY CONTAIN INTERVAL INDICIES TO SAVE COMPUTATIONAL TIME BELOW
  Pos_Index_Reference  <- processed_mapping_df %>%
    dplyr::group_by( CHROM, trait ) %>%
    dplyr::mutate( index = 1:n() ) %>%
    dplyr::mutate( peaks = cumsum(aboveBF) ) %>%
    dplyr::select( trait, CHROM, POS, index ) %>%
    dplyr::filter( index %in% c(unique(peak_df$start), unique(peak_df$end)) ) %>%
    dplyr::ungroup()
  
  Pos_Index_Reference$trait <- as.character(Pos_Index_Reference$trait)
  
  # INITIALIZE LIST TO APPEND INTERVAL POSITION DATA FOR EACH PHENOTYPE
  interval_positions <- list()
  
  # LOOP THROUGH UNIQUE PHENOTYPES TO LINK CONFIDENCE INTERVALS IN INDEX FORM TO POSITION FORM
  for( i in 1:length(peak_list)){
    
    print(paste(100*signif(i/length(peak_list),3), "%",sep=""))
    
    peak_list[[i]]$trait <- as.character(peak_list[[i]]$trait)
    
    peak_list[[i]] <- dplyr::distinct(peak_list[[i]], pID)
    
    # FILTER TO LOOK AT ONE PHENOTYPE AT A TIME
    # FILTER APPROPRIATE INTERVAL INDICIES AND CHROMOSOMES FOR THAT PHENOTYPE
    trait_i <- unique(peak_list[[i]]$trait)
    index_i <- c(peak_list[[i]]$start, peak_list[[i]]$end) 
    CHROM_i <- peak_list[[i]]$CHROM
    
    PKpos <- data.frame(Pos_Index_Reference) %>%
      dplyr::filter(trait == trait_i &
                    index %in% index_i &
                    CHROM %in%  CHROM_i) %>%
      # JOIN POSITION INFORMATION TO PHENOTYPE PEAK INFORMATION
      dplyr::left_join( ., peak_list[[i]], by= c("trait","CHROM") )%>%
      # YOU WILL GET UNWANTED SNP INDEX INFORMATION IN SITUATIONS WHERE YOU HAVE MULTIPLE PEAKS 
      # ELIMINATE THOSE BY MATCHING START AND END FROM INDEX DATAFRAME TO INDEX FROM POSITION DATAFRAME
      # FIRST FLAG
      dplyr::mutate(issues = ifelse(start == index.x | end == index.x, 1, 0))%>%
      # THEN REMOVE
      dplyr::filter(issues != 0)%>%
      # SELECT COLUMNS OF INTEREST
      dplyr::select(trait, CHROM, POS.x, POS.y, pID, log10p, index.x, index.y, start, end)%>%
      # GROUP BY PEAK IDS ORIGINALLY PRESENT IN INDEX DATAFRAME
      dplyr::group_by(CHROM, pID) %>%
      # GENERATE COLUMNS TO WITH INTERVAL POSITIONS AND PEAK POSITIONS
      dplyr::mutate(startPOS = min(POS.x),
                    peakPOS = POS.y,
                    endPOS = max(POS.x)) %>%
      # ELIMINATE REDUNDANT DATA
      dplyr::distinct(trait, CHROM, pID, peakPOS) %>%
      # SELECT COLUMNS OF NTEREST
      dplyr::select(trait, CHROM, POS = POS.y, startPOS, peakPOS, endPOS, peak_id = pID)
    
    # APPEND TO LIST
    interval_positions[[i]] <- PKpos
  }
  
  # BIND EVERYTHING
  interval_pos_df <- data.frame(data.table::rbindlist(interval_positions)) %>%
    # CALCULATE INTERVAL SIZE
    dplyr::mutate(interval_size = endPOS - startPOS)
  
  # JOIN INTERVAL POSITIONS TO DATA FRAME CONTAINING CORRELATION INFORMATION AND PHENOTYPE INFORMATION
  Final_Processed_Mappings <- dplyr::left_join( correlation_df, interval_pos_df, 
                                         by = c("trait", "CHROM", "POS"),
                                         copy = TRUE )
  
  return(Final_Processed_Mappings)
}

#' Fully process GWAS mapping output
#'
#' \code{process_mappings} takes \code{gwas_mappings} output and calculates variance explained and
#' identifies peaks and associated confidence intervals
#'
#' This function combines \code{calculate_VE}, \code{find_peaks}, and \code{identify_CI} into one function when intermediate dataframes are not wanted
#' 
#' @param mapping_df Output from \code{gwas_mappings} function
#' @param phenotype_df phenotype data frame generated by \code{process_pheno}. two element list. element 1 : traits. element 2: trait values with strains in columns
#' with each row corresponding to trait in element 1
#' @param CI_size defines the size (in # SNPs) of confidence intervals. Default is 50 and is defined in more detail below.
#' @param snp_grouping defines grouping of peaks. Defined further below, default is 200.
#' @return Outputs processed mapping dataframe that contains original mapping dataframe with appended information for significant SNPs only, including:
#' variance explained, confidence interval information, genotype information
#' @export


process_mappings <- function(mapping_df,
                             phenotype_df,
                             CI_size = 50,
                             snp_grouping = 200){
  
  # format phenotypes
  pheno <- phenotype_df[[2]]
  row.names(pheno) <- phenotype_df[[1]]
  pheno$trait <- phenotype_df[[1]]
  
  
  Processed <- mapping_df %>%
    dplyr::group_by( trait ) %>%
    dplyr::mutate( BF = -log10(.05/n()) ) %>% #  add BF threshold
    dplyr::group_by( trait ) %>%
    dplyr::mutate( aboveBF = ifelse(log10p >= BF, 1, 0) ) %>% #  label SNPs as significant
    dplyr::filter(sum(aboveBF) > 0) %>% # keep only significant mappings
    dplyr::ungroup()
  
  ## Select SNPs above BF
  snpsForVE <- Processed %>%
    dplyr::filter( aboveBF == 1 ) %>%
    dplyr::select( marker, trait )
  
  snpsForVE$trait <- as.character(snpsForVE$trait)
  
  ## Trim raw data and snp set to contain phenotypes and snps from BF mappings
  
  row.names(pheno) <- gsub("-", "\\.", row.names(pheno))
  pheno$trait <- gsub("-", "\\.", pheno$trait) 
  
  # Trim phenotypes and join to significant snps identified in mapping
  rawTr <- pheno[ row.names(pheno) %in% snpsForVE$trait,] %>%
    tidyr::gather( strain, value, -trait ) %>% # make long format
    dplyr::left_join( ., snpsForVE, 
                      by = "trait" ) # join to significant SNPs from mapping dataframe
  
  # make columns factored by dplyr into characters to minimize warnings
  rawTr$marker <- as.character(rawTr$marker)
  rawTr$strain <- as.character(rawTr$strain)
  
  # Trim snps to only contain those that are significant from mappings
  
  snp_df <- data.frame(snps)
  
  gINFO <- snp_df %>%
    dplyr::mutate( marker = paste(CHROM, POS, sep = "_")) %>%
    dplyr::filter( marker %in% snpsForVE$marker ) %>%
    tidyr::gather( strain, allele, -marker, -CHROM, -POS) 
  
  # make columns factored by dplyr into characters to minimize warnings
  gINFO$marker <- as.character(gINFO$marker)
  
  # combine genotype data, phenotype data, and significant snps from mappnings
  gINFO <- data.frame(gINFO) %>%
    dplyr::left_join( ., snpsForVE, by= "marker") %>% # join significant snps with genotypes
    dplyr::left_join( rawTr, ., by=c( "trait", "strain", "marker") ) # join to phenotypes
  
  # calculate variance explained using spearman correlation
  cors <- gINFO %>%
    # each significant snp contains genotype and phenotype information for all strains
    # so group by both to calculate variance explained for each significant snp
    dplyr::group_by( trait, marker ) %>% 
    # calculate correlation
    dplyr::mutate( var.exp = cor(value, allele, use = "pairwise.complete.obs", method="spearman")^2 )  
  
  # bring it all together, that :
  # # # genotypes
  # # # phenotypes
  # # # correlations
  # # # # # # # # # # # FOR ALL SIGNIFICANT SNPS
  
  CORmaps <- Processed %>%
    dplyr::left_join( ., cors, by=c("trait","marker","CHROM","POS"), copy=TRUE )
  
  processed_mapping_df <- Processed
  correlation_df <- CORmaps
  
  # # # Part 2
  
  
  # # # IDENTIFY PEAKS
  
  # PHENOTYPES THAT HAVE SIGNIFICANT MAPPINGS
  phenotypes <- as.character(unique(processed_mapping_df$trait))
  
  # INITIALIZE A LIST TO PUT INDIVIDUAL PHENOTYPE INTERVAL INFORMATION
  intervals <- list()
  
  # LOOP THROUGH ALL UNIQUE PHENOTYPES
  for( i in 1:length(phenotypes) ){
    
    print(paste(100*signif(i/length(phenotypes),3), "%",sep=""))
    
    # PREP DATA FRAME FOR PEAK IDENTIFICATION
    
    PeakDF <- processed_mapping_df %>%
      dplyr::filter( trait == phenotypes[i] ) %>%
      dplyr::group_by( CHROM, trait ) %>%
      dplyr::mutate( index = 1:n() ) %>% # SNP INDEX
      dplyr::mutate( peaks = cumsum(aboveBF) ) %>% # IDENTIFY PEAKS
      dplyr::filter( aboveBF == 1 )%>% # KEEP SNPS ABOVE BONFERRONI
      dplyr::group_by( CHROM, trait) %>%
      dplyr::mutate( nBF = n() ) %>% # COUNT NUMBER OF SNPS ABOVE BONFERRONI PER PHENOTYPE PER CHROMOSOME
      dplyr::group_by( CHROM, trait ) %>% 
      dplyr::arrange( CHROM, POS ) # ARRANGE DATA BY CHROMOSOME AND POSITION
    
    # generate a SNP index for SNPs on each chromosome
    SNPindex <- processed_mapping_df %>%
      dplyr::filter( trait == phenotypes[i] ) %>%
      dplyr::group_by( CHROM, trait ) %>%
      dplyr::mutate( index = 1:n() )%>%
      dplyr::distinct( CHROM, POS )%>%
      dplyr::select( CHROM, POS, index )%>%
      dplyr::filter( POS == min(POS) | POS == max(POS) )
    
    # FILTER COMPLETE DATA SET TO JUST LOOK AT ONE PHENOTYPE AT A TIME
    findPks <- PeakDF %>%
      dplyr::filter( trait == phenotypes[i] ) %>% 
      dplyr::group_by( CHROM ) %>% 
      dplyr::arrange( CHROM, POS )
    
    # IF ONLY ONE SNP PASSED SIGNIFICANCE THRESHOLD LABEL PEAK ID AS 1
    if ( findPks$nBF == 1 & length(unique(findPks$CHROM) ) == 1 ){
      
      findPks$pID <- 1
      
      # PLUS / MINUS 50 SNPS FROM PEAK SNP DEFINES CONFIDENCE INTERVAL
      findPks <- findPks %>%
        dplyr::group_by( CHROM, pID, trait ) %>%
        dplyr::mutate( start = min(index) - CI_size,
                       end = max(index) + CI_size )
      
      for( k in 1:nrow(findPks) ){
        
        tSNPs <- SNPindex %>%
          dplyr::filter( CHROM == findPks$CHROM[k] )
        
        if( findPks$start[k] < min(tSNPs$index) ){
          
          findPks$start[k] <- min(tSNPs$index)
          
        } else if( findPks$end[k] > max(tSNPs$index) ) {
          
          findPks$end[k] <- max(tSNPs$index)
          
        }
      }
      
      # APPEND TO LIST
      intervals[[i]] <- findPks %>%
        dplyr::ungroup()
      
    } 
    else 
    {
      # INITIALIZE PEAK ID COLUMN WITH 1'S :: GIVES YOU A STARTING POINT
      findPks$pID <- 1
      
      # LOOP THROUGH ROWS FOR EACH PHENOTYPE CORRESPONDING TO SNPS ABOVE BONFERRONI CORRECTION
      # START AT ROW 2 BECAUSE THERE WILL ALWAYS BE AT LEAST 1 UNIQUE PEAK
      for( j in 2:nrow(findPks) ){
        
        # IF 
        # SNP INDEX IS WITHIN A CERTAIN RANGE (snp_grouping) OF SNP FROM PREVIOUS ROW
        # AND
        # ON THE SAME CHROMOSOME AS SNP FROM PREVIOUS ROW
        # # # # CONSIDER THEM TO BE THE SAME PEAK
        
        # IF THE ABOVE CONDITIONS ARE NOT MET
        # ADD 1 TO THE PEAK ID (i.e. IDENTIFY AS A NEW PEAK)
        findPks$pID[j] <- ifelse( abs(findPks$index[j] - findPks$index[j-1]) < snp_grouping &
                                    findPks$CHROM[j] == findPks$CHROM[j-1],
                                  findPks$pID[j-1],
                                  findPks$pID[j-1]+1)
      }
      
      # PLUS / MINUS 50 SNPS FROM PEAK SNP DEFINES CONFIDENCE INTERVAL
      findPks <- findPks %>%
        dplyr::group_by( CHROM , pID, trait) %>%
        dplyr::mutate(start = min(index) - CI_size,
                      end = max(index) + CI_size)
      
      
      for( k in 1:nrow(findPks) ){
        
        tSNPs <- SNPindex %>%
          dplyr::filter( CHROM == findPks$CHROM[k] )
        
        if( findPks$start[k] < min(tSNPs$index) ){
          
          findPks$start[k] <- min(tSNPs$index)
          
        } else if( findPks$end[k] > max(tSNPs$index) ) {
          
          findPks$end[k] <- max(tSNPs$index)
          
        }
      } 
      
      
    }
    
    # APPEND TO LIST
    intervals[[i]] <- findPks %>% 
      dplyr::ungroup()
    
  }
  # BIND GENERATED LIST TOGETHER 
  intervalDF <- data.table::rbindlist(intervals)
  
  peak_df <- intervalDF
  peak_list <- intervals
  
  # FILTER COMPLETE MAPPING SET TO ONLY CONTAIN INTERVAL INDICIES TO SAVE COMPUTATIONAL TIME BELOW
  Pos_Index_Reference  <- processed_mapping_df %>%
    dplyr::group_by( CHROM, trait ) %>%
    dplyr::mutate( index = 1:n() ) %>%
    dplyr::mutate( peaks = cumsum(aboveBF) ) %>%
    dplyr::select( trait, CHROM, POS, index ) %>%
    dplyr::filter( index %in% c(unique(peak_df$start), unique(peak_df$end)) ) %>%
    dplyr::ungroup()
  
  Pos_Index_Reference$trait <- as.character(Pos_Index_Reference$trait)
  
  # INITIALIZE LIST TO APPEND INTERVAL POSITION DATA FOR EACH PHENOTYPE
  interval_positions <- list()
  
  # LOOP THROUGH UNIQUE PHENOTYPES TO LINK CONFIDENCE INTERVALS IN INDEX FORM TO POSITION FORM
  for( i in 1:length(peak_list)){
    
    print(paste(100*signif(i/length(peak_list),3), "%",sep=""))
    
    peak_list[[i]]$trait <- as.character(peak_list[[i]]$trait)
    
    peak_list[[i]] <- dplyr::distinct(peak_list[[i]], pID)
    
    # FILTER TO LOOK AT ONE PHENOTYPE AT A TIME
    # FILTER APPROPRIATE INTERVAL INDICIES AND CHROMOSOMES FOR THAT PHENOTYPE
    trait_i <- unique(peak_list[[i]]$trait)
    index_i <- c(peak_list[[i]]$start, peak_list[[i]]$end) 
    CHROM_i <- peak_list[[i]]$CHROM
    
    PKpos <- data.frame(Pos_Index_Reference) %>%
      dplyr::filter(trait == trait_i &
                      index %in% index_i &
                      CHROM %in%  CHROM_i) %>%
      # JOIN POSITION INFORMATION TO PHENOTYPE PEAK INFORMATION
      dplyr::left_join( ., peak_list[[i]], by= c("trait","CHROM") )%>%
      # YOU WILL GET UNWANTED SNP INDEX INFORMATION IN SITUATIONS WHERE YOU HAVE MULTIPLE PEAKS 
      # ELIMINATE THOSE BY MATCHING START AND END FROM INDEX DATAFRAME TO INDEX FROM POSITION DATAFRAME
      # FIRST FLAG
      dplyr::mutate(issues = ifelse(start == index.x | end == index.x, 1, 0))%>%
      # THEN REMOVE
      dplyr::filter(issues != 0)%>%
      # SELECT COLUMNS OF INTEREST
      dplyr::select(trait, CHROM, POS.x, POS.y, pID, log10p, index.x, index.y, start, end)%>%
      # GROUP BY PEAK IDS ORIGINALLY PRESENT IN INDEX DATAFRAME
      dplyr::group_by(CHROM, pID) %>%
      # GENERATE COLUMNS TO WITH INTERVAL POSITIONS AND PEAK POSITIONS
      dplyr::mutate(startPOS = min(POS.x),
                    peakPOS = POS.y,
                    endPOS = max(POS.x)) %>%
      # ELIMINATE REDUNDANT DATA
      dplyr::distinct(trait, CHROM, pID, peakPOS) %>%
      # SELECT COLUMNS OF NTEREST
      dplyr::select(trait, CHROM, POS = POS.y, startPOS, peakPOS, endPOS, peak_id = pID)
    
    # APPEND TO LIST
    interval_positions[[i]] <- PKpos
  }
  
  # BIND EVERYTHING
  interval_pos_df <- data.frame(data.table::rbindlist(interval_positions)) %>%
    # CALCULATE INTERVAL SIZE
    dplyr::mutate(interval_size = endPOS - startPOS)
  
  # JOIN INTERVAL POSITIONS TO DATA FRAME CONTAINING CORRELATION INFORMATION AND PHENOTYPE INFORMATION
  Final_Processed_Mappings <- dplyr::left_join( correlation_df, interval_pos_df, 
                                                by = c("trait", "CHROM", "POS"),
                                                copy = TRUE )
  
  return(Final_Processed_Mappings)
  
}

# variant_correlation <- function(df, 
#                                 quantile_cutoff_high = .9, 
#                                 quantile_cutoff_low = .1,
#                                 genomicTrait = F){
#   
#   
#   if(genomicTrait == T){
#     # loosely identify unique QTL identified in mappings.
#     intervals <- df %>%
#       na.omit() %>%
#       dplyr::distinct(CHROM, startPOS, endPOS ) %>%
#       dplyr::distinct(CHROM, startPOS ) %>%
#       dplyr::distinct(CHROM, endPOS ) %>%
#       dplyr::arrange(CHROM, startPOS) 
#   } else {
#     # loosely identify unique QTL identified in mappings.
#     intervals <- df %>%
#       na.omit() %>%
#       dplyr::distinct(condition, CHROM, startPOS, endPOS ) %>%
#       dplyr::distinct(CHROM, startPOS ) %>%
#       dplyr::distinct(CHROM, endPOS ) %>%
#       dplyr::arrange(CHROM, startPOS) 
#   }
#   
#   
#   # unique strains to filter snpeff output for GWAS data - doesnt matter for genomic traits
#   strains <- as.character(na.omit(unique(df$strain)))
#   # set up the database to search for gene annotations using the biomart package
#   ensembl = biomaRt::useMart("ensembl",dataset="celegans_gene_ensembl")
#   
#   # initialize a list to store gene annotations for genes most highly correlated with phenotype
#   intervalGENES <- list()
#   
#   # loop through all unique intervals
#   for( i in 1:nrow(intervals)){
#     
#     print(paste(100*signif(i/nrow(intervals),3), "%",sep=""))
#     
#     
#     nstrains <- df %>%
#       dplyr::filter( trait == intervals[i,]$trait ) %>%
#       na.omit()
#     
#     nstrains <- length(unique(nstrains$strain))
#     
#     # define chromosome and left and right bound for intervals
#     CHROM <- as.character(intervals[i,]$CHROM)
#     left <- intervals[i,]$startPOS
#     right <- intervals[i,]$endPOS
#     
#     # define region of interest for Dan's snpeff function input
#     region_of_interest <- paste0(CHROM,":",left,"-",right)
#     
#     # run variant effect prediction function
#     snpeff_output <- snpeff(region = region_of_interest, impute = F) 
#     
#     # prune snpeff outputs
#     pruned_snpeff_output <- snpeff_output %>%
#       dplyr::filter( strain %in% strains ) %>% # only keep strains used in mappings
#       dplyr::filter( !is.na(impact) ) %>% # remove rows with NA in impact (looked to be mostly splice variants)
#       dplyr::distinct( CHROM, POS, strain, effect, gene_id ) %>% # remove duplicates
#       dplyr::arrange( effect ) %>% 
#       # pull out columns of interest
#       dplyr::select( CHROM, POS, REF, ALT, GT, effect, 
#                      nt_change, aa_change, gene_name, 
#                      gene_id, feature_type, strain) %>%
#       dplyr::group_by( CHROM, POS, effect) %>% # group for individual genes and effects
#       # make numeric allele column, this makes HETs and NAs in the GT column NAs - these are excluded from the correlation analysis
#       # need to elimante hets and NAs from GT
#       dplyr::filter(!is.na(GT), GT != "HET")%>%
#       # make numeric
#       dplyr::mutate(num_allele = ifelse(GT == "REF", 0, 
#                                         ifelse(GT == "ALT", 1, NA)))%>%
#       # determine if any alleles are present in less that 5% of the population
#       dplyr::mutate(num_alt_allele = sum(num_allele, na.rm=T),
#                     num_strains = n())%>%
#       # if they are, eliminate
#       dplyr::filter(num_alt_allele / num_strains > .05) %>%
#       dplyr::filter(num_strains > nstrains*.8)
#     
#     if( nrow(pruned_snpeff_output) > 0 ){
#       # pull unique interval from processed mapping DF to recover, phenotypes, strains, log10p, phenotype value
#       # this is useful to pull out all intervals with the same confidence interval that were pruned above.
#       interval_df <- df %>%
#         dplyr::filter( CHROM == CHROM, startPOS == left, endPOS = right )%>% # filter for confidence interval of interest
#         dplyr::group_by( trait, CHROM, startPOS,endPOS ) %>% # group by unique phenotype and interval
#         dplyr::filter( log10p == max(log10p) ) %>% # pull out most significant snp to minimize redundancy
#         dplyr::distinct( trait, startPOS, endPOS, peakPOS, strain) %>% 
#         dplyr::select( trait, startPOS, endPOS, peakPOS, strain, log10p, CHROM, pheno_value = value)
#       
#       # calculate the correlation between interval variants and the phenotype 
#       # pull out only the most correlated genes 
#       
#       pheno_snpeff_df <- pruned_snpeff_output %>%
#         dplyr::left_join(., interval_df, by = "strain", copy = TRUE) %>% # join snpeff variant df to phenotype df for a particular interval
#         dplyr::distinct( strain, trait, pheno_value, gene_name) %>% # remove redundancy
#         dplyr::group_by( trait, CHROM, POS, effect, feature_type) %>% # group_by unique variant and phenotype
#         dplyr::mutate(spearman_cor = cor(pheno_value, num_allele, method = "spearman", use = "pairwise.complete.obs"))%>% # calculate correlation
#         dplyr::ungroup()%>% # ungroup to calculate quantiles of correlations
#         # dplyr::mutate(q90 = quantile(spearman_cor, probs = .9, na.rm = T) )%>%
#         # we want to keep high positively correlated and high negatively correlated variants
#         dplyr::mutate(abs_spearman_cor = abs(spearman_cor))%>%
#         dplyr::filter(abs_spearman_cor > quantile(abs_spearman_cor, probs = quantile_cutoff_high, na.rm = T) )%>%
#         dplyr::ungroup()%>%
#         # organize DF by correlation
#         dplyr::arrange(desc(abs_spearman_cor))
#       
#       
#       # get gene annotations usining biomart package
#       # attributes are the columns you want to return
#       # filters are the columns you want to filter by, in this case we want to filter by wormbase_gene - e.g. WBGene00012953
#       # values are the values you want to be present in the filter column, i.e the genes you want information from
#       # this is pulled from the highly correlated variant DF above
#       # mart is defined above as the annoted c.elegans genome
#       gene_annotations <- getBM(attributes=c('entrezgene','go_id',"external_gene_name",
#                                              "external_transcript_name","gene_biotype",
#                                              "transcript_biotype","description", "family_description",
#                                              "name_1006","wormbase_gene"), 
#                                 filters = "wormbase_gene",
#                                 values = unique(pheno_snpeff_df$gene_id),
#                                 mart=ensembl) %>%
#         dplyr::distinct(entrezgene, go_id) %>% # pu;; distinct genes
#         dplyr::rename(gene_id = wormbase_gene) # change column for joining
#       
#       # attach the correlation coefficient to gene annotation data frame to minimize looking at multiple data frames
#       gene_cors <- pheno_snpeff_df %>%
#         dplyr::select(gene_id, spearman_cor)%>%
#         dplyr::distinct(gene_id, spearman_cor) %>%
#         dplyr::left_join(gene_annotations, ., by = "gene_id") %>%
#         dplyr::arrange(desc(spearman_cor))
#       
#       # append phenotype-snpeff-correlation DF and gene annotation DF to list for every unique interval
#       intervalGENES[[i]] <- list(pheno_snpeff_df, gene_cors)
#     } 
#     else
#     {
#       intervalGENES[[i]] <- list(NA, NA)
#     }
#     
#   }
#   
#   return(intervalGENES)
# }

# called functions - from Dan Cook

# snpeff <- function(region = "II:14524173..14525111",
#                    severity = c("HIGH","MODERATE"),
#                    long = TRUE,
#                    impute = TRUE) {
#   
#   if (impute == T) {
#     vcf_file = "20150731_WI_PASS.impute.snpeff.vcf.gz"
#   } else {
#     vcf_file = "20150731_WI_PASS.snpeff.vcf.gz"
#   }
#   
#   
#   
#   if (!grepl("(I|II|III|IV|V|X|MtDNA).*", region)) {
#     gene_ids <- read_tsv("~/Dropbox/Andersenlab/WormReagents/Variation/Andersen_VCF/wb_gene.txt")
#     wb_id <- filter(gene_ids, name == region)$ID
#     wb_url <- paste0("http://api.wormbase.org/rest/field/gene/",wb_id, "/location/")
#     wb_ret <- GET(wb_url, add_headers("Content-Type"="application/json"))
#     region <- content(wb_ret)$location$genomic_position$data[[1]]$pos_string
#   } 
#   
#   # Fix region to allow wb type spec.
#   region <- gsub("\\.\\.", "-", region)
#   
#   script_dir <- list.dirs("~/Dropbox/Andersenlab/WormReagents/Variation/Andersen_VCF/")[[1]]
#   command <- paste("python",
#                    "~/Dropbox/Andersenlab/WormReagents/Variation/Andersen_VCF/query.py", 
#                    region,
#                    paste0(script_dir,vcf_file),
#                    paste(severity, collapse=","))
#   
#   tsv <- read_tsv( pipe(command), na = "None") 
#   if (long == FALSE) {
#     tsv
#   } else {
#     tsv <- gather_(tsv, "strain", "GT", names(tsv)[21:length(tsv)])  %>%
#       tidyr::separate(GT, into=c("a1","a2"), sep="/|\\|", remove=T) %>%
#       dplyr::mutate(a1=ifelse(a1 == ".", NA, a1)) %>%
#       dplyr::mutate(a2=ifelse(a2 == ".", NA, a2)) %>%
#       dplyr::mutate(GT = NA) %>%
#       dplyr::mutate(GT = ifelse(a1 == REF & a2 == REF & !is.na(a1), "REF",GT)) %>%
#       dplyr::mutate(GT = ifelse(a1 != a2 & !is.na(a1), "HET",GT)) %>%
#       dplyr::mutate(GT = ifelse(a1 == a2 & a1 != REF & !is.na(a1), "ALT",GT)) %>%
#       dplyr::select(CHROM, POS, strain, REF, ALT, a1, a2, GT, everything()) %>%
#       dplyr::arrange(CHROM, POS) 
#     
#     tsv
#   }
# }

