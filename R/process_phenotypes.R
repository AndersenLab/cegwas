#' Process Phenotype Data
#'
#' \code{process_pheno} Processes raw data file for GWAS mapping using \code{\link{gwas_mappings}} function
#'
#' This function takes raw phenotype data and eliminates outlier strains with a modified version of \code{\link{bamf_prune}} from the easysorter package.
#' Additionally it eliminates any traits that have the same values for >95% of the strains (important for binary traits).
#' If multiple strains are present that fall into the same isotype, they can be removed by setting \code{\link{remove_strains}} to 
#'
#' @param data is a dataframe containing phenotype data. The dataframe can be structured in wide or long format. %
#' \cr
#' \cr
#' \bold{long format} - The first columns should be \code{strain} whereas additional columns are traits. 
#' One row list a strain followed by all of that strains phenotypes.
#' \cr\cr
#' \code{strain}, \code{trait1}, \code{trait2}, \code{...}
#' \cr\cr
#' \bold{wide format} - The first column should be named \code{trait} and subsequent 
#' all additional columns should be strains. One row corresponding to one trait for all strains.
#' \cr\cr
#' \code{trait}, \code{strain1}, \code{strain2}, \code{...}
#' \cr\cr
#' @param remove_strains Remove strains with no known isotype. Default is TRUE.
#' @param duplicate_method Method for dealing with the presence of multiple strains falling into the same isotype. Either \code{"average"} or \code{"first"}.
#' @return Outputs a list. The first element of the list is an ordered vector of traits. 
#' The second element of the list is a dataframe containing one column for each strain, with values corresponding to traits in element 1 for rows.
#' @importFrom dplyr %>%
#' @export

process_pheno <- function(data, remove_strains = TRUE, duplicate_method = "first"){
  
  # Rename first column strain
  colnames(data)[1] <- "strain"
  
  # Reshape data from wide to long
  if (sum(row.names(cegwas::kinship) %in% data[[1]]) > 3) {
    names(data) <- stringr::str_to_lower(names(data))
    data <- data %>% tidyr::gather(trait,value,-strain) %>% tidyr::spread(strain, value)  
  }
  
  # Strain - Isotype Issues; Duplicate check
  data <- tidyr::gather(data, "strain", "val", 2:ncol(data)) %>%
          dplyr::mutate(strain = as.character(strain)) %>%
          dplyr::mutate(isotype = strain_isotype[strain]) %>%
          dplyr::mutate(warnings = strain_warnings[strain]) %>%
          dplyr::group_by(trait, isotype) %>% 
          dplyr::mutate(iso_count = n()) 

  # Warn user of potential issues with strains
  issue_warnings <- dplyr::filter(dplyr::ungroup(data), !is.na(warnings)) %>% 
                    dplyr::select(strain, warnings) %>%
                    dplyr::group_by(strain, warnings) %>%
                    unique()
  if (nrow(issue_warnings) > 0) {
    for(x in 1:nrow(issue_warnings)) {
        warn <- issue_warnings[x,]
        warn <- paste(warn$strain, ":", warn$warnings)
        warning(warn, call. = F)
    }
  }
  
  # See if any strains with no known isotypes are used and drop.
  if (sum(is.na(data$isotype)) > 0) { 
    if (remove_strains == T) {
      warning("Missing isotypes for the following strains: ",
              paste(unique(data$strain[is.na(data$isotype)]), collapse = ", "),
              "\nNo known isotype! Removing strains.")
      data <- dplyr::filter(data, !is.na(isotype))
    } else {
      stop("Missing isotypes for the following strains: ",
           paste(unique(data$strain[is.na(data$isotype)]), collapse = ", "),
           "\nNo known isotype! Please remove strain(s) or set remove_strains = TRUE.")
    }
  }

  # Handle duplicates
  repeat_isotypes <- dplyr::filter(dplyr::ungroup(data), iso_count > 1) %>% 
  dplyr::select(isotype, strain) %>%
  dplyr::distinct()
  if (nrow(repeat_isotypes) > 0) {
    warning("Strains were phenotyped that belong to the same isotype:", immediate. = T)
    write.table(repeat_isotypes, "", quote = F, row.names = F, sep = "\t")
    if (duplicate_method == "first") {
      data <- dplyr::group_by(data, trait, isotype) %>% 
        dplyr::mutate(first = row_number()) %>%
        dplyr::filter(first == 1) %>%
        dplyr::select(-first)
      message("Using the first strain in each isotype group.")
    } else if (duplicate_method == "average") {
      data <- dplyr::group_by(data, trait, isotype) %>% 
        dplyr::mutate(val = mean(val)) %>%
        dplyr::distinct()
      message("Taking average of strains belonging to the same isotype group.")
    }
  }
  
  
  
  # Return data frame to previous state
  data <- dplyr::select(data, -strain, -warnings, -iso_count) %>%
  dplyr::rename(strain = isotype) %>%
  tidyr::spread(strain, val)
  
  # identify any traits that only have 1 unique value
  pheno <- data.frame(data.frame(data)[,1:ncol(data)], 
                      uniq = data.frame(uniq = c(apply(data.frame(data)[,2:ncol(data)], 1, function(x) length(unique(x))))),
                      row.names = data$trait)
  
  # remove identified traits
  phen <- pheno %>%
    dplyr::filter(uniq != 1)%>%
    dplyr::select(-uniq)
  
  # removes binary phenotypes where one phenotype is in less than 5% of strains
  phen1 <- remove_lowFreq_phenotypes(phen)
  # run modified version of bamf_prune from easy sorter package
  # does the same thing as bamf_prune, just different input data structure
  phen2 <- mod_bamf_prune(phen1)
  # removes binary phenotypes where one phenotype is in less than 5% of strains
  phen3 <- remove_lowFreq_phenotypes(phen2, wide = FALSE)
  
  # make into long formated data
  phen4 <- phen3 %>%
    tidyr::spread(strain,value)
  
  # generate traits list
  traits <- phen4$trait
  
  # remove trait column from phenotype data frame
  phen5 <- data.frame(phen4) %>% 
    dplyr::select(-trait)
  
  # make phenotypes numeric
  phen6 <- data.frame(lapply(phen5, function(x) as.numeric(as.character(x))))
  
  # generate list data structure to feed into mapping function
  output <- list(traits, phen6)
  
  return(output)
}

# called functions
# function to remove phenotypes with low frequencies (only important for binary traits)
remove_lowFreq_phenotypes <- function(data, wide = TRUE){
  
  if(wide == T){
    output <- data.frame(data)%>%
      tidyr::gather(strain, value, -trait)%>% # make long
      dplyr::group_by(trait)%>% # group
      dplyr::mutate(nst = n(), # number strain per trait
             nuniq = length(unique(value)))%>% # number of unique values
      dplyr::mutate(u1 = ifelse(nuniq == 2, unique(value)[1], -9999), # if two unique values, add first unique value
             u2 = ifelse(nuniq == 2, unique(value)[2], -9999))%>% # if two unique values, add second unique value
      dplyr::mutate(count1 = ifelse(u1 != -9999 & value == unique(value)[1], 1, 0), # add 1 if value = first unique (requires only 2 uniques)
             count2 = ifelse(u2 != -9999 & value == unique(value)[2], 1, 0))%>% # add 1 if value = second unique (requires only 2 uniques)
      dplyr::mutate(freq1 = sum(count1)/nst, # get frequency of values from traits where there are only two unique trait values
             freq2 = sum(count2)/nst)%>% # get frequency of values from traits where there are only two unique trait values
      dplyr::mutate(cuts = ifelse(u1 != -9999 & 
                             u2 != -9999 & 
                             (freq1 < .05 | freq1 > .95), 1,0))%>% # add ID for traits that have less than 5% of strains with one trait
      dplyr::filter(cuts != 1)%>% # remove traits identified by cut1
      dplyr::select(trait, strain, value) # remove columns used for cut
  } 
  else
  {
    output <- data.frame(data)%>%
      dplyr::group_by(trait)%>% # group
      dplyr::mutate(nst = n(), # number strain per trait
             nuniq = length(unique(value)))%>% # number of unique values
      dplyr::mutate(u1 = ifelse(nuniq == 2, unique(value)[1], -9999), # if two unique values, add first unique value
             u2 = ifelse(nuniq == 2, unique(value)[2], -9999))%>% # if two unique values, add second unique value
      dplyr::mutate(count1 = ifelse(u1 != -9999 & value == unique(value)[1], 1, 0), # add 1 if value = first unique (requires only 2 uniques)
             count2 = ifelse(u2 != -9999 & value == unique(value)[2], 1, 0))%>% # add 1 if value = second unique (requires only 2 uniques)
      dplyr::mutate(freq1 = sum(count1)/nst, # get frequency of values from traits where there are only two unique trait values
             freq2 = sum(count2)/nst)%>% # get frequency of values from traits where there are only two unique trait values
      dplyr::mutate(cuts = ifelse(u1 != -9999 & 
                             u2 != -9999 & 
                             (freq1 < .05 | freq1 > .95), 1,0))%>% # add ID for traits that have less than 5% of strains with one trait
      dplyr::filter(cuts != 1)%>% # remove traits identified by cut1
      dplyr::select(trait, strain, value) # remove columns used for cut
  }
  
  return(output)
}

# modified version of bamf_prune
# works the same as from the easysorter package, but with different input format
mod_bamf_prune <- function(data){
  
  napheno <- data[is.na(data$value), ] %>%
    dplyr::mutate(bamfoutlier1 = NA, bamfoutlier2 = NA, bamfoutlier3 = NA)
  
  datawithoutliers <- data %>%
    
    # Filter out all of the wash and/or empty wells
    
    dplyr::filter(!is.na(strain)) %>%
    
    # Group on condition and trait, the, calculate the first and third
    # quartiles for each of the traits
    
    dplyr::group_by(trait) %>%
    dplyr::summarise(iqr = IQR(value, na.rm = TRUE),
                     q1 = quantile(value, probs = .25, na.rm = TRUE),
                     q3 = quantile(value, probs = .75,
                                   na.rm = TRUE)) %>%
    
    # Add a column for the boundaries of each of the bins
    
    dplyr::mutate(cut1h = q3 + (iqr * 2),
                  cut1l =q1 - (iqr * 2),
                  cut2h = q3 + (iqr * 3),
                  cut2l =q1 - (iqr * 3),
                  cut3h = q3 + (iqr * 4),
                  cut3l =q1 - (iqr * 4),
                  cut4h = q3 + (iqr * 5),
                  cut4l =q1 - (iqr * 5),
                  cut5l = q1 - (iqr * 7),
                  cut5h = q3 + (iqr * 7),
                  cut6l = q1 - (iqr * 10),
                  cut6h = q3 + (iqr * 10)) %>%
    
    # Join the bin boundaries back to the original data frame
    
    dplyr::left_join(data, ., by=c("trait")) %>%
    
    # Add columns tallying the total number of points in each of the bins
    
    dplyr::mutate(onehs = ifelse( cut2h > value & value >= cut1h,
                                  1, 0),
                  onels = ifelse( cut2l < value & value <= cut1l,
                                  1, 0),
                  twohs = ifelse( cut3h > value & value >= cut2h,
                                  1, 0),
                  twols = ifelse( cut3l < value & value <= cut2l,
                                  1, 0),
                  threehs = ifelse(cut4h > value & value >= cut3h,
                                   1, 0),
                  threels = ifelse(cut4l < value & value <= cut3l,
                                   1, 0),
                  fourhs = ifelse(cut5h > value &  value >= cut4h,
                                  1, 0),
                  fourls = ifelse(cut5l < value &  value <= cut4l,
                                  1, 0),
                  fivehs = ifelse(cut6h > value & value >= cut5h,
                                  1, 0),
                  fivels = ifelse(cut6l < value & value <= cut5l,
                                  1, 0),
                  sixhs = ifelse(value >= cut6h, 1, 0),
                  sixls = ifelse(value <= cut6l, 1, 0)) %>%
    
    # Group on condition and trait, then sum the total number of data points
    # in each of the IQR multiple bins
    
    dplyr::group_by(trait) %>%
    dplyr::mutate(s1h = sum(onehs, na.rm = TRUE),
                  s2h = sum(twohs, na.rm = TRUE),
                  s3h = sum(threehs, na.rm = TRUE),
                  s4h = sum(fourhs, na.rm = TRUE),
                  s5h = sum(fivehs, na.rm = TRUE),
                  s1l = sum(onels, na.rm = TRUE),
                  s2l = sum(twols, na.rm = TRUE),
                  s3l = sum(threels, na.rm = TRUE),
                  s4l = sum(fourls, na.rm = TRUE),
                  s5l = sum(fivels, na.rm = TRUE),
                  s6h = sum(sixhs, na.rm = TRUE),
                  s6l = sum(sixls, na.rm = TRUE))%>%
    
    # Group on condition and trait, then check to see if the number of
    # points in each bin is more than 5% of the total number of data points
    
    dplyr::group_by(trait) %>%
    dplyr::mutate(p1h = ifelse(sum(onehs, na.rm = TRUE) / n() >= .05,1,0),
                  p2h = ifelse(sum(twohs, na.rm = TRUE) / n() >= .05,1,0),
                  p3h = ifelse(sum(threehs, na.rm = TRUE) / n() >= .05,1,0),
                  p4h = ifelse(sum(fourhs, na.rm = TRUE) / n() >= .05,1,0),
                  p5h = ifelse(sum(fivehs, na.rm = TRUE) / n() >= .05,1,0),
                  p6h = ifelse(sum(sixhs, na.rm = TRUE) / n() >= .05,1,0),
                  p1l = ifelse(sum(onels, na.rm = TRUE) / n() >= .05,1,0),
                  p2l = ifelse(sum(twols, na.rm = TRUE) / n() >= .05,1,0),
                  p3l = ifelse(sum(threels, na.rm = TRUE) / n() >= .05,1,0),
                  p4l = ifelse(sum(fourls, na.rm = TRUE) / n() >= .05,1,0),
                  p5l = ifelse(sum(fivels, na.rm = TRUE) / n() >= .05,1,0),
                  p6l = ifelse(sum(sixls,
                                   na.rm = TRUE) / n() >= .05,1,0)) %>%
    
    # Count the number of observations in each condition/trait combination
    
    dplyr::mutate(numst = n()) %>%
    
    # Group on condition and trait, then filter out NAs in any of the added
    # columns
    
    dplyr::group_by(trait) %>%
    dplyr::filter(!is.na(trait), !is.na(value), !is.na(iqr), !is.na(q1),
                  !is.na(q3), !is.na(cut1h), !is.na(cut1l), !is.na(cut2h),
                  !is.na(cut2l), !is.na(cut3h), !is.na(cut3l),
                  !is.na(cut4h), !is.na(cut4l), !is.na(cut5l),
                  !is.na(cut5h), !is.na(cut6l), !is.na(cut6h),
                  !is.na(onehs), !is.na(onels), !is.na(twohs),
                  !is.na(twols), !is.na(threehs), !is.na(threels),
                  !is.na(fourhs), !is.na(fourls), !is.na(fivehs),
                  !is.na(fivels), !is.na(sixhs), !is.na(sixls),
                  !is.na(s1h), !is.na(s2h), !is.na(s3h), !is.na(s4h),
                  !is.na(s5h), !is.na(s1l), !is.na(s2l), !is.na(s3l),
                  !is.na(s4l), !is.na(s5l), !is.na(s6h), !is.na(s6l),
                  !is.na(p1h), !is.na(p2h), !is.na(p3h), !is.na(p4h),
                  !is.na(p5h), !is.na(p6h), !is.na(p1l), !is.na(p2l),
                  !is.na(p3l), !is.na(p4l), !is.na(p5l), !is.na(p6l),
                  !is.na(numst)) %>%
    
    # Add three columns stating whether the observation is an outlier
    # based the three outlier detection functions below
    
    dplyr::ungroup() %>%
    
    dplyr::mutate(cuts = categorize1(.),
                  cuts1 = categorize2(.),
                  cuts2 = categorize3(.))
  
  
  output <- datawithoutliers %>%
    dplyr::rename(bamfoutlier1 = cuts, 
                  bamfoutlier2 = cuts1,
                  bamfoutlier3 = cuts2)%>%
    dplyr::filter(!bamfoutlier1 & !bamfoutlier2 & !bamfoutlier3)%>%
    dplyr::select(trait, strain, value)
  
  return(output)
}

categorize1 <- function(data) {
  with(data,
       (sixhs >= 1 & ( (s6h + s5h + s4h ) / numst) <= .05
        & (s5h == 0 | s4h == 0))
       | (sixls >= 1 & ( (s6l + s5l + s4l) / numst) <= .05
          & (s5l == 0 | s4l == 0))
  )
}

# If the 5 innermost bins are discontinuous by more than a 1 bin gap, the
# observation is in the fifth bin (between 7 and 10x IQR outside the
# distribution), and the four outermost bins make up less than 5% of the
# population, mark the observation an outlier

categorize2 <- function(data) {
  with(data,
       ( ( ! ( (s5h >= 1 & s4h >= 1 & s3h >= 1 & s2h >= 1 & s1h >= 1)
               | (s5h >= 1 & s3h >= 1 & s2h >= 1 & s1h >= 1)
               | (s5h >= 1 & s4h >= 1 & s2h >= 1 & s1h >= 1)
               | (s5h >= 1 & s4h >= 1 & s3h >= 1 & s1h >= 1)
               | (s5h >= 1 & s4h >= 1 & s3h >= 1 & s2h >= 1)))
         & (fivehs == 1 & ( (s6h + s5h + s4h + s3h) / numst) <= .05))
       | ( ( ! ( (s5h >= 1 & s4l >= 1 & s3l >= 1 & s2l >= 1 & s1l >= 1)
                 | (s5h >= 1 & s3l >= 1 & s2l >= 1 & s1l >= 1)
                 | (s5h >= 1 & s4l >= 1 & s2l >= 1 & s1l >= 1)
                 | (s5h >= 1 & s4l >= 1 & s3l >= 1 & s1l >= 1)
                 | (s5h >= 1 & s4l >= 1 & s3l >= 1 & s2l >= 1)))
           & (fivels == 1 & ( (s6l + s5l + s4l + s3l) / numst) <= .05))
  )
}

# If the 4 innermost bins are discontinuous by more than a 1 bin gap, the
# observation is in the fourth bin (between 5 and 7x IQR outside the
# distribution), and the four outermost bins make up less than 5% of the
# population, mark the observation an outlier

categorize3 <- function(data) {
  with(data,
       ( ( ! ( (s4h >= 1 & s3h >= 1 & s2h >= 1 & s1h >= 1)
               | (s4h >= 1 & s2h >= 1 & s1h >= 1)
               | (s4h >= 1 & s3h >= 1 & s1h >= 1)
               | (s4h >= 1 & s3h >= 1 & s2h >= 1)))
         & (fourhs == 1 & (s5h + s4h + s3h + s2h) / numst <= .05))
       | ( ( ! ( (s4l >= 1 & s3l >= 1 & s2l >= 1 & s1l >= 1)
                 | (s4l >= 1 & s2l >= 1 & s1l >= 1)
                 | (s4l >= 1 & s3l >= 1 & s1l >= 1)
                 | (s4l >= 1 & s3l >= 1 & s2l >= 1)))
           & (fourls == 1 & (s5l + s4l + s3l + s2l) / numst <= .05))
  )
}