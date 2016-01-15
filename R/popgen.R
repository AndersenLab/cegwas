#' Tajima's D
#'
#' \code{tajimas_d} uses the PopGenome package to calculate Tajima's D across an interval. 
#'
#' This is the detail section if you want to fill out in the future
#'
#' @param vcf_path character value of the directory path where tabix-indexed VCF file is
#' @param vcf_name character value corresponding to tabix-indexed VCF file name
#' @param chromosome character value corresponding to chromosome of interest. This should be the same chromosome name that is present in the VCF file
#' @param interval_start numeric value start of interval of interest
#' @param interval_end numeric value end of interval of interest
#' @param window_size size of window to calculate Tajima's D in SNPs (default = 300)
#' @param slide_distance number of SNPs to shift window (default = 100)
#' @param samples character vector of strains to use for analysis (default is 152 strain list present in cegwas package)
#' @param outgroup character value corresponding to strain (or strains) to serve as reference (default value is "N2")
#' @param site_of_interest numeric value, if specified adds a line to the plot at this position. 
#' @return Outputs a list, where the first element is a dataframe containing Tajima's D calculations and the second is a ggplot2 object visualizing Tajima's D across the interval
#' @importFrom dplyr %>%
#' @export

tajimas_d <- function(vcf_path = paste0(path.package("cegwas"),"/"),
                      vcf_name = "WI.20160106.impute.vcf.gz", 
                      chromosome = "II", 
                      interval_start = 11021073, 
                      interval_end = 12008179, 
                      window_size = 300, 
                      slide_distance = 100, 
                      samples = colnames(snps[,3:ncol(snps)]),
                      outgroup = "N2",
                      site_of_interest = 11875145){
  
  setwd(vcf_path)
                        
  gen <- PopGenome::readVCF(vcf_name, 
                            numcols = 10000, 
                            tid = chromosome, 
                            frompos = interval_start, 
                            topos = interval_end, 
                            samplenames = samples, 
                            approx = TRUE) 
  
  gen1 <- PopGenome::set.populations(gen, list(samples), diploid = FALSE)
  
  gen2 <- PopGenome::set.outgroup(gen1,outgroup, diploid = FALSE)
  
  s_gen <- PopGenome::sliding.window.transform(gen2, width = window_size, jump = slide_distance, whole.data = FALSE)
  
  test <- data.frame(snps = 1:length(as.numeric(colnames(as.data.frame(s_gen@BIG.BIAL[[1]][,])))), position = as.numeric(colnames(as.data.frame(s_gen@BIG.BIAL[[1]][,]))))
  
  slides <- list()
  for(i in 1:length(s_gen@SLIDE.POS)){
    slides[[i]] <- data.frame(window = i, snps = s_gen@SLIDE.POS[[i]])
  }
  
  slides <- dplyr::rbind_all(slides) %>%
    dplyr::left_join(data.frame(test), ., by = "snps")
  
  
  s_gen_stats <- PopGenome::neutrality.stats(s_gen)
  
  td <- data.frame(window = 1:length(s_gen@SLIDE.POS), Td = s_gen_stats@Tajima.D) %>%
    dplyr::left_join(slides, ., by = "window")%>%
    dplyr::group_by(window)%>%
    dplyr::mutate(swindow = min(as.numeric(position)),
                  ewindow= max(as.numeric(position)))%>%
    dplyr::mutate(midwindow = (swindow+ewindow)/2)%>%
    dplyr::rename(Td = pop.1)%>%
    dplyr::distinct(Td, window)
  
  tajimas_d_plot <- ggplot2::ggplot(td)+
    ggplot2::aes(x = position/1e6, y = Td)+
    ggplot2::geom_point(size = 2)+
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=24, face="bold", color="black"),
                   axis.text.y = ggplot2::element_text(size=24, face="bold", color="black"),
                   axis.title.x = ggplot2::element_text(size=24, face="bold", color="black", vjust=-.3),
                   axis.title.y = ggplot2::element_text(size=24, face="bold", color="black"),
                   strip.text.x = ggplot2::element_text(size=24, face="bold", color="black"),
                   strip.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
                   plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
                   legend.position="none",
                   panel.background = ggplot2::element_rect( color="black",size=1.2),
                   strip.background = ggplot2::element_rect(color = "black", size = 1.2)) +
    ggplot2::labs(x = "Genomic Position (Mb)",
                  y = "Tajima's D")
  
  if(!is.null(site_of_interest)){
    tajimas_d_plot <- tajimas_d_plot + 
      ggplot2::geom_vline(ggplot2::aes(xintercept = site_of_interest/1e6), color = "red", alpha = .7, size = 2)
  }
  
  return(list(td, tajimas_d_plot))
}

#' Global distribution of allele
#'
#' \code{global_allele_distribution} uses the ggplot2 and ggmap packages to plot the global distribution of an allele of interest
#'
#' This is the detail section if you want to fill out in the future
#'
#' @param gene character value gene of interest. Can be gene name "top-2" or gene id "WBGene00010785"
#' @param variant character value corresponding amino acid change you are interested in plotting. Currently only supports visualization of SNVs that alter an amino acid residue.
#' @param colors character vector containing two colors. The first element of the vector will be the ALT genotype and the second will be the REF genotype. Default values are c("purple","salmon").
#' @return Outputs a ggplot object containing a visualization of the global distribution of allele of interest
#' @importFrom dplyr %>%
#' @export

global_allele_distribution <- function(gene = "top-2", variant = "797", colors = c("purple", "salmon")){
  ###################################################################################################
  # Recentre worldmap (and Mirrors coordinates) on longitude 160
  ### Code by Claudia Engel  March 19, 2012, www.stanford.edu/~cengel/blog
  
  ### Recenter ####
  center <- 160 # positive values only
  
  # shift coordinates to recenter CRAN Mirrors
  isotype_location$long.recenter <- ifelse(isotype_location$longitude < center - 180 , isotype_location$longitude + 360, isotype_location$longitude)
  isotype_location$lat <- isotype_location$latitude
  
  
  allele_info <- snpeff(gene) %>%
    dplyr::select(CHROM, POS, isotype = strain, aa_change,GT) %>%
    dplyr::filter(grepl(variant,aa_change))%>%
    dplyr::distinct(isotype)%>%
    dplyr::left_join(., isotype_location, by = "isotype")
  
  # shift coordinates to recenter worldmap
  worldmap <- ggplot2::map_data("world")
  worldmap$long.recenter <- ifelse(worldmap$long < center - 180 , worldmap$long + 360, worldmap$long)
  
  # now regroup
  worldmap.rg <- plyr::ddply(worldmap, .(group), RegroupElements, "long.recenter", "group")
  
  # close polys
  worldmap.cp <- plyr::ddply(worldmap.rg, .(group.regroup), ClosePolygons, "long.recenter", "order") # use the new grouping var
  #############################################################################
  
  # Plot worldmap using data from worldmap.cp
  
  worldmap = ggplot2::ggplot(aes(x = long.recenter, y = lat), data = worldmap.cp) + 
    ggplot2::geom_polygon(aes(group = group.regroup), fill="#f9f9f9", colour = "grey65") + 
    ggplot2::scale_y_continuous(limits = c(-60, 85)) + 
    ggplot2::coord_equal() +  
    ggplot2::theme_bw() + 
    ggplot2::theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(), 
      plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1),
      panel.border = element_rect(colour = "black"))
  
  # plot_title <- paste0("Global Distribution of ", gene," ", unique(allele_info$aa_change), " allele")
  
  # Plot the CRAN Mirrors
  worldmap = worldmap + ggplot2::geom_point(data = allele_info, aes(long.recenter, latitude, color = GT),
                                            pch = 19, size = 2, alpha = .7) +
    ggplot2::scale_color_manual(values = c(colors[1], colors[2]),
                                name = "Genotype")
  
  return(worldmap)
}

### Function to regroup split lines and polygons
# Takes dataframe, column with long and unique group variable, returns df with added column named group.regroup
RegroupElements <- function(df, longcol, idcol){
  g <- rep(1, length(df[,longcol]))
  if (diff(range(df[,longcol])) > 300) { # check if longitude within group differs more than 300 deg, ie if element was split
    d <- df[,longcol] > mean(range(df[,longcol])) # we use the mean to help us separate the extreme values
    g[!d] <- 1 # some marker for parts that stay in place (we cheat here a little, as we do not take into account concave polygons)
    g[d] <- 2 # parts that are moved
  }
  g <- paste(df[, idcol], g, sep=".") # attach to id to create unique group variable for the dataset
  df$group.regroup <- g
  df
}

### Function to close regrouped polygons
# Takes dataframe, checks if 1st and last longitude value are the same, if not, inserts first as last and reassigns order variable
ClosePolygons <- function(df, longcol, ordercol){
  if (df[1,longcol] != df[nrow(df),longcol]) {
    tmp <- df[1,]
    df <- rbind(df,tmp)
  }
  o <- c(1: nrow(df)) # rassign the order variable
  df[,ordercol] <- o
  df
}



