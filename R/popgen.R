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

#' Geographic distribution of alleles
#'
#' \code{allele_distribution} uses the ggplot2 package to plot the global distribution of an allele of interest.
#'
#'
#' @param locus A gene or position. Examples: \code{top-2}, \code{II:14425175}, or "WBGene00010785"
#' @param variant character value corresponding amino acid change you are interested in plotting. Currently only supports visualization of SNVs that alter an amino acid residue.
#' @param location region to restrict plot to. Uses google api for lookup. Examples: \code{europe}.
#' @param colors character vector containing two colors. The first element of the vector will be the ALT genotype and the second will be the REF genotype. Default values are c("purple","salmon").
#' @param alpha Set marker alpha.
#' @return Outputs a ggplot object containing a visualization of the global distribution of allele of interest
#' @importFrom dplyr %>%
#' @export

allele_distribution <- function(locus = "top-2", variant = 797, location = NA,  colors = c("purple", "salmon"), alpha = 0.8){
  
  query <- snpeff(locus, severity = "ALL", elements = "ALL") 
  
  if (!grepl(":", locus)) {
    allele_info <- query %>%
      dplyr::select(CHROM, POS, isotype = strain, aa_change,GT) %>%
      dplyr::filter(grepl(paste0("[^0-9]", variant, "[^0-9]"),aa_change)) 
  } else {
    chrom_pos <- stringr::str_split(locus, ":")[[1]]
    allele_info <- snpeff(locus) %>%
      dplyr::select(CHROM, POS, isotype = strain, aa_change,GT) %>%
      dplyr::filter(CHROM == chrom_pos[1], POS == chrom_pos[2]) 
  }
  allele_info <- allele_info %>% dplyr::distinct(isotype) %>%
    dplyr::left_join(., get_strain_isotype(), by = "isotype")
  
  rplot <- ggplot(allele_info) + 
    geom_path(data = map_data("world"), aes(x=long, y=lat, group=group)) +
    geom_point(data = subset(allele_info, GT == "REF"), 
               aes(x = longitude, y = latitude, color = colors[1], ordered = TRUE), 
               size = 3, position = "jitter", alpha = alpha) +
    geom_point(data = subset(allele_info, GT == "ALT"), 
               aes(x = longitude, y = latitude, color = colors[2], ordered = TRUE), 
               size = 3, position = "jitter", alpha = alpha)+ 
    scale_color_manual(values = colors, 
                       labels = c("REF", "ALT"),
                       name = "Genotype") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.title = element_text(face = "bold"))
  if(!is.na(location)) {
    bound <- ggmap::geocode(location, "more")
    rplot <- rplot + scale_x_continuous(limits = c(bound$west, bound$east)) +
      scale_y_continuous(limits = c(bound$south, bound$north))
  }
  rplot
}
