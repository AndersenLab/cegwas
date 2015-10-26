#' Manhattan Plot for GWAS Mapping Data
#'
#' \code{manplot} generates a manhattan plot using \code{ggplot2}
#'
#'
#' @param plot_df the output from the \code{gwas_mappings} function. 
#' @return Ouput is a ggplot object facetted by chromosome. SNPs above bonferroni corrected p-value (gray line) are colored blue.
#' Confidence interval for a given peak is highlighted in red.
#' @export


manplot <- function(plot_df) {
    plot_traits <- unique(plot_df$trait)
    plots <- lapply(plot_traits, function(i) {
        plot_df %>%
        dplyr::filter(trait == i) %>%
        dplyr::distinct(marker) %>%
        dplyr::mutate(n_snps = n())%>%
        ggplot2::ggplot(.) +
        ggplot2::aes(x = POS/1e6, y = log10p) +
        ggplot2::scale_color_manual(values = c("black","blue","red")) +
        ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6, 
                               xmax = endPOS/1e6, 
                               ymin = 0, 
                               ymax = Inf, 
                               fill = "blue", 
                               alpha=.1), 
                           color = NA)+
        ggplot2::geom_hline(ggplot2::aes(yintercept = BF),
                            color = "gray", 
                            alpha = .75,  
                            size = 1) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(.05/n_snps)),
                            color = "gray", 
                            alpha = .75,  
                            size = 1)+
        ggplot2::geom_point( ggplot2::aes(color= factor(aboveBF)) ) +
        ggplot2::facet_grid( . ~ CHROM, scales = "free_x" ) +
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
                      y = expression(-log[10](p)),
                      title = i)
    })
    plots
  }



#' PxG plot
#'
#' \code{pxg_plot} generates a boxplot of phenotypes split by genotype at QTL peak position using \code{ggplot2}
#'
#'
#' @param plot_df the output from the \code{gwas_mappings} function. 
#' @param loc custom location to output the pxg plot. Specified as CHROM:POS (e.g. II:1023423)
#' @return Ouput is a ggplot object facetted by peak SNP (if there are multiple peaks in a mapping). 
#' Genotypes are encoded as REF or ALT and are on the x-axis. Phenotypes are on y-axis.
#' @export

pxg_plot <- function(plot_df, loc = NA){
    plot_traits <- unique(plot_df$trait)
    plots <- lapply(plot_traits, function(x) {
      plot_peak <- plot_df %>%
        na.omit() %>%
        dplyr::filter(trait == x) %>%
        dplyr::distinct(strain, value, peakPOS) %>%
        dplyr::select(strain, value, CHROM, peakPOS, -allele) %>%
        dplyr::mutate(chr_pos = paste(CHROM, peakPOS, sep="_"))
      
      if (is.na(loc)) {
      loc <- plot_peak %>% dplyr::select(CHROM, peakPOS) %>% 
                    dplyr::distinct() %>% 
                    dplyr::transmute(CHROM_POS = paste0(CHROM, ":",peakPOS))
      }
      
      snpeff(loc[[1]], severity = c("MODIFIER", "LOW", "MODERATE", "HIGH")) %>%
      dplyr::select(CHROM, POS, strain, GT) %>%
      dplyr::distinct() %>%
      dplyr::left_join(plot_peak) %>%
      dplyr::filter(!is.na(GT)) %>%
      ggplot2::ggplot(., aes(x = GT, y = value, fill = GT)) + 
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::geom_boxplot() +
        ggplot2::theme_bw() +
        ggplot2::geom_jitter(alpha = .7) +
        ggplot2::facet_grid(.~chr_pos, scales = "free") +
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
        ggplot2::labs(y = "Phenotype", x = "Genotype", title = x)
    })
  plots
}


#' Plot variants for gene
#'
#' \code{gene_variants} generates a plot to visualize presence of variants for a particular gene of interest
#'
#'
#' @param gene is the gene of interest in gene name (e.g. "top-2", "pot-2") or wormbase gene ID (e.g. "WBGene00010785") format. 
#' @return Ouput is a list of ggplot objects with strains on the Y axis and variants for gene of interest on X axis. Tiles are colored by variant or reference call. 
#' @examples  gene_variants(gene = c("top-2","pot-2"))
#' @examples  test[[1]]
#' @examples test[[2]]
#' @export

gene_variants <- function(gene){
  
  gene_variant_plot <- list()
  
  for(i in 1:length(gene)){
    
    gene_info <- snpeff(gene[i])
    
    gene_to_plot <- gene_info %>%
      dplyr::arrange(POS, GT) %>%
      dplyr::mutate(fac_aa = factor(aa_change, 
                                    ordered = is.ordered(aa_change),
                                    levels = aa_change,
                                    labels = aa_change)) %>%
      dplyr::filter(!is.na(GT), GT != "HET") %>%
      tidyr::complete(CHROM, POS, strain) %>%
      tidyr::fill(fac_aa) %>%
      dplyr::mutate(GT = ifelse(is.na(GT), "ZMISSING", GT))
  
    
    gene_variant_plot[[i]] <- ggplot2::ggplot(gene_to_plot) +
      ggplot2::aes(x = fac_aa, 
                   y = factor(strain, 
                              ordered = is.ordered(GT),
                              levels = strain,
                              labels = strain), 
                   fill = GT) +
      ggplot2::scale_fill_manual(values = c("#F97848","#FFFFB2","white"), breaks=c("REF", "ALT")) +
      ggplot2::geom_tile(color = "black") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=12, face="bold", color="black", angle = 60, vjust= 1.2, hjust = 1.2),
                     axis.text.y = ggplot2::element_text(size=12, face="bold", color="black"),
                     axis.title.x = ggplot2::element_text(size=18, face="bold", color="black", vjust=-.3),
                     axis.title.y = ggplot2::element_text(size=18, face="bold", color="black"),
                     strip.text.x = ggplot2::element_text(size=14, face="bold", color="black"),
                     strip.text.y = ggplot2::element_text(size=14, face="bold", color="black"),
                     plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
                     panel.background = ggplot2::element_blank(),
                     strip.background = ggplot2::element_rect(color = "black", size = 1.2)) +
      ggplot2::labs(y = "Strain", x = "Variant") +
      ggplot2::scale_x_discrete(expand = c(0,0))
  }
  
  return(gene_variant_plot)
}
