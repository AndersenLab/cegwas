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
  
  if(length(unique(plot_df$trait)) == 1){
    plot_df %>%
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
                    y = expression(log[10](p)),
                    title = unique(plot_df$trait))
  } 
  else
  {
    plot_traits <- unique(plot_df$trait)
    for(i in length(plot_traits)){
      plot_df %>%
        dplyr::filter(trait == plot_traits[i]) %>%
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
                      y = expression(log[10](p)),
                      title = plot_traits[i])
    }
    
  }
}


#' PxG plot
#'
#' \code{pxg_plot} generates a boxplot of phenotypes split by genotype at QTL peak position using \code{ggplot2}
#'
#'
#' @param plot_df the output from the \code{gwas_mappings} function. 
#' @return Ouput is a ggplot object facetted by peak SNP (if there are multiple peaks in a mapping). 
#' Genotypes are encoded as 1 or 0 and are on the x-axis. Phenotypes are on y-axis.
#' @export

pxg_plot <- function(plot_df){
  
  if(length(unique(plot_df$trait)) == 1){
  plot_df %>%
    na.omit()%>%
    dplyr::distinct(strain, value, peakPOS)%>%
    dplyr::select(strain, value, CHROM, peakPOS, allele)%>%
    dplyr::mutate(chr_pos = paste(CHROM, peakPOS, sep="_"))%>%
    ggplot2::ggplot(.)+
    ggplot2::aes(x = factor(allele), y = value)+
    ggplot2::scale_fill_brewer(palette = "Set1")+
    ggplot2::geom_boxplot( ggplot2::aes(fill = factor(allele)))+
    ggplot2::theme_bw()+
    ggplot2::geom_jitter(alpha = .7)+
    ggplot2::facet_grid(.~chr_pos, scales = "free")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=24, face="bold", color="black"),
                   axis.text.y = ggplot2::element_text(size=24, face="bold", color="black"),
                   axis.title.x = ggplot2::element_text(size=24, face="bold", color="black", vjust=-.3),
                   axis.title.y = ggplot2::element_text(size=24, face="bold", color="black"),
                   strip.text.x = ggplot2::element_text(size=24, face="bold", color="black"),
                   strip.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
                   plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
                   legend.position="none",
                   panel.background = ggplot2::element_rect( color="black",size=1.2),
                   strip.background = ggplot2::element_rect(color = "black", size = 1.2))+
    ggplot2::labs(y = "Phenotype", x = "Genotype", title = unique(plot_df$trait))
  }
  else
  {
    plot_traits <- unique(plot_df$trait)
    for(i in length(plot_traits)){
      plot_df %>%
        na.omit()%>%
        dplyr::filter(trait == plot_traits[i]) %>%
        dplyr::distinct(strain, value, peakPOS) %>%
        dplyr::select(strain, value, CHROM, peakPOS, allele) %>%
        dplyr::mutate(chr_pos = paste(CHROM, peakPOS, sep="_")) %>%
        ggplot2::ggplot(.) + 
        ggplot2::aes(x = factor(allele), y = value) +
        ggplot2::scale_fill_brewer(palette = "Set1") +
        ggplot2::geom_boxplot( ggplot2::aes(fill = factor(allele))) +
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
        ggplot2::labs(y = "Phenotype", x = "Genotype", title = plot_traits[i])
    }
  }
}


