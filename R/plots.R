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
      ggplot2::geom_rect(aes(xmin = startPOS/1e6, 
                             xmax = endPOS/1e6, 
                             ymin = 0, 
                             ymax = Inf, 
                             fill = "blue", 
                             alpha=.1), 
                         color = NA)+
      ggplot2::geom_hline(aes(yintercept = BF),
                          color = "gray", 
                          alpha = .75,  
                          size = 1)+
      ggplot2::geom_point( aes(color= factor(aboveBF)) ) +
      ggplot2::facet_grid( . ~ CHROM, scales = "free_x" ) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = element_text(size=24, face="bold", color="black"),
                     axis.text.y = element_text(size=24, face="bold", color="black"),
                     axis.title.x = element_text(size=24, face="bold", color="black", vjust=-.3),
                     axis.title.y = element_text(size=24, face="bold", color="black"),
                     strip.text.x = element_text(size=24, face="bold", color="black"),
                     strip.text.y = element_text(size=16, face="bold", color="black"),
                     plot.title = element_text(size=24, face="bold", vjust = 1),
                     legend.position="none",
                     panel.background = element_rect( color="black",size=1.2),
                     strip.background = element_rect(color = "black", size = 1.2)) +
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
        ggplot2::geom_rect(aes(xmin = startPOS/1e6, 
                               xmax = endPOS/1e6, 
                               ymin = 0, 
                               ymax = Inf, 
                               fill = "blue", 
                               alpha=.1), 
                           color = NA)+
        ggplot2::geom_hline(aes(yintercept = BF),
                            color = "gray", 
                            alpha = .75,  
                            size = 1)+
        ggplot2::geom_point( aes(color= factor(aboveBF)) ) +
        ggplot2::facet_grid( . ~ CHROM, scales = "free_x" ) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = element_text(size=24, face="bold", color="black"),
                       axis.text.y = element_text(size=24, face="bold", color="black"),
                       axis.title.x = element_text(size=24, face="bold", color="black", vjust=-.3),
                       axis.title.y = element_text(size=24, face="bold", color="black"),
                       strip.text.x = element_text(size=24, face="bold", color="black"),
                       strip.text.y = element_text(size=16, face="bold", color="black"),
                       plot.title = element_text(size=24, face="bold", vjust = 1),
                       legend.position="none",
                       panel.background = element_rect( color="black",size=1.2),
                       strip.background = element_rect(color = "black", size = 1.2)) +
        ggplot2::labs(x = "Genomic Position (Mb)",
                      y = expression(log[10](p)),
                      title = plot_traits[i])
    }
    
  }
}