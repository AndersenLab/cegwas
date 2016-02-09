#' Manhattan Plot for GWAS Mapping Data
#'
#' \code{manplot} generates a manhattan plot using \code{ggplot2}
#'
#'
#' @param plot_df the output from the \code{gwas_mappings} function. 
#' @param bf_line_color Set color of bonferroni line.
#' @return Ouput is a ggplot object facetted by chromosome. SNPs above bonferroni corrected p-value (gray line) are colored blue.
#' Confidence interval for a given peak is highlighted in red.
#' @export


manplot <- function(plot_df, bf_line_color = "gray") {
    plot_traits <- unique(plot_df$trait)
    plots <- lapply(plot_traits, function(i) {
        plot_df %>%
        dplyr::filter(trait == i) %>%
        dplyr::distinct(marker) %>%
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
                            color = bf_line_color, 
                            alpha = .75,  
                            size = 1) +
        ggplot2::geom_point( ggplot2::aes(color= factor(aboveBF)) ) +
        ggplot2::facet_grid( . ~ CHROM, scales = "free_x" ) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=14, face="bold", color="black"),
                       axis.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
                       axis.title.x = ggplot2::element_text(size=20, face="bold", color="black", vjust=-.3),
                       axis.title.y = ggplot2::element_text(size=20, face="bold", color="black"),
                       strip.text.x = ggplot2::element_text(size=20, face="bold", color="black"),
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
#' @param use_base Show base at position rather than REF/ALT.
#' @param color_strains character vector containing strains to color in plot. Default is c("N2","CB4856")
#' @return Ouput is a ggplot object facetted by peak SNP (if there are multiple peaks in a mapping). 
#' Genotypes are encoded as REF or ALT and are on the x-axis. Phenotypes are on y-axis.
#' @export

pxg_plot <- function(plot_df, loc = NA, use_base = F, color_strains = c("N2","CB4856")){
  plot_traits <- unique(plot_df$trait)
  plots <- lapply(plot_traits, function(x) {
     plot_peak <- plot_df %>%
      na.omit() %>%
      dplyr::filter(trait == x) %>%
      dplyr::distinct(strain, value, peakPOS) %>%
      dplyr::select(strain, value, CHROM, POS = peakPOS, -allele) %>%
      dplyr::mutate(chr_pos = paste(CHROM, POS, sep="_"))
    
     if (is.na(loc)) {   
      loc <- plot_peak %>% dplyr::select(CHROM, POS) %>% 
        dplyr::distinct() %>% 
        dplyr::transmute(chr_pos = paste0(CHROM, ":",POS))
     }
      
      to_plot <- snpeff(loc[1], severity = "ALL", elements = "ALL") %>%
      dplyr::select(strain, CHROM, POS, GT, REF, ALT) %>%
      dplyr::distinct() %>%
      dplyr::mutate(chr_pos = paste0(CHROM, ":", POS))

      to_plot <- dplyr::left_join(to_plot, plot_peak) %>%
          dplyr::mutate(chr_pos = paste(CHROM, POS, sep=":"))

      to_plot <- dplyr::filter(to_plot, !is.na(value)) %>%
      dplyr::distinct(strain, value, POS) %>%
      dplyr::filter(!is.na(GT)) %>%
      dplyr::group_by(GT) %>%
      dplyr::mutate(GT2 = ifelse(use_base, ifelse(GT == "REF", REF, ALT), GT )) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(GT = GT2)
      
      # Handle Ref/Alt Diffs
      

    
    if (!unique(is.na(color_strains))) {
      to_plot <- to_plot %>%
        dplyr::mutate(colors = ifelse(strain %in% color_strains, strain, "aa_ignore"))%>%
        dplyr::arrange(colors)
    } else {
      to_plot <- to_plot %>%
        dplyr::mutate(colors = "aa_ignore")
    }
    
    to_plot %>%
      ggplot2::ggplot(., ggplot2::aes(x = GT, y = value, fill = GT)) + 
      ggplot2::scale_fill_brewer(palette = "Set1") +
      ggplot2::geom_boxplot(outlier.colour = NA) +
      ggplot2::theme_bw() +
      ggplot2::geom_jitter(ggplot2::aes(color = colors,
                                        size = colors)) +
      ggplot2::scale_color_manual(values = c("black", "#2474FF", "#EE8F03", "red", "green"),
                                  labels = c("Other", 
                                             unique(to_plot$colors)[2], 
                                             unique(to_plot$colors)[3], 
                                             unique(to_plot$colors)[4], 
                                             unique(to_plot$colors)[5])) +
      ggplot2::scale_size_manual(values = c(2, 12, 12, 12, 12),
                                 labels = c("Other", 
                                            unique(to_plot$colors)[2], 
                                            unique(to_plot$colors)[3], 
                                            unique(to_plot$colors)[4], 
                                            unique(to_plot$colors)[5])) +
      # ggplot2::scale_color_brewer(palette = "Dark2", name = "Specified\nStrains") +
      ggplot2::facet_grid(.~chr_pos, scales = "free") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=24, face="bold", color="black"),
                     axis.text.y = ggplot2::element_text(size=24, face="bold", color="black"),
                     axis.title.x = ggplot2::element_text(size=24, face="bold", color="black", vjust=-.3),
                     axis.title.y = ggplot2::element_text(size=24, face="bold", color="black"),
                     strip.text.x = ggplot2::element_text(size=24, face="bold", color="black"),
                     strip.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
                     plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
                     # legend.position="none",
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

#' Plot LD values for significant SNPs in a mapping
#'
#' \code{plot_peak_ld} generates a plot to visualize LD (r) for significant SNPs in a mapping 
#'
#'
#' @param plot_df is the output from the process_mappings function. 
#' @param trait is a string object corresponding to a trait of interest if plot_df has multiple traits in it.
#' @return returns an LDHeatmap object of SNPs that are above the Bonferroni corrected p-value
#' @examples  plot_peak_ld(processed_mapping_df) # for a one trait mapping data frame
#' @examples plot_peak_ld(plot_df = all_maps, trait = "amsacrine_f.l1") # for a multiple trait mapping data frame
#' @export

plot_peak_ld <- function(plot_df, trait = NULL){
  
  
  if (is.null(trait)) {
    snp_df <- plot_df %>% na.omit()
  }
  else {
    snp_df <- dplyr::filter(plot_df, trait == trait) %>% 
      na.omit()
  }
  ld_snps <- dplyr::filter(snps, CHROM %in% snp_df$CHROM, POS %in% 
                             snp_df$POS)
  ld_snps <- data.frame(snp_id = paste(ld_snps$CHROM, ld_snps$POS, 
                                       sep = "_"), data.frame(ld_snps)[, 3:ncol(ld_snps)])
  sn <- list()
  for (i in 1:nrow(ld_snps)) {
    sn[[i]] <- genetics::genotype(as.character(gsub(1, "T/T", 
                                                    gsub(-1, "A/A", ld_snps[i, 2:ncol(ld_snps)]))))
  }
  test <- data.frame(sn)
  colnames(test) <- (ld_snps$snp_id)
  if (ncol(test) == 1) {
    print("Only one significant SNP, not calculating LD")
  }
  else {
    
    ldcalc <- t(genetics::LD(test)[[3]])
    diag(ldcalc) <- 1
    
    
    LDs <- tbl_df(data.frame(ldcalc) %>%
          dplyr::add_rownames(var = "SNP1")) %>%
          tidyr::gather(SNP2, Dprime, -SNP1) %>%
          dplyr::arrange(SNP1) %>%
          tidyr::separate(SNP1, sep = "_", into = c("CHROM1", "POS1"), remove = F) %>%
          dplyr::arrange(CHROM1, as.numeric(POS1))
    
    ldplot <- ggplot2::ggplot(LDs)+
      ggplot2::aes(x = factor(SNP1, levels = SNP1, ordered = T), y = factor(SNP2, levels = SNP1, ordered = T)) +
      ggplot2::geom_tile(ggplot2::aes(fill = Dprime)) +
      ggplot2::geom_text(ggplot2::aes(label = signif(Dprime,3)), fontface = "bold", size = 12)+
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=24, face="bold", color="black"),
                     axis.text.y = ggplot2::element_text(size=24, face="bold", color="black"),
                     axis.title.x = ggplot2::element_text(size=0, face="bold", color="black", vjust=-.3),
                     axis.title.y = ggplot2::element_text(size=0, face="bold", color="black"),
                     strip.text.x = ggplot2::element_text(size=24, face="bold", color="black"),
                     strip.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
                     plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
                     legend.position="none",
                     panel.background = ggplot2::element_rect( color="black",size=1.2),
                     strip.background = ggplot2::element_rect(color = "black", size = 1.2)) +
                     scale_x_discrete(labels = function(x) { gsub("_", ":", x)}, expand = c(0,0)) +
                     scale_y_discrete(labels = function(x) { gsub("_", ":", x)}, expand = c(0,0)) +
                     scale_fill_continuous(high = "#FF0000", low = "white", na.value = "white")
    
    
    
    ldplot <- cowplot::ggdraw(cowplot::switch_axis_position(ldplot, 'y'))
    #     rgb.palette <- grDevices::colorRampPalette(rev(c("blue", 
    #                                                      "orange", "red")), space = "rgb")
    #     ld_outs <- LDheatmap::LDheatmap(test, LDmeasure = "r", 
    #                                     SNP.name = colnames(test), color = rgb.palette(18))
    #     LD.grob1 <- grid::editGrob(ld_outs$LDheatmapGrob, gPath("heatMap", 
    #                                                             "title"), gp = gpar(cex = 1.25, col = "black"))
    #     LD.grob2 <- grid::editGrob(LD.grob1, gPath("geneMap", 
    #                                                "title"), gp = gpar(cex = 0, col = "orange"))
    #     LD.grob3 <- grid::editGrob(LD.grob2, gPath("Key", "title"), 
    #                                gp = gpar(cex = 1.25, col = "black"))
    #     grid::grid.newpage()
    #     grid::grid.draw(LD.grob3)
    return(ldplot)
  }
}


na_to_1 <- function(x) { ifelse(is.na(x), 1, x)}

#' QQ-plot implemented in ggplot2
#'
#' \code{qq_plot} generates a QQ plot given a vector of log10(p) values 
#'
#'
#' @param log10p vector of log10(p) values 
#' @return returns a ggplot2 object 
#' @export

qq_plot <- function(log10p){
  # following four lines from base R's qqline()
  y <- quantile(log10p[!is.na(log10p)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  
  d <- data.frame(resids = log10p)
  
  qqpl <- ggplot2::ggplot(d, ggplot2::aes(sample = resids)) + 
    ggplot2::stat_qq() + 
    ggplot2::geom_abline(slope = slope, intercept = int, color = "red")+
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 24,face = "bold", color = "black"), 
                   axis.text.y = ggplot2::element_text(size = 24, face = "bold", color = "black"), 
                   axis.title.x = ggplot2::element_text(size = 24, face = "bold", color = "black", vjust = -0.3), 
                   axis.title.y = ggplot2::element_text(size = 24,face = "bold", color = "black"), 
                   strip.text.x = ggplot2::element_text(size = 24,face = "bold", color = "black"), 
                   strip.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"), 
                   plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1), 
                   legend.position = "none", 
                   panel.background = ggplot2::element_rect(color = "black", size = 1.2), 
                   strip.background = ggplot2::element_rect(color = "black", size = 1.2)) +
    labs(x = "Theoretical", y = "Observed")
 
  return(qqpl) 
}
