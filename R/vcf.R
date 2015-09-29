#' Generate Kinship
#'
#' \code{generate_kinship} uses the rrBLUP package function \code{A.mat} to generate a kinship matrix for mapping.
#' Requires bcftools. Heterozygous calls and missing genotypes are ignored.
#'
#' @param vcf - a bcf/vcf/vcf.gz file
#' @return A kinship matrix. 
#' @export


generate_kinship <- function(vcf) {
  gt_file <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = "")
  command <- paste0("bcftools view  -m2 -M2  ", vcf," | ",
                    "bcftools query --print-header -f '%CHROM\\t%POS[\\t%GT]\n' |",
                    "sed 's/0|0/0/g'   | ",
                    "sed 's/1|1/1/g'   | ",
                    "sed 's/0|1/NA/g'  | ",
                    "sed 's/1|0/NA/g'  | ",
                    "sed 's/.|./NA/g'  | ",
                    "sed 's/0\\/0/0/g'  | ",
                    "sed 's/1\\/1/1/g'  | ",
                    "sed 's/0\\/1/NA/g' | ",
                    "sed 's/1\\/0/NA/g' | ",
                    "sed 's/.\\/./NA/g'  > ", gt_file)
  system(command, intern = T)
  df <- tbl_df(fread(gt_file))
  names(df) <- gsub(":GT","",gsub("[\\# ]{0,2}\\[[0-9]+\\]","",names(df)))
  
  # Fix position and unite columns.
  df <- mutate(df, POS=as.character(POS)) %>%
    unite(CHROM_POS, CHROM, POS, sep = "_")
  row.names(df) <- df$CHROM_POS
  df <- select(df, -CHROM_POS)
  kinship <- A.mat(t(data.matrix(df)))
}




