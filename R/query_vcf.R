
get_vcf <- function(remote = F, version = vcf_version) {

  # Set vcf path; determine whether local or remote
  vcf_path <- paste0("~/Dropbox/Andersenlab/Reagents/WormReagents/_SEQ/WI/WI-", vcf_version, "/vcf/WI.", vcf_version, ".snpeff.vcf.gz")
  # Use remote if not available.
  local_or_remote <- "locally"
  if (!file.exists(vcf_path) | remote == T) {
    vcf_path <- paste0("http://storage.googleapis.com/elegansvariation.org/releases/", vcf_version, "/WI.", vcf_version, ".soft-filtered.vcf.gz")
    message("Using remote vcf")
  } else {
    system(paste0("touch ", vcf_path,".csi"))
    message("Using local vcf")
  }
  vcf_path
}



#' Query VCF Data
#'
#' \code{query_vcf} enables you to query variants within a VCF
#'
#' @param ... Gene names, regions, or wormbase identifiers to query.
#' @param info Info columns to output. If an \code{ANN} (annotation) column is available it is automatically fetched. [\strong{Default} \code{c()}]
#' @param format Format columns to output. A \code{"GT"} or \code{"TGT"} column must be specified to retrieve genotypes. \itemize{
#' @param impact A vector of impact levels to filter on (LOW, MODERATE, HIGH, MODIFIER). "ALL" can be used to return ALL variants. [\strong{Default} \code{c('MODERATE', 'HIGH')}]
#'     \item \code{GT} uses a numeric represetnation (0=REF, 1=ALT) and outputs g1, g2, and genotype (0=REF homozygous, 1=HET, 2=ALT homozygous).
#'     \item \code{TGT} uses the base representation (ATGC) and outputs two columns: a1, a2.
#' }
#' [\strong{Default} \code{c("TGT")}]
#' @param samples A set of samples to subset on [\strong{default:} \code{"ALL"}]
#' @param long Return dataset in long or wide format. [\strong{Default} \code{TRUE}]
#' @param vcf Use a custom VCF.
#' @return Dataframe with variant data
#'
#' @examples query_vcf("pot-2","II:1-10000","WBGene00010785")
#' @importFrom dplyr %>%
#' @export

query_vcf <- function(...,
                      info = c(),
                      format = c("TGT"),
                      impact = c("MODERATE", "HIGH"),
                      samples="ALL",
                      long = TRUE,
                      vcf = get_vcf()) {

  regions <- unlist(list(...))

  impact_options <- c(NA,
                      "LOW",
                      "MODERATE",
                      "HIGH",
                      "MODIFIER")

  # Allow user to specify 'ALL'
  if (impact[[1]] == "ALL") {
    impact <- impact_options
  } else {
    # Check that impact is properly specified
    assertthat::assert_that(all(impact %in% impact_options),
                            msg = "Invalid impact specified. Must be LOW, MODERATE, HIGH, or MODIFIER")
  }

  # Ensure that bcftools is available:
  conn <- pipe("bcftools --version")
  bcftools_version <- as.double(stringr::str_extract(readLines(conn)[1], "[0-9]+\\.[0-9]+"))
  close(conn)
  assertthat::assert_that(!(is.na(bcftools_version) | bcftools_version < 1.2),
                          msg = "bcftools 1.2+ required for this query_vcf")


  # Samples
  sample_query <- glue::glue("bcftools query -l {vcf}")
  sample_names <- readr::read_lines(suppressWarnings(pipe(sample_query)))

  # Check that VCF has annotation field
  vcf_header_query <- glue::glue("bcftools view -h {vcf}")
  vcf_header <- readr::read_lines(suppressWarnings(pipe(vcf_header_query)))

  # Pull out info columns
  info_set <- stringr::str_match(vcf_header, '##INFO=<ID=([^,>]+).*Type=([^,>]+).*Description=\"([^,>]+)\"') # nolint
  colnames(info_set) <- c("g", "ID", "Type", "Description")
  info_columns <- purrr::discard(info_set[, 2], is.na)
  info_column_types <- purrr::discard(info_set[, 3], is.na)

  # Add ANN to info if present but not specified.
  if ("ANN" %in% info_columns) {
    info <- c(info, "ANN")
  }

  # Pull out format columns
  format_set <- stringr::str_match(vcf_header, '##FORMAT=<ID=([^,>]+).*Type=([^,>]+).*Description=\"([^,>]+)\"') # nolint
  colnames(format_set) <- c("g", "ID", "Type", "Description")
  format_columns <- c(purrr::discard(format_set[, 2], is.na), "TGT")
  format_column_types <- c(purrr::discard(format_set[, 3], is.na), "TGT")

  if (is.null(regions)) {
    # Begin Exclude Linting
    cat(crayon::bold("\nVCF:"), vcf, "\n")
    cat(crayon::bold("VCF Samples:"), length(sample_names), "\n\n")
    cat(crayon::bold("Info\n"))
    utils::write.table(
      dplyr::tbl_df(
        info_set[, 2:4]) %>%
        dplyr::filter(stats::complete.cases(.)),
      sep = "\t",
      row.names = F,
      quote = F)
    cat("\n")
    cat(crayon::bold("Format"), "\n")
    utils::write.table(
      dplyr::tbl_df(
        format_set[, 2:4]) %>%
        dplyr::filter(stats::complete.cases(.)),
      sep = "\t",
      row.names = F,
      quote = F)
    cat("\n")
    return(invisible(NULL))
    # End Exclude Linting
  }


  assertthat::assert_that(all(samples %in% sample_names) | (samples[[1]] == "ALL"),
                          msg = {
                            unknown_samples <- paste(samples[!samples %in% sample_names], collapse = ", ")
                            glue::glue("Specified samples are not in the VCF: {unknown_samples}")
                          })


  assertthat::assert_that(all(info %in% info_columns),
                          msg = {
                            missing_info_columns <- paste0(info[!(info %in% info_columns)], collapse = ",")
                            glue::glue("VCF does not have specified info fields: {missing_info_columns}")
                          })

  assertthat::assert_that(all(format %in% format_columns),
                          msg = {
                            missing_format_columns <- paste0(format[!(format %in% format_columns)], collapse = ",")
                            glue::glue("VCF does not have specified format fields: {missing_format_columns}")
                          })


  results <- suppressWarnings(lapply(regions, function(query) {
    # Save region as query

    # Fix region specifications
    query <- gsub("\\.\\.", "-", query)
    query <- gsub(",", "", query)

    # Resolve region names
    if (!grepl("(I|II|III|IV|V|X|MtDNA).*", query)) {
      elegans_gff <- get_db(table = "wormbase_gene")
      # Pull out regions by element type.
      region <- dplyr::collect(dplyr::filter(elegans_gff,
                                             locus == query |
                                               gene_id == query |
                                               transcript_id == query) %>%
                                 dplyr::select(chrom, start, end, gene_id, locus,  exon_id, transcript_id, transcript_biotype) %>%
                                 dplyr::distinct(.keep_all = TRUE)) %>%
        dplyr::summarize(chrom = chrom[1], start = min(start), end = max(end)) %>%
        dplyr::mutate(region_format = paste0(chrom, ":", start, "-", end)) %>%
        dplyr::select(region_format) %>%
        dplyr::distinct(.keep_all = TRUE)

      region <- paste(region$region_format, collapse = ",")

      query_type <- "locus"

      if (stringr::str_length(regions[[1]]) == 0) {
        message(paste0(query, " not found."))
        region <- NA
      }
    } else {
      region <- query
      query_type <- "region"
    }

    base_header <- c("CHROM",
                     "POS",
                     "REF",
                     "ALT",
                     "FILTER")

    if ("ANN" %in% info) {
      ann_header <- c("allele",
                      "effect",
                      "impact",
                      "gene_name",
                      "gene_id",
                      "feature_type",
                      "feature_id",
                      "transcript_biotype",
                      "exon_intron_rank",
                      "nt_change",
                      "aa_change",
                      "cdna_position_or_len",
                      "protein_position",
                      "distance_to_feature",
                      "error",
                      "extra")
    } else {
      ann_header <- c()
    }

    info_query <- paste0(info, collapse = "\\t%")
    format_query <- paste0(format, collapse = "!%")
    # If using long format provide additional information.
    query_string <- glue::glue("'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%{info_query}[\\t%{format_query}]\\n'")
    if (long == T) {
      query_string <- glue::glue("'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%{info_query}[\\t%{format_query}]\\n'")
    }

    if (samples != "ALL") {
      sample_query <- glue::glue("--samples ", paste(samples, collapse = ","))
    } else {
      sample_query <- ""
      samples <- sample_names
    }

    # Grep impacts to speed up intake
    if (impact != "ALL" & !is.na(impact)) {
      impact_grep <- paste(purrr::discard(impact, is.na), collapse = "|")
      impact_grep <- glue::glue("| egrep \"({impact_grep})\" - ")
    } else {
      impact_grep <- ""
    }

    output_file <- tempfile()
    command <- paste("bcftools",
                     "query",
                     sample_query,
                     "--regions",
                     region,
                     "-f",
                     query_string,
                     vcf,
                     impact_grep,
                     ">",
                     output_file)
    if (!is.na(region)) {
      message(glue::glue("Query: {query}"))
      conn <- system(command)
      result <- try(dplyr::tbl_df(data.table::fread(output_file,
                                                    col.names = c(base_header,
                                                                  info,
                                                                  samples),
                                                    sep = "\t",
                                                    na.strings = c("", "."))),
                    silent = FALSE)
      try(file.remove(output_file))
      if (!grepl("^Error.*", result[[1]][1])) {
        tsv <- result %>%
          dplyr::mutate(REF = ifelse(REF == TRUE, "T", REF), # T nucleotides are converted to 'true'
                        ALT = ifelse(ALT == TRUE, "T", ALT))
      } else {
        tsv <- as.data.frame(NULL)
      }

      # If no results are returned, stop.
      if (typeof(tsv) == "character" | nrow(tsv) == 0) {
        warning("No Variants")
        NA
      } else {
        # If no ANN column exists; Create a dummy one to unnest on.
        if (!("ANN" %in% colnames(tsv))) {
          tsv <- dplyr::mutate(tsv, ANN = "")
        }
        tsv <- tsv %>%
          tidyr::unnest(ANN, .sep = ",") %>%
          {
            if (!is.null(ann_header))
              tidyr::separate(., ANN, into = ann_header, sep = "\\|") %>%
              dplyr::select(dplyr::one_of(c(base_header, ann_header)), dplyr::everything(), -extra)
            else
              dplyr::select(., dplyr::one_of(c(base_header)), dplyr::everything())
          } %>%
          dplyr::mutate(query = query, region = region) %>%
          dplyr::select(CHROM, POS, query, region, dplyr::everything())

        # For locus queries, filter out non-matching genes.
        if (query_type == "locus") {
          tsv <- dplyr::filter(tsv,
                               (query == gene_id) | (query == gene_name) | (query == feature_id))
        }

        # Filter impact (only if ANN present)
        if ("impact" %in% names(tsv)) {
          tsv <- tsv[tsv$impact %in% impact, ]
        } else if (impact != "ALL" | !is.na(impact)) {
          warning("Warning: No ANN column specified; Variants will not be filtered on impact.")
        }

        if (nrow(tsv) == 0) {
          message(paste("No Results for", region, "after filtering"))
        }

        if (long == FALSE) {
          tsv
        } else {
          tsv <- tidyr::gather_(tsv, "SAMPLE", "FORMAT_COLUMN", samples)  %>%
            tidyr::separate(FORMAT_COLUMN,
                            into = format,
                            sep = "\\!",
                            convert = TRUE,
                            remove = T) %>%
                            {
                              if ("DP" %in% format) dplyr::mutate(., DP = as.integer(ifelse( (DP == ".") | is.na(DP), 0, DP))) else .
                            } %>%
                            {
                              if ("TGT" %in% format)
                                tidyr::separate(.,
                                                TGT,
                                                sep = "\\/|\\|",
                                                into = c("a1", "a2"),
                                                convert = T,
                                                remove = T)
                              else
                                .
                            } %>%
                            {
                              if ("GT" %in% format)
                                tidyr::separate(.,
                                                GT,
                                                sep = "\\/|\\|",
                                                into = c("g1", "g2"),
                                                convert = T,
                                                remove = T) %>%
                                dplyr::mutate_at(as.integer, .vars = c("g1", "g2")) %>%
                                # Why this weird way? It's faster.
                                dplyr::mutate(genotype = as.integer(rowSums(.[, c("g1", "g2")])))
                              else
                                .
                            }

          # Fix NAs and convert columns
          tsv <- tsv %>%
            dplyr::mutate_all(function(x) { ifelse((x == ".") | (x == ""), NA, x)}) %>%
            # Info columns
            dplyr::mutate_at(.vars = info_columns[info_column_types == "Integer" & info_columns %in% info],
                             as.integer) %>%
            dplyr::mutate_at(.vars = info_columns[info_column_types == "Float" & info_columns %in% info],
                             as.numeric) %>%
            dplyr::mutate_at(.vars = info_columns[info_column_types == "Flag" & info_columns %in% info],
                             as.logical) %>%
            # Format Columns
            dplyr::mutate_at(.vars = format_columns[format_column_types == "Integer" & format_columns %in% format],
                             as.integer) %>%
            dplyr::mutate_at(.vars = format_columns[format_column_types == "Float" & format_columns %in% format],
                             as.numeric) %>%
            dplyr::mutate_at(.vars = format_columns[format_column_types == "Flag" & format_columns %in% format],
                             as.logical)

          column_order <- c("CHROM",
                            "POS",
                            "REF",
                            "ALT",
                            "SAMPLE",
                            "FILTER",
                            "FT",
                            "a1",
                            "a2",
                            "g1",
                            "g2",
                            "genotype",
                            "query",
                            "region")
          column_order_use <- c(column_order[column_order %in% names(tsv)],
                                names(tsv)[!names(tsv) %in% column_order])

          tsv <- tsv %>%
            dplyr::select_at(.vars = column_order_use) %>%
            dplyr::arrange(CHROM, POS)

          # Remove ANN field if its empty
          if (!("ANN" %in% info_columns)) {
            tsv <- tsv %>% dplyr::select(-ANN)
          }
        }
        tsv
      }
    }
  }
  ))
  results <- do.call(rbind, results)
  if (!"CHROM" %in% names(results)){
    warning("No Results")
    return(NA)
  }
  results <- results %>% dplyr::filter(!is.na(CHROM))
  results
}
