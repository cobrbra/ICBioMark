## code to prepare `ensembl_gene_lengths` dataset goes here

ensembl_mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes <- biomaRt::getBM(attributes = "hgnc_symbol", filters = "chromosome_name", values = c(paste(1:22), "X", "Y"), mart = ensembl_mart, useCache = FALSE)

ensembl_gene_lengths <- biomaRt::getBM(attributes = c("chromosome_name", "cds_length", "ensembl_gene_id"), filters = "hgnc_symbol", values = genes$hgnc_symbol, mart = ensembl_mart, useCache = FALSE) %>%
  {dplyr::full_join(., biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"), filters = "ensembl_gene_id", values = .$ensembl_gene_id, mart = ensembl_mart, useCache = FALSE), by = "ensembl_gene_id")} %>%
  dplyr::filter(chromosome_name %in% c(paste(1:22), "X", "Y")) %>%
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::mutate(max_cds = max(cds_length, na.rm = TRUE)) %>%
  dplyr::filter(is.finite(max_cds)) %>%
  dplyr::filter(!is.na(max_cds)) %>%
  dplyr::ungroup() %>%
  dplyr::select(hgnc_symbol, max_cds, chromosome_name) %>%
  dplyr::distinct() %>%
  dplyr::mutate(Hugo_Symbol = hgnc_symbol) %>%
  dplyr::mutate(Chromosome = chromosome_name) %>%
  dplyr::select(Hugo_Symbol, max_cds, Chromosome)

rownames(ensembl_gene_lengths) <- ensembl_gene_lengths$Hugo_Symbol

usethis::use_data(ensembl_gene_lengths, overwrite = TRUE)
