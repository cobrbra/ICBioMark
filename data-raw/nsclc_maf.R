## code to prepare `nsclc_maf` dataset goes here

nsclc_maf <- fetch_all_tcgadata(case_id = "nsclc_tcga_broad_2016_sequenced", gprofile_id = "nsclc_tcga_broad_2016_mutations", glist = ensembl_gene_lengths$Hugo_Symbol, mutations = TRUE) %>%
  dplyr::select(gene_symbol, mutation_type, case_id, chr, start_position, end_position) %>%
  dplyr::mutate(Tumor_Sample_Barcode = case_id,
                Hugo_Symbol = gene_symbol,
                Variant_Classification = mutation_type,
                Start_Position = start_position,
                End_Position = end_position,
                Chromosome = chr) %>%
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification, Chromosome, Start_Position, End_Position)

usethis::use_data(nsclc_maf, overwrite = TRUE)
