## code to prepare `example_tables` dataset goes here
example_tables <- get_tables(example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100), gene_list = example_maf_data$gene_lengths$Hugo_Symbol)

usethis::use_data(example_tables, overwrite = TRUE)
