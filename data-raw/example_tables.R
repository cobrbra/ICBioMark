## code to prepare `example_tables` dataset goes here
example_tables <- get_mutation_tables(example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100), gene_list = paste0("GENE_", 1:20))

usethis::use_data(example_tables, overwrite = TRUE)
