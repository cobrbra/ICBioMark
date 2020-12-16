## code to prepare `example_tables` dataset goes here
example_tables <- get_tables(example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100))

usethis::use_data(example_tables, overwrite = TRUE)
