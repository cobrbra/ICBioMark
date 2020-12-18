## code to prepare `example_tib_tables` dataset goes here
example_tib_tables <- get_biomarker_tables(maf = example_maf_data$maf, biomarker = "TIB", sample_list = paste0("SAMPLE_", 1:100))

usethis::use_data(example_tib_tables, overwrite = TRUE)
