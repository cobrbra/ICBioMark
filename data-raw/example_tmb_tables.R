## code to prepare `example_tmb_tables` dataset goes here
example_tmb_tables <- get_biomarker_tables(maf = example_maf_data$maf, biomarker = "TMB", sample_list = paste0("SAMPLE_", 1:100))

usethis::use_data(example_tmb_tables, overwrite = TRUE)
