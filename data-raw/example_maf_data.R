## code to prepare `example_maf_data` dataset goes here
example_maf_data <- generate_maf_data(n_samples = 20, n_genes = 50)

usethis::use_data(example_maf_data, overwrite = TRUE)
