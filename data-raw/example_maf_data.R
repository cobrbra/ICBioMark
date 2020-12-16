## code to prepare `example_maf_data` dataset goes here
example_maf_data <- generate_maf_data(n_samples = 100, n_genes = 20)

usethis::use_data(example_maf_data, overwrite = TRUE)
