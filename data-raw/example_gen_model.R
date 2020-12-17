## code to prepare `example_gen_model` dataset goes here
example_gen_model <- fit_gen_model(example_maf_data$gene_lengths, table = example_tables$train)

usethis::use_data(example_gen_model, overwrite = TRUE)
