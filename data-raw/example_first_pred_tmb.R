## code to prepare `example_first_pred_tmb` dataset goes here
example_first_pred_tmb <- pred_first_fit(example_gen_model, lambda = exp(seq(-9, -14, length.out = 100)),
                                   training_matrix = example_tables$train$matrix,
                                   gene_lengths = example_maf_data$gene_lengths)

usethis::use_data(example_first_pred_tmb, overwrite = TRUE)
