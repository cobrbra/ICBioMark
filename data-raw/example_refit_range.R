## code to prepare `example_refit_range` dataset goes here
example_refit_range <- pred_refit_range(pred_first = example_first_pred_tmb, gene_lengths = example_maf_data$gene_lengths)

usethis::use_data(example_refit_range, overwrite = TRUE)
