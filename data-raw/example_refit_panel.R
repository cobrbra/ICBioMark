## code to prepare `example_refit_panel` dataset goes here
example_refit_panel <- pred_refit_panel(pred_first = example_first_pred_tmb, gene_lengths = example_maf_data$gene_lengths, genes = paste0("GENE_", 1:10))

usethis::use_data(example_refit_panel, overwrite = TRUE)
