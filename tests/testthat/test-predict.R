context("Functions for fitting predictive models")

test_that("pred_first_fit gives reference value", {
  output <- pred_first_fit(example_gen_model, lambda = exp(seq(-9, -14, length.out = 100)),
                           training_matrix = example_tables$train$matrix,
                           gene_lengths = example_maf_data$gene_lengths)
  expect_equal(output, example_first_pred_tmb)
})
