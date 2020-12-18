context("Inputting MAFs, outputting mutation matrices.")

test_that("A simple call to get_table_from_maf() has the right structure and value", {
  output <- get_table_from_maf(maf = example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100))

  expect_is(output, "list")
  expect_equal(length(names(output)), 5)
  expect_equal(names(output), c("matrix", "sample_list", "gene_list", "mut_types_list", "col_names"))
  expect_known_output(output, file = "datatest_table.rds")
})

test_that("A simple call to get_mutation_tables() has the right structure and value", {
  output <- get_mutation_tables(maf = example_maf_data$maf, sample_list = paste0("SAMPLE_", 1:100), gene_list = example_maf_data$gene_lengths$Hugo_Symbol)

  expect_is(output, "list")
  expect_equal(length(names(output)), 3)
  expect_equal(names(output), c("train", "val", "test"))

  train_data <- output$train
  expect_is(train_data, "list")
  expect_equal(length(names(train_data)), 5)
  expect_equal(names(train_data), c("matrix", "sample_list", "gene_list", "mut_types_list", "col_names"))

  expect_is(train_data$matrix, "dgCMatrix")

  expect_equal(output, example_tables)
})
