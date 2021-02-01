## code to prepare `example_predictions` dataset goes here
example_predictions <- get_predictions(example_refit_range, new_data = example_tables$val)

usethis::use_data(example_predictions, overwrite = TRUE)
