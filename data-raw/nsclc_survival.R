## code to prepare `nsclc_survival` dataset goes here
nsclc_survival <- get_clinical_data("nsclc_tcga_broad_2016_sequenced")

usethis::use_data(nsclc_survival, overwrite = TRUE)
