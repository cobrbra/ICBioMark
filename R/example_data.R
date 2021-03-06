#' Simulated MAF Data
#'
#' @description An example dataset generated by the function generate_maf_data(), with n_sample = 100 and n_genes = 20.
#'
#' @format A list with two entries:
#' \describe{
#'   \item{maf}{An annotated mutation dataframe with 3 columns and 1346 rows: \describe{
#'     \item{Tumor_Sample_Barcode}{A sample id for each mutation.}
#'     \item{Hugo_Symbol}{The name of the gene location for each mutation.}
#'     \item{Variant_Classification}{The mutation type for each mutation.}}}
#'   \item{gene_lengths}{A data frame with two rows: \describe{
#'     \item{Hugo_Symbol}{The name of each gene.}
#'     \item{max_cds}{The length of each gene, as defined by maximum coding sequence.}}}
#' }
"example_maf_data"

#' Mutation Matrices from Simulated Data
#'
#' @description Mutation data extracted from the pre-loaded example mutation data example_maf_data, using the function get_mutation_tables().
#'
#' @format A list with three entries:
#' \describe{
#'   \item{train}{An object 'train'.}
#'   \item{val}{An object 'val'.}
#'   \item{test}{An object 'test'.}
#' }
#' Each of these three objects is a list with the following entries (for more detail see the documentation for the function get_table_from_maf()):
#' \describe{
#'   \item{matrix}{A sparse matrix of mutations.}
#'   \item{sample_list}{A character vector of sample IDs, corresponding to the rows of the mutation matrix.}
#'   \item{gene_list}{A character vector of gene names.}
#'   \item{mut_types_list}{A character vector of mutation types.}
#'   \item{colnames}{A character vector of gene name/mutation type combinations (in each case separated by the character "_"), corresponding to the columns of the mutation matrix.}
#' }
"example_tables"

#' Generative Model from Simulated Data
#'
#' @description An example of the output produced by fit_gen_model() on simulated data.
#'
#' @format A list with two entries:
#' \describe{
#'   \item{fit}{A glmnet fit object.}
#'   \item{dev}{A table containing the average deviance of each cross-validation fold, for each penalisation factor in fit$lambda.}
#'   \item{s_min}{The index of the regularisation penalty minimising average deviance across folds.}
#' }
"example_gen_model"

#' Tumour Mutation Burden of Example Train, Validation and Test Data.
#'
#' @description An example output produced by using the function get_biomarker_tables(), applied to the example MAF data pre-loaded in example_maf_data$maf.
#'
#' @format A list with threeobjects: 'train', 'val' and 'test'. Each is a dataframe with two columns:
#' \describe{
#'   \item{Tumor_Sample_Barcode}{A unique ID for each sample.}
#'   \item{TMB}{The value of Tumour Mutation Burden for that sample.}
#' }
"example_tmb_tables"

#' Tumour Indel Burden of Example Train, Validation and Test Data.
#'
#' @description An example output produced by using the function get_biomarker_tables(), applied to the example MAF data pre-loaded in example_maf_data$maf.
#'
#' @format A list with threeobjects: 'train', 'val' and 'test'. Each is a dataframe with two columns:
#' \describe{
#'   \item{Tumor_Sample_Barcode}{A unique ID for each sample.}
#'   \item{TIB}{The value of Tumour Indel Burden for that sample.}
#' }
"example_tib_tables"

#' First-Fit Predictive Model Fitting on Example Data
#'
#' @description An example output from the function pred_first_fit(), applied to pre-loaded example mutation data.
#'
#' @format A list with six entries:
#' \describe{
#'   \item{fit}{A gglasso fit.}
#'   \item{panel_genes}{A matrix where each row corresponds to a gene, each column to an iteration of the group
#'   lasso with a different penalty factor, and the elements booleans specifying whether that gene was selected
#'   to be included in that iteration.}
#'   \item{panel_lengths}{A vector giving total panel length for each gglasso iteration.}
#'   \item{p}{The vector of weights used in the optimisation procedure.}
#'   \item{K}{The bias penalty factor used in the optimisation procedure.}
#'   \item{names}{Gene and mutation type information as used when fitting the generative model.}
#'
#' }
"example_first_pred_tmb"

#' Refitted Predictive Model Fitted on Example Data
#'
#' @description An example output from use of the function pred_refit_panel(), applied to example gene length data and generative model fit.
#'
#' @format A list with three entries:
#' \describe{
#'   \item{fit}{A list with a single element 'beta', a matrix with prediction weights.}
#'   \item{panel_genes}{A matrix (in this case with a single column) where each row corresponds to a gene,
#'   and each entry corresponds to whether the gene is included in the panel.}
#'   \item{panel_lengths}{A vector of length 1 giving total panel length.}
#' }
"example_refit_panel"

#' Refitted Predictive Models Fitted on Example Data
#'
#' @description An example output from use of the function pred_refit_range(), applied to example gene length data and generative model fit.
#'
#' @format A list with six entries:
#' \describe{
#'   \item{fit}{A list with a single element 'beta', a matrix with prediction weights.}
#'   \item{panel_genes}{A matrix where each row corresponds to a gene,
#'   and each entry corresponds to whether the gene is included in the panel.}
#'   \item{panel_lengths}{A vector giving total panel length.}
#' }
"example_refit_range"

#' Example Predictions
#'
#' @description An example output from use of the function get_predictions(), applied to the pre-loaded datasets example_refit_range and example_tables$val .
#'
#' @format A list with two entries:
#' \describe{
#'   \item{predictions}{A a matrix containing a row for each sample and a column for each panel.}
#
#'   \item{panel_lengths}{A vector giving total panel lengths.}
#' }
"example_predictions"
