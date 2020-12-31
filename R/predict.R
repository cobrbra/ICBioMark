#' Construct Optimisation Parameters.
#'
#' From the learned generative model and training data, produces a vector of weights p to be used in
#' the subsequent group lasso optimisation, alongside a biomarker-dependent normalisation quantity p_norm.
#'
#' @param gen_model (list)
#' A generative mutation model, fitted by fit_gen_model().
#' @param training_matrix (sparse matrix)
#' A sparse matrix of mutations in the training dataset, produced by get_mutation_tables().
#' @param marker_mut_types (character)
#' A character vector listing which mutation types (of the set specified in the generative model
#' attribute 'names') constitute the biomarker in question.
#' @param gene_lengths (dataframe)
#' A table with two columns: Hugo_Symbol and max_cds, providing the lengths of the genes to be modelled.
#'
#' @return
#' A list with three entries:
#'  * A vector p, with an entry corresponding to each combination of gene and mutation type specified
#'   in the generative model fitted. Each component is a non-negative value corresponding to a weighting
#'   p to be supplied to a group lasso optimisation.
#'  * A numeric p_norm, giving the factor between p_{gs} and þ_{0gs} (see paper for details).
#'  * A vector biomarker_columns, detailing which of the elements of p correspond to gene/mutation type
#'   combinations contributing to the biomarker in question.
#'
#' @export
#'
#' @examples
#' p <- get_p(example_gen_model, example_tables$train$matrix,
#'            marker_mut_types = c("I"), gene_lengths = example_maf_data$gene_lengths)
#' print(p$p[1:5])
#' print(p$p_norm)
#' print(p$bc[1:5])


get_p <- function(gen_model, training_matrix, marker_mut_types, gene_lengths) {

  n_samples <- nrow(training_matrix)
  n_genes <- length(gen_model$names$gene_list)
  n_mut_types <- length(gen_model$names$mut_types_list)

  mutation_vector <- training_matrix
  dim(mutation_vector) <- c(n_samples*n_genes*n_mut_types, 1)
  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol

  biomarker_ids <- purrr::map(marker_mut_types, ~which(gen_model$names$mut_types_list == .))
  biomarker_columns <- sort(unlist(purrr::map(biomarker_ids, ~seq(., n_genes*n_mut_types, n_mut_types))))

  t_s_getter <- Matrix::sparseMatrix(i = rep(1:n_mut_types, times = n_genes, each = n_samples),
                                     j = 1:(n_samples*n_genes*n_mut_types),
                                     dims = c(n_mut_types, n_samples*n_genes*n_mut_types))
  t_s <- as.vector(t_s_getter %*% mutation_vector)
  t_s_getter <- NULL

  weights <- rep(t_s, times = n_genes)
  weights <- weights*rep(gene_lengths[gen_model$names$gene_list,]$max_cds, each = n_mut_types)
  weights <- weights/mean(weights)

  i <- 1:(n_genes*n_mut_types)
  j <- rep(1:(n_genes*n_mut_types))
  j_ <- rep(seq(1, n_genes*n_mut_types, n_mut_types), each = n_mut_types)

  param_getter <- Matrix::sparseMatrix(c(i, i), c(j, j_))
  param_est <- as.vector(param_getter %*% gen_model$fit$beta[, gen_model$s_min])
  names(param_est) <- gen_model$names$col_names
  param_getter <- NULL

  p_norm <- sum((weights * exp(gen_model$fit$a0[gen_model$s_min] + param_est))[biomarker_columns])
  p <- weights * exp(gen_model$fit$a0[gen_model$s_min] + param_est) / p_norm
  names(p) <- gen_model$names$col_names

  return(list(p = p, p_norm = p_norm, bc = biomarker_columns))
}

#' Construct Bias Penalisation
#'
#' @param gen_model (list)
#' A generative mutation model, fitted by fit_gen_model().
#' @param p_norm (numeric)
#' Scaling factor between coefficients of p and parameters of generative model (see paper for details).
#' @param training_matrix (sparse matrix)
#' A sparse matrix of mutations in the training dataset, produced by get_mutation_tables().
#' @param marker_training_values (dataframe)
#' A dataframe containing training values for the biomarker in question.
#' @param method (function)
#' How to select a representative biomarker value from the training dataset. Defaults to max().
#'
#' @return
#' A numerical value, to be used as a penalty weighting in the subsequent group lasso optimisation.
#' @export
#'
#' @examples
#' K <- get_K(example_gen_model, 1, example_tables$train$matrix)
#' print(K)

get_K <- function(gen_model, p_norm, training_matrix, marker_training_values = NULL, method = max) {

  if (is.null(marker_training_values)) {
    marker_training_values <- data.frame(Tumor_Sample_Barcode = gen_model$names$sample_list,
                                         value = Matrix::rowSums(training_matrix),
                                         stringsAsFactors = FALSE)
  }

  K <- exp(gen_model$fit$a0[[gen_model$s_min]])*method(marker_training_values$value)*p_norm
  return(K)

}

#' Extract Panel Details from Group Lasso Fit
#'
#' @param gene_lengths (dataframe)
#' A table with two columns: Hugo_Symbol and max_cds, providing the lengths of the genes to be modelled.
#' @param fit (list)
#' A fit from the group lasso algorithm, produced by the function gglasso (package: gglasso).
#' @param gene_list (character)
#' A character vector of genes listing the genes (in order) included in the model pred_fit.
#' @param mut_types_list (character)
#' A character vector listing the mutation type groupings (in order) included in the model pred_fit.
#'
#' @return
#' A list of two elements:
#'  * panel_genes: A matrix where each row corresponds to a gene, each column to an iteration of the group
#'   lasso with a different penalty factor, and the elements booleans specifying whether that gene was selected
#'   to be included in that iteration.
#'  * panel_lengths:
#' @export
#'
#' @examples
#' panels <- get_panels_from_fit(example_maf_data$gene_lengths, example_pred_fit$fit,
#' example_gen_model$names$gene_list, mut_types_list = example_gen_model$names$mut_types_list)
#'
#' print(panels$fit)

get_panels_from_fit <- function(gene_lengths, fit, gene_list, mut_types_list) {

  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol
  n_genes <- length(gene_list)
  n_mut_types <- length(mut_types_list)
  col_names <- paste0(rep(gene_list, each = n_mut_types), "_", rep(mut_types_list, times = n_genes))
  pred_fit_abs <- as.data.frame(abs(fit$beta))
  panel_inclusion <- dplyr::group_by(pred_fit_abs, gene = factor(rep(gene_list,
                                            each = length(mut_types_list)),
                                        levels = gene_list))
  panel_genes <- dplyr::select(dplyr::summarise(tbl = panel_inclusion, dplyr::across(TRUE, sum)), setdiff(colnames(panel_inclusion), "gene")) != 0
  rownames(panel_genes) <- gene_list
  panel_gene_lengths <- gene_lengths[gene_list,]$max_cds

  panel_lengths <- as.vector(t(panel_gene_lengths) %*% panel_genes)

  return(list(panel_genes = panel_genes, panel_lengths = panel_lengths))

}

#' First-Fit Predicitve Model with Group Lasso
#'
#' This function implements the first-fit
#'
#' @param gen_model (list)
#' A generative mutation model, fitted by fit_gen_model().
#' @param lambda (numeric)
#' A vector of penalisation weights for input to the group lasso optimiser gglasso.
#' @param biomarker (character)
#' The biomarker in question. If "TMB" or "TIB", then automatically defines the subsequent
#' varaible marker_mut_types.
#' @param marker_mut_types (character)
#' The set of mutation type groupings constituting the biomarker being estimated. Should be
#' a vector comprising of elements of the mut_types_list vector in the 'names' attribute of
#' gen_model.
#' @param training_matrix (sparse matrix)
#' A sparse matrix of mutations in the training dataset, produced by get_mutation_tables().
#' @param gene_lengths (dataframe)
#' A table with two columns: Hugo_Symbol and max_cds, providing the lengths of the genes to be modelled.
#' @param marker_training_values (dataframe)
#' A dataframe containing two columns: 'Tumor_Sample_Barcode', containing the sample IDs for the training
#' dataset, and a second column containing training values for the biomarker in question.
#' @param K_method (function)
#' How to select a representative biomarker value from the training dataset. Defaults to max().
#' @param free_genes (character)
#' Which genes should escape penalisation (for example when augmenting a pre-existing panel).
#'
#' @return
#' A list of six elements:
#'  * fit: Output of call to gglasso.
#'  * panel_genes: A matrix where each row corresponds to a gene, each column to an iteration of the group
#'   lasso with a different penalty factor, and the elements booleans specifying whether that gene was selected
#'   to be included in that iteration.
#'  * panel_lengths: A vector giving total panel length for each gglasso iteration.
#'  * p: The vector of weights used in the optimisation procedure.
#'  * K: The bias penalty factor used in the optimisation procedure.
#'  * names: Gene and mutation type information as used when fitting the generative model.
#'
#' @export
#'
#' @examples
#' example_first_fit <- pred_first_fit(example_gen_model, lambda = exp(seq(-9, -14, length.out = 100)),
#'                                     training_matrix = example_tables$train$matrix,
#'                                     gene_lengths = example_maf_data$gene_lengths)
#'
pred_first_fit <- function(gen_model, lambda = exp(seq(-16,-24, length.out = 100)), biomarker = "TMB",
                           marker_mut_types = c("NS", "I"), training_matrix, gene_lengths, marker_training_values = NULL,
                           K_method = max, free_genes = c()) {

  if (biomarker == "TIB") {
    marker_mut_types <- c("I")
  }

  wrong_mutation_types <- setdiff(marker_mut_types, gen_model$names$mut_types_list)
  if (length(wrong_mutation_types) > 0) {
    stop(paste0("Mutation types ", paste(wrong_mutation_types, collapse = ", ", " not in generative model.")))
  }

  n_samples <- nrow(training_matrix)
  n_genes <- length(gen_model$names$gene_list)
  n_mut_types <- length(gen_model$names$mut_types_list)

  p <- get_p(gen_model = gen_model, training_matrix = training_matrix, marker_mut_types = marker_mut_types,
             gene_lengths = gene_lengths)

  K <- get_K(gen_model = gen_model, p_norm = p$p_norm, training_matrix = training_matrix,
             marker_training_values = marker_training_values, method = K_method)

  rownames(gene_lengths) <- gene_lengths$Hugo_Symbol
  gene_lengths[free_genes,'max_cds'] <- 0
  pf <- gene_lengths[gen_model$names$gene_list,]$max_cds

  X <- matrix(0, n_mut_types * n_genes + 1, n_mut_types * n_genes)
  colnames(X)

  X[1,] <- sqrt(K)*p$p
  diag(X[2:(n_mut_types * n_genes + 1),]) <- sqrt(p$p)
  Y <- c(sqrt(K), sqrt(p$p))
  Y[setdiff(2:(n_mut_types * n_genes + 1), 1 + p$bs)] <- 0

  fit <- gglasso::gglasso(x = X, y = Y, loss = "ls", lambda = lambda,
                          group = rep(1:n_genes, each = n_mut_types), intercept = FALSE, pf = pf)
  rownames(fit$beta) <- gen_model$names$col_names

  panels <- get_panels_from_fit(gene_lengths = gene_lengths, fit = fit,
                                gene_list = gen_model$names$gene_list, mut_types_list = gen_model$names$mut_types_list)

  return(list(fit = fit, panel_genes = panels$panel_genes, panel_lengths = panels$panel_lengths, p = p$p, K = K, names = gen_model$names))
}

pred_refit_panel <- function(gen_model) {

}

pred_refit_range <- function(gen_model, lamda = exp(seq(-16, -24, length.out = 100))) {

}

pred_refit_size <- function(gen_model, size = 1, pred_range = NULL) {


}
