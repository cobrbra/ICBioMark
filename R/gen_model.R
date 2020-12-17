#' Fit Generative Model
#'
#' @param gene_lengths (dataframe)
#' A table with two columns: Hugo_Symbol and max_cds, providing the lengths of the genes to be modelled.
#' @param matrix (Matrix::sparseMatrix)
#' A mutation matrix, such as produced by the function get_table_from_maf().
#' @param sample_list (character)
#' The set of samples to be modelled.
#' @param gene_list (character)
#' The set of genes to be modelled.
#' @param mut_types_list (character)
#' The set of mutation types to be modelled.
#' @param col_names (character)
#' The column names of the 'matrix' parameter.
#' @param table (list)
#' Optional parameter combining matrix, sample_list, gene_list, mut_types_list, col_names, as is produced by the function get_tables().
#' @param nlambda (numeric)
#' The length of the vector of penalty weights, passed to the function glmnet::glmnet().
#' @param n_folds (numeric)
#' The number of cross-validation folds to employ.
#' @param maxit (numeric)
#' Technical parameter passed to the function glmnet::glmnet().
#' @param seed_id (numeric)
#' #' Input value for the function set.seed().
#'
#'
#' @return A list comprising three objets:
#' * An object 'fit', a fitted glmnet model.
#' * A table 'dev', giving average deviances for each regularisation penalty factor and cross-validation fold.
#' * An integer 's_min', the index of the regularsisation penalty minimising cross-validation deviance.
#' @export
#'
#' @examples

fit_gen_model <- function(gene_lengths, matrix = NULL, sample_list = NULL, gene_list = NULL, mut_types_list = NULL, col_names = NULL, table = NULL, nlambda = 100, n_folds = 10, maxit = 1e9, seed_id = 1234) {

  set.seed(seed_id)

  if (is.null(tables) & any(is.null(matrix), is.null(sample_list), is.null(gene_list), is.null(mut_types_list), is.null(col_names))) {
    stop("If not providing a full tables object, must provide the inputs sample_list, gene_list, mut_types_list and col_names")
  }

  if(!is.null(tables)) {
    matrix <- table$matrix
    sample_list <- table$sample_list
    gene_list <- table$gene_list
    mut_types_list <- table$mut_types_list
    col_names <- table$col_names
  }

  n_samples <- length(sample_list)
  n_genes <- length(gene_list)
  n_mut_types <- length(mut_types_list)

  if (any(dim(matrix) != c(n_samples, n_genes * n_mut_types))) {
    stop(paste0("Matrix has dimension ", dim(matrix), ", should have dimension ", c(n_samplesm, n_genes * n_mut_types)))
  }

  # Making output vector
  mutation_vector <- matrix
  dim(mutation_vector) <- c(n_samples * n_genes * n_mut_types, 1)

  # Constructing design matrix
  i <- 1:(n_samples*n_genes*n_mut_types)
  j <- rep(1:(n_genes*n_mut_types), each = n_samples)
  j_ <- rep(seq(1, n_genes*n_mut_types, n_mut_types), each = n_samples*n_mut_types)

  design_matrix <- Matrix::sparseMatrix(c(i, i), c(j, j_))

  t_i_getter <- Matrix::sparseMatrix(i = rep(1:n_samples, n_genes*n_mut_types), j = 1:(n_samples*n_genes*n_mut_types), dims = c(n_samples, n_samples*n_genes*n_mut_types))
  t_i <- as.vector(t_i_getter %*% mutation_vector)
  t_i_getter <- NULL

  t_s_getter <- Matrix::sparseMatrix(i = rep(1:n_mut_types, times = n_genes, each = n_samples), j = 1:(n_samples*n_genes*n_mut_types), dims = c(n_mut_types, n_samples*n_genes*n_mut_types))
  t_s <- as.vector(t_s_getter %*% mutation_vector)
  t_s_getter <- NULL

  # We use weights (rather than offsets) to fit the model, as it makes glmnet more efficient.
  weights <- rep(t_s, each = n_samples, times = n_genes)* rep(t_i, times = n_mut_types*n_genes)
  weights <- weights*rep(gene_lengths$max_cds, each = n_samples*n_mut_types)
  weights <- weights/mean(weights)

  weighted_observations <- ifelse(weights == 0, 0, as.vector(mutation_vector)/weights)

  # First (full dataset) run
  print("First glmnet run (full dataset)")
  fit <- glmnet::glmnet(x = design_matrix, y = weighted_observations, nlambda = nlambda, weights = weights, family = "poisson", trace.it = 1, maxit = maxit)

  # Cross-validation
  print("Cross-validation:")
  partitions <- sample(rep(1:n_folds, n_samples * n_genes * n_mut_types / n_folds), n_samples * n_genes * n_mut_types)
  dev <- matrix(0, n_folds, length(fit$lambda))

  for (fold in 1:n_folds) {
    part.design <- design_matrix[partitions != fold, ]
    part.response <- weighted_observations[partitions != fold]
    part.weights <- weights[partitions != fold]

    writeLines(paste("\nFitting glmnet on fold", fold))
    part.fit <- glmnet::glmnet(part.design, y = part.response, family = "poisson", lambda = fit$lambda, weights = part.weights, trace.it = 1, maxit = maxit)
    part.weights <- NULL; part.design <- NULL; part.response <- NULL;

    print("Computing statistics")
    ### This code should be made much simpler - it was originally written to squeeze the maximum amount out of a laptop without crashing.
    pb <- utils::txtProgressBar(max = 15, width = 100, style = 3)

    part.test.design <- design_matrix[partitions == fold,]
    utils::setTxtProgressBar(pb, 1)

    part.test.product <- part.test.design %*% part.fit$beta
    utils::setTxtProgressBar(pb, 2)

    part.test.design <- NULL;

    part.test.product <- part.test.product + part.fit$a0
    utils::setTxtProgressBar(pb, 3)

    part.fit <- NULL;

    part.test.product <- exp(part.test.product)
    utils::setTxtProgressBar(pb, 4)

    part.test.weights <- Matrix::Matrix(rep(weights[partitions == fold], length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda))
    utils::setTxtProgressBar(pb, 5)

    part.test.predict <- part.test.weights * part.test.product
    utils::setTxtProgressBar(pb, 6)

    part.test.product <- NULL;  part.test.weights <- NULL

    part.test.response <- Matrix::Matrix(rep(mutation_vector[1:(n_samples*n_mut_types*n_genes),][partitions == fold], length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda), sparse = TRUE)
    utils::setTxtProgressBar(pb, 7)

    part.test.residuals <- part.test.response - part.test.predict
    utils::setTxtProgressBar(pb, 8)

    part.test.logpredict <- ifelse(part.test.predict == 0, 0, log(part.test.predict))
    utils::setTxtProgressBar(pb, 9)

    part.test.predict <- NULL

    log_response <- mutation_vector[partitions == fold, ]
    log_response <- ifelse(log_response == 0, 0, log(log_response))
    utils::setTxtProgressBar(pb, 10)

    part.test.logresponse <- Matrix::Matrix(rep(log_response, length(fit$lambda)), n_samples*n_genes*n_mut_types/n_folds, length(fit$lambda), sparse = TRUE)
    utils::setTxtProgressBar(pb, 11)

    log_response <- NULL

    log_div <- part.test.logresponse - part.test.logpredict
    utils::setTxtProgressBar(pb, 12)

    part.test.logresponse <- NULL; part.test.logpredict <- NULL;

    log_div <- part.test.response*log_div
    utils::setTxtProgressBar(pb, 13)

    part.test.response <- NULL

    deviances <- log_div - part.test.residuals
    utils::setTxtProgressBar(pb, 14)

    part.test.residuals <- NULL; log_div <- NULL

    dev[fold,] <- 2*Matrix::colMeans(deviances)
    utils::setTxtProgressBar(pb, 15)

    deviances <- NULL

  }

  s_min <- max(which(colMeans(dev) == min(colMeans(dev))))
  return(list(fit = fit, dev = dev, s_min = s_min))

}
