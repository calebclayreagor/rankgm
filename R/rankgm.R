#' Rank gene-expression modules
#'
#' @description
#' Infers gene expression-module rankings by fitting negative binomial mixture models to gene-count data in single cells.
#' Gene-module rankings are then used to uniquely rank-order genes for each cell in the dataset.
#'
#' @param object Loom file, provided as an H5File-class object, with read/write access. Must contain gene-count data
#'   in the primary `matrix`. See `LoomPy` and `hdf5r` for more details.
#' @param chunk.size Specifies the maximum number of single cells to load into memory at once. Default is 500.
#' @param verbose Specifies whether to print update messages during the fitting of mixture models. Default is TRUE.
#' @param parallel Specifies whether to parallelize the fitting of mixture models for each batch. Default is FALSE.
#' @param overwrite Specifies whether to overwrite any existing results if they are found. Default is FALSE.
#' @param BPPARAM Specifies the back-end to use for parallelization. See `bpparam` in `BiocParallel` for more details.
#'
#' @return Returns an H5File-class object with an updated Loom file containing a `rankgm` layer with the inferred single-
#'   cell-specific gene rankings, and a column attribute `rankgm_success` indicating the success (1) or failure (0) of the
#'   model-fitting procedure for each cell. The ranks are returned as all zeros if the fitting procedure is unsuccessful.
#'
#' @examples
#' \dontrun{
#' pbmc.h5 <- H5File$new("pbmc3k.loom", mode="r+")
#' pbmc.h5 <- RankGeneModules(pbmc.h5, overwrite = TRUE)
#' pbmc.h5$close_all()
#' }
#'
#' @export RankGeneModules
RankGeneModules <- function(object, chunk.size = 500, verbose = TRUE, parallel = FALSE,
                            overwrite = FALSE, BPPARAM = BiocParallel::bpparam()) {

  # TO-DO:
  #   -checks for chunk.size, verbose, parallel, overwrite, BPPARAM

  # 'object' checks: H5File-class with r+ access, contains required groups/datasets, overwrite (y/n)
  if (class(object)[1] != 'H5File') {
    stop("Object must be an H5File-class object. See package hdf5r for more details.")}
  if (object$mode != 'r+') {
    stop("H5File must have read/write access. Use with mode 'r+'. See package hdf5r for more details.")}
  if (!("matrix" %in% names(object))) {
    stop("Loom file must contain a 'matrix' dataset. See LoomPy for more details.")}
  if (!("layers" %in% names(object))) {
    stop("Loom file must contain a group 'layers'. See LoomPy for more details.")}
  if (!("col_attrs" %in% names(object))) {
    stop("Loom file must contain a group 'col_attrs'. See LoomPy for more details.")}
  if (("rankgm" %in% names(object[["layers"]])) & (!overwrite)) {
    stop("Loom file already contains a 'rankgm' layer. To overwrite, set flag to TRUE.")}
  if (("rankgm_success" %in% names(object[["col_attrs"]])) & (!overwrite)) {
    stop("Loom file already contains a 'rankgm_success' attribute. To overwrite, set flag to TRUE.")}

  # register BP back-end for parallelization
  if (parallel) BiocParallel::register(BPPARAM)

  # select groups: layers, col_attrs
  layers <- object[["layers"]]
  col_attrs <- object[["col_attrs"]]

  # get dimensions of 'matrix'
  dims <- object[["matrix"]]$dims

  # create dataset 'rankgm' in group 'layers': unique gene rankings
  if ("rankgm" %in% names(layers)) layers$link_delete("rankgm")
  sp_rankgm <- hdf5r::H5S$new(dims = dims, maxdims = c(Inf, dims[2]))
  layers$create_dataset(name = "rankgm", space = sp_rankgm,
                        dtype = hdf5r::h5types$H5T_NATIVE_INT)

  # create dataset 'rankgm_success' in group 'col_attrs': model fitting success (0/1)
  if ("rankgm_success" %in% names(col_attrs)) col_attrs$link_delete("rankgm_success")
  sp_rankgm_success <- hdf5r::H5S$new(dims = dims[1], maxdims = Inf)
  col_attrs$create_dataset(name = "rankgm_success", space = sp_rankgm_success,
                           dtype = hdf5r::h5types$H5T_NATIVE_INT)

  # row (cell) index boundaries for chunking, used like (lower, upper]
  chunk.indices <- c(seq(from = 0, to = dims[1], by = chunk.size), dims[1])

  # iterate over data chunks (n_bounds-1)
  n.iter <- (length(chunk.indices) - 1)
  if (verbose) message("Fitting models to find the unique gene rankings for each cell...")
  success.total <- 0
  for (i in seq(n.iter)) {

    # load data chunk of specified size and convert to sparse counts matrix
    chunk.data <- object[["matrix"]][(chunk.indices[i]+1):chunk.indices[i+1],]
    n.cells <- dim(chunk.data)[1]
    counts.dgCMat <- Matrix::t(as(chunk.data, 'dgCMatrix'))

    # check that all non-zero elements in chunk.data are positive integers
    if (any(counts.dgCMat@x < 0)) {
      stop("Negative values found in 'matrix'. All gene counts must be positive integers.")}
    if (!all(counts.dgCMat@x == floor(counts.dgCMat@x))) {
      stop("Non-integer values found in 'matrix'. All counts must be positive integers.")}

    # convert non-zero counts to nested list for (possibly parallel) lapply
    intervals <- findInterval(seq(Matrix::nnzero(counts.dgCMat)), counts.dgCMat@p, left.open = TRUE)
    counts.list <- split(counts.dgCMat@x, intervals)

    # (possibly parallelized) inference of gene-expression modules and gene rankings from gene counts
    if (parallel) {
      results <- BiocParallel::bplapply(counts.list, inferRanks, n.genes.total = nrow(counts.dgCMat))
    } else {
      results <- lapply(counts.list, inferRanks, n.genes.total = nrow(counts.dgCMat))}

    # separate gene rankings from binary success indicator (0/1)
    ranks.list <- lapply(results, function(x) utils::head(x, -1))
    ranks <- unsplit(ranks.list, intervals)
    fits.success <- unlist(lapply(results, function(x) utils::tail(x, 1)))
    success.n <- sum(fits.success)
    success.total <- success.total + success.n

    # convert back to matrix (sparse -> dense), then set gene ranks for cells in chunk
    ranks.dgCMat <- Matrix::sparseMatrix(i = (counts.dgCMat@i+1), p = counts.dgCMat@p,
                                         x = ranks, dims = counts.dgCMat@Dim)
    ranks.Mat <- Matrix::t(as(Matrix::drop0(ranks.dgCMat), 'Matrix'))
    #ranks.dgTMat <- Matrix::t(as(Matrix::drop0(ranks.dgCMat), 'TsparseMatrix'))
    #chunk.ranks <- layers[["rankgm"]][(chunk.indices[i]+1):chunk.indices[i+1],]
    #chunk.ranks[cbind((ranks.dgTMat@i+1),(ranks.dgTMat@j+1))] <- ranks.dgTMat@x
    layers[["rankgm"]][(chunk.indices[i]+1):chunk.indices[i+1],] <- as(ranks.Mat, "integer")

    # set binary success indicator column attribute ('rankgm_success') for cells in chunk
    col_attrs[["rankgm_success"]][(chunk.indices[i]+1):chunk.indices[i+1]] <- as.integer(fits.success)

    if (verbose) {
      message(paste("Finished batch", i, "of", n.iter, "- successful fits for", success.n, "of", n.cells, "cells"))
    }
  }

  if (verbose) message(paste("Done. Successful fits for", success.total, "of", dims[1], "total cells."))

  # return H5File
  return(object)}


#' Infer gene expression-module rankings
#'
#' @description
#' Fits a series of negative binomial mixture models to the gene counts from a single cell to first infer and rank the
#' underlying gene-expression modules, then rank-order genes according to their most likely module.
#'
#' @param y.nz List of all non-zero gene counts from a given cell
#' @param n.genes.total Total number of genes in the given cell, including zeros
#' @param thresh.p The p-value threshold for adding new distributions to the mixture model, see `lkrTest` for details
#' @param max.n The maximum number of distributions to include in the final mixture model
#' @param seed Fixed seed for reproducibility, default is 1234
#'
#' @return Returns a single numeric vector containing the single-cell-specific gene rankings for the non-zero gene
#'   counts, followed by a binary indicator of the success (1) or failure (0) of the model-fitting procedure. The
#'   ranks are returned as all zeros if the fitting procedure is unsuccessful.
#'
#' @keywords internal
#' @export
inferRanks <- function(y.nz, n.genes.total, thresh.p = 0.05, max.n = 10, seed = 1234) {
  set.seed(seed)

  # append omitted zeros to the vector of gene counts
  y.nz <- as.numeric(unlist(y.nz))
  y <- c(y.nz, rep(0, (n.genes.total - length(y.nz))))

  # fit mixture models to infer the gene-expression modules
  for (i in seq(2, max.n)) {

    # set the initial values for the model parameters
    dim.theta <- ((3*i) - 1)
    theta <- rep(0., dim.theta)
    widx <- seq(1, dim.theta-2, 3)
    nidx <- c(seq(2, dim.theta-2, 3), dim.theta-1)
    pidx <- c(seq(3, dim.theta-2, 3), dim.theta)
    theta[widx] <- (9*(10^(-seq(i-1))))
    theta[nidx] <- seq(-3, 1, 4/(i-1))
    theta[pidx] <- seq(.9, .1, -.8/(i-1))

    # ui matrix, for constrained optimization
    n.constr <- ((6*i) - 3)
    n.lower <- ((2*i) - 1)
    n.upper <- (i + 1)
    n.relate <- (3 * (i-1))
    ui <- matrix(0, n.constr, dim.theta)

    # lower bounds
    ii <- seq(n.lower)
    jj <- c(widx, pidx)
    ui[cbind(ii,jj)] <- 1

    # upper bounds
    ii <- n.lower + c(seq(n.upper-1), rep(n.upper, length(widx)))
    jj <- c(pidx, widx)
    ui[cbind(ii,jj)] <- -1

    # p_k > p_k+1
    ii <- rep(n.lower + n.upper + seq(n.relate/3), each = 2)
    if (i==2) { jj <- pidx
    } else { jj <- c(pidx[1], rep(pidx[2:(length(pidx)-1)], each=2), pidx[length(pidx)]) }
    ui[cbind(ii,jj)] <- rep(c(1,-1), n.relate/3)

    # n_k < n_k+1
    ii <- rep(n.lower + n.upper + n.relate/3 + seq(n.relate/3), each = 2)
    if (i==2) { jj <- nidx
    } else { jj <- c(nidx[1], rep(nidx[2:(length(nidx)-1)], each=2), nidx[length(nidx)]) }
    ui[cbind(ii,jj)] <- rep(c(-1,1), n.relate/3)

    # w_k > w_k+1
    ii <- (n.lower + n.upper + 2*n.relate/3 + seq(n.relate/3))
    if (i==3) { ii <- rep(ii, each = 2)
    } else if (i>3) { ii <- c(rep(ii[1:(length(widx)-1)], each = 2), rep(ii[length(widx)], length(widx))) }
    if (i==2) { jj <- widx
    } else if (i==3) { jj <- c(widx, widx)
    } else { jj <- c(widx[1], rep(widx[2:(length(widx)-1)], each=2), widx[length(widx)], widx) }
    ui[cbind(ii,jj)] <- c(rep(c(1,-1), (n.relate/3) - 1), c(rep(1, length(widx)-1), 2))

    # ci values, for constrained optimization
    ci <- rep(-1, n.constr)
    ii <- c(seq(n.lower), n.lower + n.upper + seq(n.relate-1))
    ci[ii] <- sqrt(.Machine$double.eps)
    ci[n.lower + n.upper] <- (-1 + sqrt(.Machine$double.eps))
    ci[n.constr] <- (1 + sqrt(.Machine$double.eps))

    # constrained optimization to infer the gene-expression modules
    p.val <- 0; b <- rep(1, 2)
    if (i==2) {
      f1 <- stats::constrOptim(theta, negLL, grLL, ui, ci, y=y)
    } else {
      f2 <- stats::constrOptim(theta, negLL, grLL, ui, ci, y=y)
      p.val <- lkrTest(f1, f2) # likelihood ratio test
      y.rank <- dnbmix(f2$par, y, FALSE, FALSE) # equal priors
      y.rank <- apply(y.rank, 2, which.max)
      b <- matrixStats::binCounts(y.rank, bx=seq(max(y.rank)+1)) }
    if ((p.val>thresh.p)|(sum(b>0)<i)) break
    else if (i>2) f1 <- f2 }

  # rank-order genes by assigning them to their most likely gene-expression modules
  y.rank <- dnbmix(f1$par, y, FALSE, FALSE) # equal priors
  y.rank <- apply(y.rank, 2, which.max)

  # model fitting is successful if optimization converges and >=1 gene in each module
  b <- matrixStats::binCounts(y.rank, bx=seq(max(y.rank)+1))
  fit.success <- as.numeric(f1$convergence == 0) * as.numeric(sum(b>0) == length(b))

  # gene rankings, for non-zero counts only
  if (fit.success) {
    y.rank <- (y.rank-1)
    y.rank <- y.rank[seq(length(y.nz))]
  } else {
    # ranks are all zero if unsuccessful
    y.rank <- rep(0, length(y.nz)) }

  # return rankings, success (0/1)
  return(c(y.rank, fit.success)) }


#' Negative binomial mixture probability density
#'
#' @description
#' Probability density for the mixture of \eqn{k} negative binomial distributions: \eqn{w_1 NB(n_1, p_1) + w_2
#' NB(n_2, p_2) + ... + (1 - \Sigma w) NB(n_k, p_k)}
#'
#' @param theta Numeric vector of parameter values \eqn{w_1, n_1, p_1, w_2, n_2, p_2, ..., p_{k-1}, n_k, p_k}
#' @param y Numeric vector of gene counts from a single cell, including zeros
#' @param dsum Logical, whether to sum the probabilities across the distributions, default is TRUE
#' @param weights Logical, whether to weight the probabilities according to the parameters \eqn{w_1, w_2, ...,
#'   w_{k-1}}, default is TRUE
#'
#' @return Returns a matrix of probability densities for the given gene counts, possibly weighted and/or summed
#'   across the mixture of distributions.
#'
#' @keywords internal
#' @export
dnbmix <- function(theta, y, dsum = TRUE, weights = TRUE) {

  # theta defines the number of distributions in the mixture
  n.mix <- ((length(theta) + 1) / 3)
  dmix <- matrix(nrow = n.mix, ncol = length(y))

  # select the parameters for dist_i (w_i, n_i, p_i)
  for (i in seq(n.mix)) {

    ## parameters if only one distribution in the mixture
    #if (n.mix == 1) {
    #  w_i <- 1.; n_i <- (10^theta[1])
    #  p_i <- theta[2]

    # parameters for distributions i=1 to k-1
    if (i < n.mix) {
      w_i <- theta[((i-1)*3)+1]
      n_i <- (10^theta[((i-1)*3)+2])
      p_i <- theta[((i-1)*3)+3]

    # parameters for the final distribution, k
    } else {
      widx <- seq(1, length(theta)-2, 3)
      w_i <- (1 - sum(theta[widx]))
      n_i <- (10^theta[length(theta)-1])
      p_i <- theta[length(theta)] }

    # find the (possibly weighted) probability densities
    dmix[i,] <- stats::dnbinom(y, n_i, p_i)
    if (weights == TRUE) dmix[i,] <- w_i * dmix[i,]
    }

  # sum across distributions, if specified
  if (dsum == TRUE) dmix <- colSums(dmix)

  # return densities
  return(dmix) }


#' Negative log-likelihood
#'
#' @description
#' Computes \eqn{-\ell(\theta|y)} for the gene counts \eqn{y} and parameters \eqn{\theta}
#'
#' @param theta Numeric vector of parameter values, see `dnbmix` for details
#' @param y Numeric vector of gene counts from a single cell, including zeros
#'
#' @return Returns the negative log-likelihood value for \eqn{\theta|y}
#'
#' @keywords internal
#' @export
negLL <- function(theta, y) {
  return(-sum(log(dnbmix(theta, y))))}


#' Log-likelihood gradient
#'
#' @description
#' Computes \eqn{\frac{\partial \ell}{\partial \theta}} for the gene counts \eqn{y} and parameters \eqn{\theta}
#'
#' @param theta Numeric vector of parameter values, see `dnbmix` for details
#' @param y Numeric vector of gene counts from a single cell, including zeros
#' @param clip.scale Maximum allowed magnitude of the gradient, relative to the current parameter value
#'
#' @return Returns clipped values of the log-likelihood gradient: \eqn{\nabla \ell = \frac{\partial \ell}
#'   {\partial w_1}, \frac{\partial \ell}{\partial n_1}, \frac{\partial \ell}{\partial p_1}, \frac{\partial
#'   \ell}{\partial w_2}, \frac{\partial \ell}{\partial n_2}, \frac{\partial \ell}{\partial p_2}, ...,
#'   \frac{\partial \ell}{\partial p_{k-1}}, \frac{\partial \ell}{\partial n_k}, \frac{\partial \ell}
#'   {\partial p_k}}
#'
#' @keywords internal
#' @export
grLL <- function(theta, y, clip.scale = 2e4) {
  N <- length(y)

  # theta defines the number of distributions in the mixture
  dim.theta <- length(theta)
  n.mix <- ((dim.theta+1) / 3)
  widx <- seq(1, dim.theta-2, 3)
  nidx <- c(seq(2, dim.theta-2, 3), dim.theta-1)
  theta[nidx] <- (10^theta[nidx])

  # compute the gradient values
  for (i in seq(1, n.mix)) {
    if (i < n.mix) {
      # for distributions i=1 to k-1: return the gradient values for w_i, n_i, and p_i
      gr_i <- c((N/theta[(3*(i-1))+1]) - (N/(1 - sum(theta[widx]))),
                sum(digamma(y + theta[(3*(i-1))+2])) - N*(digamma(theta[(3*(i-1))+2]) + log(theta[(3*(i-1))+3])),
                -sum(y/(1 - theta[(3*(i-1))+3])) + (N*theta[(3*(i-1))+2]/theta[(3*(i-1))+3]))
    } else {
      # for the last distribution: return the gradient values for n_k and p_k (no w_k)
      gr_i <- c(sum(digamma(y + theta[length(theta)-1])) - N*(digamma(theta[length(theta)-1]) + log(theta[length(theta)])),
                -sum(y/(1 - theta[length(theta)])) + (N*theta[length(theta)-1]/theta[length(theta)])) }
    if (i == 1) gr <- gr_i else gr <- c(gr, gr_i) }

  # return the clipped gradient values
  return(pmin(clip.scale * theta, abs(gr)) * sign(gr)) }


#' Likelihood-ratio test
#'
#' @description
#' Uses a chi-square test to compute the probability of rejecting the null hypothesis that the simpler ('f1') and more
#' complex ('f2') models both fit the data equally well, given the ratio of the models' likelihoods.
#'
#' @param f1 Fitted model, see `optim` in `stats` for more details
#' @param f2 See f1
#'
#' @return Returns the probability of rejecting the null hypothesis from a one-tailed chi-squared distribution
#'
#' @keywords internal
#' @export
lkrTest <- function(f1, f2) {
  lam <- (2 * (f1$value - f2$value)) # ratio
  df <- (length(f2$par) - length(f1$par)) # d.o.f.
  return((1 - stats::pchisq(lam, df))) }
