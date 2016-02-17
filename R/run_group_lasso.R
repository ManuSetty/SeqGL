

#' Hierarchical clustering of features
#'
#' @param features Feature matrix for clustering. Rows are sequences, columns are kmers
#' @param no.groups Number of desired groups
#' @param no.cores Number of cores for parallel processing
#' @description This function performs hierarchical clustering of features based on their correlation
#' @return List containing \code{groups} which show the feature assignments and \code{linkage}
#' @seealso \code{\link{build.train.test.data}}
#' @export

run.clustering <- function (features, no.groups=200, no.cores=1) {

	time.start <- get.time ()
	show ('Running clustering...')

	# Determine correlation distance
	if (no.cores > 1)
	    allowWGCNAThreads (no.cores)
	dist <- as.dist (1 - WGCNA::cor (as.matrix (features)))

	# Hierarchical clustering
	clustering <- hclust (dist, method='complete')

	# Groups
	groups <- cutree (clustering, k= no.groups)
    time.end <- get.time ()
    show (sprintf ("Time for running clustering: %.2f minutes", (time.end - time.start)/60 ))

	return (list (linkage=clustering, groups=groups)) 
}



#' Function to run group lasso with the specified parameters
#'
#' @param features Feature matrix
#' @param labels Example labels (Should be +1/-1)
#' @param groups Feature assignment to groups. Each feature should belong to exactly one group
#' @param lambda1, lambda2 Regularization parameters
#' @param no.cores Number of cores for parallel processing
#' @description Using the spams toolbox to run group lasso
#' @return The w vector, weights of each feature
#' @seealso \code{\link{run.group.lasso}}
#' @export

group.lasso <- function (features, labels, groups, lambda1, lambda2, no.cores=1) {

	# Run group lasso with appropriate commands
	res <- spams.fistaFlat(as.matrix (labels), features, as.matrix (rep (0, ncol (features))),
		TRUE, numThreads = no.cores, verbose = TRUE, lambda1= lambda1, lambda2=lambda2, 
		loss = 'logistic',regul = 'sparse-group-lasso-l2', tol=1e-4, max_it=1000, 
		groups=matrix (as.integer (groups, nrow=1)))

	# Extract and return w vector
	wvec <- as.vector (res[[1]])
	return (wvec)
}


#' Group lasso cross validation
#'
#' @param train.features Training feature matrix
#' @param train.labels Training labels (Should be +1/-1)
#' @param test.features Test feature matrix
#' @param test.labels Test labels (Should be +1/-1)
#' @param groups Feature assignment to groups. Each feature should belong to exactly one group
#' @param lambdas Vector of regularization parameters
#' @param no.cores Number of cores for parallel processing
#' @description Using the spams toolbox to run group lasso at various parameters to identify 
#' parameter combination with best test auc.
#' @return List containing \code{lambdas} and a matrix of test aucs named \code{auc.matrix}
#' @seealso \code{\link{run.group.lasso}}
#' @export

group.lasso.eval.parameters <- function (train.features, train.labels,
	test.features, test.labels, groups,
	lambdas = c(1e-2, 5e-3, 1e-3, 5e-4, 1e-4 ), no.cores=1) {

	time.start <- get.time ()
	show ('Running crossvalidation ...')

		
	# Run group lasso for different combination of lambdas
	auc.matrix <- matrix (0, length (lambdas), length (lambdas))
	for (i in 1:length (lambdas)) {
		for (j in 1:length (lambdas)) {
			show (sprintf ("i: %d of %d, j: %d of %d", i, length (lambdas), j, length (lambdas)))

			# Run group lasso and determine test auc
			w <- group.lasso (train.features, train.labels, groups, 
				lambdas[i], lambdas[j], no.cores)
			auc.matrix[i,j] = get.auc (test.labels, as.vector (test.features %*% w))$auc

		}
	}
	gc ()
    time.end <- get.time ()
    show (sprintf ("Time for running crossvalidation: %.2f minutes", (time.end - time.start)/60 ))

	return (list (lambdas=lambdas, aucs.matrix=auc.matrix))
}



#' Run group lasso for peaks
#'
#' @param train.features Training feature matrix
#' @param train.labels Training labels (Should be +1/-1)
#' @param test.features Test feature matrix
#' @param test.labels Test labels (Should be +1/-1)
#' @param groups Feature assignment to groups. Each feature should belong to exactly one cluster.
#' @param lambad1 L1 regularization parameter
#' @param lambda2 L2 regularization parameter
#' @param no.cores Number of cores for parallel processing
#' @description Using the spams toolbox to run group lasso at regularization parameters selected 
#' using \code{\link{group.lasso.eval.parameters}}
#' @return List of results containing \code{w}, \code{test.preds}, \code{aucs}
#' @seealso \code{\link{run.group.lasso}}
#' @export


run.group.lasso <- function (train.features, train.labels,
	test.features, test.labels, groups, lambda1, lambda2, no.cores=1) {

	time.start <- get.time ()
	show ('Running group lasso...')

	# Run group lasso and determine test auc
	results <- list ()
	results$w <- group.lasso (train.features, train.labels, groups, 
		lambda1, lambda2, no.cores)
	results$test.preds <- as.vector (test.features %*% results$w)
	results$test.auc <- get.auc (test.labels, results$test.preds)
	names (results$w) <- colnames (train.features)
	
	# Save parameters
	results$params <- list ()
	results$params$lambda1 <- lambda1
	results$params$lambda2 <- lambda2
	gc ()
    time.end <- get.time ()
    show (sprintf ("Time for running group lasso: %.2f minutes", (time.end - time.start)/60 ))

	return (results)
}


