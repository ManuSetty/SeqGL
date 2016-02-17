
#' Function to determine logistic residuals
#'
#' @param X Feature matrix
#' @param labels Labels (Should be +1/-1)
#' @param w Group lasso model
#' @return Vector of residuals

logistic.errors <- function (X, labels, w) {
	return (log (1 + exp (-labels * (X %*% w))))
}


#' Group scores for each group and example
#'
#' @param train.features Training feature matrix
#' @param labels Training labels (Should be +1/-1)
#' @param group.w Group lasso model
#' @param groups Feature assignment to groups. Each feature should belong to exactly one group.
#' @return Matrix of size no. of sequences x no. of groups.
#' @seealso \code{\link{run.group.lasso}}
#' @export

determine.peak.scores <- function (train.features, labels, group.w, groups){

	# Residuals when each group is set to 0
	peak.scores <- matrix (0, nrow (train.features), length (unique (groups)))
	for (group in unique (groups)) {	
		inds <- which (groups == group)	
		new.w <- group.w;  new.w[inds] <- 0
		peak.scores[,group] <- as.vector (logistic.errors (train.features, labels, new.w))
	}

	# Error changes
	residuals <- logistic.errors (train.features, labels, group.w)
	peak.scores <- peak.scores - as.vector (residuals)

	return (peak.scores)
}


#' Ranking of groups based on predictive power
#'
#' @param train.features Training feature matrix
#' @param labels Training labels (Should be +1/-1)
#' @param group.w Group lasso model
#' @param groups Feature assignment to groups. Each feature should belong to exactly one group.
#' @param no.cores No of cores for parallel processing
#' @return Matrix with a row for each column. Each group will have the following fields:
#' class (+1/-1), score, group.size (no. of kmers in group), class.size (no. non zero kmers in group), group
#' @seealso \code{\link{run.group.lasso}}
#' @export

determine.group.scores <- function (train.features, labels, group.w, groups,
	no.cores=1) {

	group.score <- function (group) {

		## Prediction error without features
		new.w <- group.w
		feature.index <- which (groups == group)
		new.w[feature.index]  <- 0
		peak.scores <- logistic.errors (train.features, labels, new.w) - residuals
		means <- tapply (peak.scores, labels, sum)

		## score class and value
		ind <- which (means == max (means))[1]
		ret.vector <- c(as.numeric (names (means))[ind], unname (means[ind]))
		ret.vector <- c(ret.vector, length (feature.index), 
			length (which (sign (group.w[feature.index]) == ret.vector[1])), group)
		return (ret.vector)
	}

	# Determine residuals
	time.start <- get.time ()
	residuals <- logistic.errors (train.features, labels, group.w)

	# Determine scores for all groups
	group.scores <- t(simplify2array (mclapply (1:length (unique (groups)), group.score, mc.cores=no.cores)))
	colnames (group.scores) <- c('class', 'score', 'group.size', 'class.size', 'group')

	time.end <- get.time ()
	show (sprintf ("Time for determining group scores: %.2f", (time.end - time.start)/60))

	return (group.scores)

}


#' Determine group membership
#'
#' @param labels Training labels (Should be +1/-1)
#' @param group.scores Matrix of group scores
#' @param peak.scores Matrix of example group scores
#' @param test.classes Classes to determine group members
#' @param fdr.cutoff q-value from examples of other class
#' @param no.cores No of cores for parallel processing
#' @return \item{group.members}{Binary matrix of size no. of sequences x no. of groups indicating group membership}
#'  \item{group.pvals}{FDR-corrected p-value matrix of size no. of sequences x no. of groups}
#' @seealso \code{\link{run.group.lasso}}
#' @export

determine.group.members <- function (labels, group.scores, 
	peak.scores, test.classes=c(1, -1), fdr.cutoff=0.05, no.cores=1) {

	# Set up matrix
	time.start <- get.time ()

	# Function to find group members
	find.members <- function (group) {
		time.start <- get.time ()

		# Set up group values
		group.class <- group.scores[group, 'class']
		gcs <- peak.scores[labels == group.class, group]
		if (!group.class %in% test.classes || length (which (gcs != 0)) == 0)
			return (rep (1, nrow (peak.scores)))

		# Set up empirical null
		group.emp.null <- peak.scores[labels == -group.class, group]
	    group.emp.null <- sort (group.emp.null[group.emp.null != 0])

	    # Check if no distributoin exists
	    if (length (group.emp.null) == 0){
		    pvals <- rep (1, length (gcs)); pvals[gcs > 0] <- 0
	    } else {
			# Subset necessary examples
			test.examples <- which (gcs >= 1e-2 & gcs < max (group.emp.null))
			indexes <- rep (1, length (gcs))			

			# Find indexes
			if (length (test.examples) > 0 ) {
				indexes[test.examples] <- sapply (gcs[test.examples], 
					function (x) { which (x >= group.emp.null[-length (group.emp.null)] & x < group.emp.null[-1])
					})
			}
			indexes[gcs >= max (group.emp.null)] <- length (group.emp.null)

			# Pvalues and adjust them
			pvals <- 1 - indexes / length (group.emp.null)
			pvals [gcs == 0] <- 1
	    }

		# Set up members and return
		members <- rep (NA, nrow (peak.scores))
		members[labels == group.class] <- pvals
		return (members)

	}

	# Group nominal pvalues
	group.pvals <- simplify2array (mclapply (1:ncol (peak.scores), 
		find.members, mc.cores=no.cores))

	# Adjust them and find members
	group.members <- apply (group.pvals, 2, 
		function (x) { 
			x[which (x != 1)] <- p.adjust (x[which (x != 1)], method='BH')
			members <- rep (0, length (x))
			members[which (x < fdr.cutoff)] <- 1;
			return (members)
		})


	time.end <- get.time ()
	show (sprintf ("Time for determining group members: %.2f", (time.end - time.start)/60))

	return (list (group.members=group.members, group.pvals=group.pvals))

}



