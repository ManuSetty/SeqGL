

#' Function to extract sequences and build feature matrices
#'
#' @param pos.regions Genomic Ranges object representing positive class regions
#' @param neg.regions Genomic Ranges object representing negative class regions
#' @param dictionary.file Dictionary file for the kernel build using \code{ChIPKernels}
#' @param org Name of the organism. Default is hg19. This assumes that the corresponding BSGenome library is installed
#' @param select.top.features Logical indicating whether feature selection should be performed
#' @param ... Additional arugments to \code{select.top.features}
#' @return List containging feature matrix, labels and indexes for both train and test.
#' @return Matrix of positions indicating the position in the sequence with maximum kmer score.
#' @seealso \link{\code{select.top.features}}
#' @export

build.train.test.data <- function (pos.regions, neg.regions, 
	dictionary.file, org="hg19",
	select.top.features=TRUE, ...) {

	start.time <- get.time ()
	# Sample examples for train and test
	pos.length <- length (pos.regions); neg.length <- length (neg.regions)
	pos.train.inds <- sample (pos.length)[1:ceiling ((pos.length / 2))]
	neg.train.inds <- sample (neg.length)[1:ceiling ((neg.length / 2))]
	pos.test.inds <- setdiff (1:pos.length, pos.train.inds)
	neg.test.inds <- setdiff (1:neg.length, neg.train.inds)

	# Extract sequences
    show ('Extract sequences for all training and test examples...')
    time.start <- get.time ()
    org <- load.bsgenome (org)
	seqs <- get.seqs (org, c(pos.regions, neg.regions))
	gc ()

	# Train and test sequences
    pos.seqs <- seqs[1:pos.length]; neg.seqs <- seqs[(pos.length + 1):length (seqs)]
	rm (seqs)
    gc ()

    # Check and reomove Ns in sequences
    exclude <- union (grep ('N', pos.seqs), grep ('N', neg.seqs))
    # Remove peaks which cross the chromosome limits
    exclude <- c(exclude, 
    	which (pmax (end (pos.regions), end (neg.regions)) > seqlengths (org)[as.character (seqnames (pos.regions))]))
    pos.train.inds <- setdiff (pos.train.inds, exclude); neg.train.inds <- setdiff (neg.train.inds, exclude);
    pos.test.inds <- setdiff (pos.test.inds, exclude); neg.test.inds <- setdiff (neg.test.inds, exclude);

    # Train and test seqs
    train.regions <- c(pos.regions[pos.train.inds], neg.regions[neg.train.inds])
    test.regions <- c(pos.regions[pos.test.inds], neg.regions[neg.test.inds])
	train.seqs <- c(pos.seqs[pos.train.inds], neg.seqs[neg.train.inds])
	test.seqs <- c(pos.seqs[pos.test.inds], neg.seqs[neg.test.inds])
    time.end <- get.time ()
    show (sprintf ("Time for extracting sequences: %.2f minutes", (time.end - time.start)/60 ))

	# Train and test labels
	train.labels <- c(rep (1, length (pos.train.inds)), rep (-1, length (neg.train.inds)))
	test.labels <- c(rep (1, length (pos.test.inds)), rep (-1, length (neg.test.inds)))

	# Build train and test features
	show ('Building features from dictionary...')
	features <- build.features.kernels  (dictionary.file, c(train.seqs, test.seqs), NULL, FALSE)$features
	train.features <- features[1:length (train.seqs), ]
	test.features <- features[(length (train.seqs)+1):nrow (features),]

	# Select top features
	feature.inds <- NULL
	if (select.top.features) {
		show ('Selecting top features...')
		time.start <- get.time ()
		feature.inds <- select.top.features (train.features, train.labels, ...)
		train.features <- train.features[,feature.inds]
		test.features <- test.features[,feature.inds]
	    time.end <- get.time ()
	    show (sprintf ("Time for selecting features: %.2f minutes", (time.end - time.start)/60 ))
	}

	# Package and return
	show ('Package and return...')
	results <- list (train.features = train.features, test.features = test.features,
		train.inds = list (pos=pos.train.inds, neg=neg.train.inds), 
		test.inds = list (pos=pos.test.inds, neg=neg.test.inds), feature.inds = feature.inds,
		train.labels = train.labels, test.labels = test.labels, 
		train.regions = train.regions, test.regions = test.regions,
		dictionary.file = dictionary.file, train.seqs=train.seqs, test.seqs=test.seqs)

    time.end <- get.time ()
    show (sprintf ("Total time for constructing data: %.2f minutes", (time.end - start.time)/60 ))

	return (results)

}


#' Function to select top features
#'
#' @param features Feature matrix 
#' @param labels Feature labels (Should be +1/-1)
#' @param feature.count Number of features to select. Default is 10000.
#' @param min.example.fraction Minimal fraction of examples that a feature should defined in. Default is 1\%.
#' @description The differences are determined by weighted means. The number of examples in which a particular
#' feature is defined is used to estimate the means in both classes.
#' @return Vector of selected feature indexes.
#' @export

select.top.features <- function (features, labels, 
	feature.count = 10000, min.example.fraction=0.01) {

	# Based on non zero means 
	# Determine summary of sparse matrices for faster computation
	pos.summary <- as.matrix (summary (features[labels == 1,]))
	neg.summary <- as.matrix (summary (features[labels == -1,]))

	# Differences of means
	pos.counts <- neg.counts <- pos.means <- neg.means <- rep (0, ncol (features))
	means <- tapply (pos.summary[,3], pos.summary[,2], mean)
	pos.means[as.numeric (names (means))] <- means

	means <- tapply (neg.summary[,3], neg.summary[,2], mean)
	neg.means[as.numeric (names (means))] <- means
	diffs <- abs (pos.means - neg.means)

	# Eliminate examples without significant presence
	counts <- table (pos.summary[,2])
	pos.counts[as.numeric (names (counts))] <- as.vector (counts)
	counts <- table (neg.summary[,2])
	neg.counts[as.numeric (names (counts))] <- as.vector (counts)
	inds <- which (pos.counts < min.example.fraction * length (which (labels == 1)) &
		neg.counts < min.example.fraction * length (which (labels == -1)))
	diffs[inds] <- 0

	# Sort and pick features
	feature.count <- min (feature.count, length (which (diffs > 0)))
	return (sort (diffs, index.return=TRUE, decreasing=TRUE)$ix[1:feature.count])
}
