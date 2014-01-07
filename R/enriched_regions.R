
#' Function to extract high scoring regions for a particular group
#'
#' @param test.group Group to be tested
#' @param test.class Class of the group being tested
#' @param seqs Sequences of examples belonging to group
#' @param regions Genomic location of examples belonging to group
#' @param group.model Logistic regression model for the particular group
#' @param dictionary.file Dictionary file for the kernel build using \code{ChIPKernels}
#' @param bed.file Bed file for writing extracted regions
#' @param half.binding.width Half of the enriched region size
#' @return GRanges object of enriched regions

extract.high.scoring.locations <- function (test.group, test.class, seqs, regions, group.model,
	dictionary.file, bed.file=NULL, half.binding.width = 25) {

	# Set up all variables
	show (sprintf ("Set up all data (%d)...", test.group))
	load (dictionary.file)

	# Sequences and features
	kmers <- names (group.model)
	kmer.features <- sapply (strsplit (kmers, '.', fixed=TRUE), function (x) { x[length (x)]})

	# Lengths
	kmer.len <- nchar (colnames (pairwise.kmers)[1])
	seq.len <- width (seqs[1])

	# Determine weighted scores
	show (sprintf ("Determine weighted scores (%d)...", test.group))
	weighted.scores <- sparseMatrix (length (seqs), seq.len, x=0)

	for (i in 1:(seq.len - kmer.len + 1)) {
		sub <- as.character (substring (seqs, i, i + kmer.len -1))
		weighted.scores[,i] <- ChIPKernels:::kmer.pairwise.features (pairwise.kmers, sub, 
			kmer.column.mapping[kmer.features]) %*% group.model[kmers]
	}
	# Adjust for class
	weighted.scores <- weighted.scores * test.class

	# Extract local scoring regions
	# max.postions closest to mid points
	show (sprintf ('Determine positions and examples (%d)...', test.group))
	mid.point <- ncol (weighted.scores) /2 
	max.positions <- apply (weighted.scores, 1, function (x) { which ( x == max (x))}) 
	max.positions <- sapply (max.positions, function (x) {y <- abs (x - mid.point); x[y==min(y)][1] })
	# Spans around mid point
	spans <- cbind (pmax (1, max.positions - half.binding.width),
		pmin (seq.len, max.positions + half.binding.width))

	# Granges to bed for Homer
	p <- regions
	end (p) <- start (p) + spans[,2]  - 1
	start (p) <- start (p) + spans[,1] 
	
	# Output bed and fa
	df <- cbind (as.character (seqnames (p)), start (p), end (p), sprintf ("Region%d", 1:length (p)))
	if (!is.null (bed.file))
		write.table (df, file=bed.file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	
	invisible (p)
}



