

#' Wrapper for running Seqgl. 
#'
#' @param pos.regions GRanges object of positive regions
#' @param neg.regions GRanges object of negative regions.  
#' @param org Organism from which data was derived
#' @param res.dir Output directory
#' @param dictionary.file The positional wildcard dictionary file built using ChIPKernel
#' @param no.groups Number of groups to use for group lasso
#' @param lambdas Regularization parameters for cross validation
#' @param no.cores Number of cores for parallel processing
#' @param motifs Logical indicating if motifs should be computed
#' @param test.classes Which of the test classes to find motifs for
#' @param fdr.cutoff FDR cutoff for identifying examples best predicted by each group
#' @param ... Additional arguments to bulid.train.test.data
#' @details The positive and negative regions will be split evenly to create training and test sets
#' @return GRanges object of enriched regions
#' @export

run.seqGL.wrapper <- function (pos.regions, neg.regions, org='hg19', res.dir,
	dictionary.file, no.groups=500, lambdas=c(1e-2, 5e-3, 1e-3, 5e-4, 1e-4),
	no.cores = 1, motifs=FALSE, test.classes=c(1, -1), fdr.cutoff=0.05, ...) {
 
	# Determine data
	train.test.data <- build.train.test.data (pos.regions, neg.regions, dictionary.file, org=org, ...)
	saveRDS (train.test.data, file=sprintf ("%s/train_test_data.Rds", res.dir))

	# Run clustering
	clustering.results <- run.clustering (train.test.data$train.features, no.groups, no.cores)
	saveRDS (clustering.results, file=sprintf ("%s/clustering_results.Rds", res.dir))

	# Group lasso parameter evaluation
	param.eval <- group.lasso.eval.parameters (train.test.data$train.features, train.test.data$train.labels,
		train.test.data$test.features, train.test.data$test.labels, clustering.results$groups, 
		lambdas=lambdas, no.cores = no.cores)
	saveRDS (param.eval, file=sprintf ("%s/param_eval.Rds", res.dir))

	# Run group lasso
    inds <- which (param.eval$aucs.matrix == max (param.eval$aucs.matrix), arr.ind=TRUE)
	group.lasso.results <- run.group.lasso (train.test.data$train.features, train.test.data$train.labels,
		train.test.data$test.features, train.test.data$test.labels, clustering.results$groups,
		param.eval$lambdas[inds[1]], param.eval$lambdas[inds[2]], no.cores=no.cores)

	# Determine group error changes and scores
	show ('Determing scores and memberes...')
	group.error.changes <- determine.group.error.changes (train.test.data$train.features, 
		train.test.data$train.labels, group.lasso.results$w, clustering.results$groups)
	group.scores <- determine.group.scores (train.test.data$train.features, 
		train.test.data$train.labels, group.lasso.results$w, clustering.results$groups, no.cores)
	group.members <- determine.group.members  (train.test.data$train.labels, 
		group.scores, group.error.changes, fdr.cutoff=fdr.cutoff,
		no.cores=no.cores, test.classes=test.classes)$group.members
	save (group.lasso.results, group.error.changes, group.scores, group.members,
		file=sprintf ("%s/group_lasso_results.Rdata", res.dir))

	# Run group lasso motifs
	if (motifs) {
		group.motifs (res.dir, dictionary.file,
			no.cores=no.cores, org=org, test.classes=1)
	}
}





#' SeqGL pipeline 
#'
#' @param peaks.file Bed file containing the peaks. See Details.
#' @param out.dir Path to the output directory
#' @param data.type ChIP or DNase: indicating the experiment though which the peaks were derived
#' @param org IUPAC organism code. hg19, hg18, mm10, mm9, mm8 are supported. Note that the corresponding 
#' BSGenome package has to be installed.
#' @param span Width of the peaks used for analysis. Default is 150
#' @param max.exmamples Maximum examples for training. Note that group membership will be determined for all examples
#' @param no.cores Number of cores for parallel processing
#' @param no.groups Number of groups to use for group lasso
#' @param dictionary.file The positional wildcard dictionary file built using ChIPKernel.
#' A dicionary will be built if it not specified.
#' @details SeqGL results will be available in \code{out.dir}
#' @export

run.seqGL <- function (peaks.file, out.dir, data.type, org,
	span=150, max.examples=ifelse (data.type == 'ChIP', 5e3, 4e4),  no.cores = 1,
	no.groups=ifelse (data.type == 'ChIP', 20, 100),
	dictionary.file=system.file( "extdata/wildcard_dict_kmer8_mismatches2_alpha5_consecutive_mis.Rdata", package="SeqGL" )) {

	start.time <- get.time ()

	# Set up constants
	if (data.type == 'ChIP') {
		feature.count <- 10000
		fdr.cutoff <- 0.1
	} 
	if (data.type == 'DNase') {
		feature.count <- 25000
		fdr.cutoff <- 0.01
	}

	show ('Determining training and test sets for SeqGL... ')
	time.start <- get.time ()
	# Read peaks, convert to GRanges and sort
	regions <- read.table (peaks.file, stringsAsFactors=FALSE, header=TRUE)
	all.regions <- GRanges (regions[,'chrom'], IRanges (regions[,'chromStart'], regions[,'chromEnd']),
		score=regions[,'score'], summit=regions[,'summit'], name=regions[,'name'])

	# Error checks
	if (any (is.na (all.regions$name)))
	   stop ("Please specify peak names for all peaks using the 'name' column in peaks file...")
	if (any (is.na (all.regions$summit)))
	   stop ("Please specify summit positions for all peaks using the 'summit' column in peaks file...")
	if (any (is.na (all.regions$score)))
	   stop ("Please specify scores for all peaks using the 'summit' column in peaks file...")

	# Center peaks and find negatives
	all.regions <- all.regions[sort (all.regions$score, index.return=TRUE, decreasing=TRUE)$ix]
	start (all.regions) <- end (all.regions) <- start (all.regions) + all.regions$summit - 1
	all.regions <- resize (all.regions, fix='center', span)

	# Positive and negative examples
	pos.regions <- all.regions[1:min (max.examples, length (all.regions))]
	neg.regions <- shift (pos.regions, span * 2)

	# Remove overlaps
	use.inds <- which (countOverlaps (neg.regions, pos.regions)  == 0)
	pos.regions <- pos.regions[use.inds]; neg.regions <- neg.regions[use.inds]
	time.end <- get.time ()
	show (sprintf ("Time for building training and test sets: %.2f", (time.end - time.start) / 60))

	# Build dictionary
	show ('Building dicitionary used for determining features...')
	if (is.null (dictionary.file)) {
		build.wildcard.dictionary (8, 2, out.dir, c('A', 'C', 'G', 'T', 'N'), TRUE)
		dictionary.file <- sprintf ("%s/wildcard_dict_kmer8_mismatches2_alpha5_consecutive_mis.Rdata", out.dir)
	}

	# Run group lasso
	show ('Running SeqGL...')
	run.seqGL.wrapper (pos.regions, neg.regions, org, out.dir,
		dictionary.file, no.groups, no.cores = no.cores, fdr.cutoff=fdr.cutoff,
		motifs=TRUE, test.classes=1, feature.count=feature.count)
	show ('SeqGL complete')

	# Find group membership for all examples
	show ('Determining group membership for all examples...')
	time.start <- get.time ()

    # Load everything
    train.test.data <- readRDS (sprintf ("%s/train_test_data.Rds", out.dir))
    clustering.results <- readRDS (sprintf ("%s/clustering_results.Rds", out.dir))
	load (sprintf ("%s/group_lasso_results.Rdata", out.dir))
	# Plot  auc
	get.auc (train.test.data$test.labels, group.lasso.results$test.preds,
		sprintf ("%s/test_auc.pdf", out.dir))

    # Sequences and labels
	all.negs <- shift (all.regions, span*2)
	all.features <- build.features.kernels (dictionary.file, get.seqs (load.bsgenome (org),
                                  c(all.regions, all.negs)), kmers=colnames (train.test.data$train.features), verbose=FALSE)$features
	all.labels <- rep (c(1, -1), each=length (all.regions))

	# Group members setups etc
	group.error.changes <- determine.group.error.changes (all.features, all.labels, 
		group.lasso.results$w, clustering.results$groups)
	group.scores <- determine.group.scores (all.features, all.labels,
		group.lasso.results$w, clustering.results$groups, no.cores)
	group.members <- determine.group.members  (all.labels,
		group.scores, group.error.changes, fdr.cutoff=fdr.cutoff,
		no.cores=no.cores, test.classes=1)$group.members
	colnames (group.members) <- sprintf ("Group%d", 1:ncol (group.members))

	# Write to files
	group.dir <- sprintf ("%s/group_members", out.dir)
	dir.create (group.dir)
	test.groups <- dir (sprintf ("%s/group_motifs", out.dir))
	for (group in test.groups) 
			granges.to.bed (all.regions[which (group.members[,group] == 1)], sprintf ("%s/%s_members.bed", group.dir, group))
	time.end <- get.time ()
	show (sprintf ("Time for determining group membership: %.2f", (time.end - time.start) / 60))

	end.time <- get.time ()
	show (sprintf ("Processing complete. Total time: %.2f", (end.time - start.time) / 60))
	show (sprintf ("SeqGL results are in %s/group_motifs.html", out.dir))
	show (sprintf ("Group membership information is in %s", group.dir))

}












