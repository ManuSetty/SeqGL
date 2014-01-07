

#' Function for finding motifs associate with each group
#'
#' @param res.dir Directory with group lasso results
#' @param dictionary.file Dictionary file for the kernel build using \code{ChIPKernels}
#' @param no.cores Number of cores for parallel processing
#' @param min.examples Minimum examples for a group to be considered
#' @param test.classes Group classes to be considered
#' @param org Organism from which data was derived
#' @details The motifs for all groups will be found under folder group_motifs. The summary of all 
#' groups and motifs will be written to group_motifs.html
#' @return GRanges object of enriched regions
#' @export

group.motifs <- function (res.dir, dictionary.file,
	no.cores=1, min.examples = 25, test.classes=c(1, -1), org='hg19') {

	# Test to check if the homer tool is in path
	test <- system ("findMotifsGenome.pl", ignore.stdout=TRUE, ignore.stderr=TRUE, intern=FALSE)
	if (test != 0) {
		stop ("Homer motif finding tool not found. Please install and add to your path.")
	}

	# Load all the group lasso data and results
    train.test.data <- readRDS (sprintf ("%s/train_test_data.Rds", res.dir))
    clustering.results <- readRDS (sprintf ("%s/clustering_results.Rds", res.dir))
	load (sprintf ("%s/group_lasso_results.Rdata", res.dir))
	

	# Create output directory
	group.motif.dir <- sprintf ("%s/group_motifs/", res.dir)
	dir.create (group.motif.dir)

	# Get test groups
	test.groups <- group.scores[colSums (group.members) > min.examples & group.scores[,'class'] %in% test.classes, 'group']

	group.motif <- function (group) {
		show (group)
		
		# Create directory for results
		homer.dir <- sprintf ("%s/Group%d", group.motif.dir, group)
		dir.create (homer.dir)

		# Determine high scoring locations
		bed.file <- sprintf ("%s/seqs.bed", homer.dir)
		examples <- which (group.members[,group] == 1)
		res <- extract.high.scoring.locations (group, group.scores[group, 'class'], 
			train.test.data$train.seqs[examples], train.test.data$train.regions[examples], 
			group.lasso.results$w[clustering.results$groups == group],
			dictionary.file, bed.file)

		# Run homer
		homer.cmd <- sprintf ("findMotifsGenome.pl %s %s %s", bed.file, org, homer.dir)
		homer.cmd <- sprintf ("%s -p 1 -size given -len 6,8,10,12 -noknown -mset vertebrates", homer.cmd)
		system (homer.cmd)
	}

	res <- mclapply (test.groups, group.motif, mc.cores=no.cores, mc.preschedule=FALSE)

	# Create report 
	lines <- c("<HTML><HEAD><TITLE>Motifs</TITLE></HEAD>",
		"<BODY>",
		"<TABLE border=\"1\" cellpading=\"0\" cellspacing=\"0\">",
		"<TR><TD>Group</TD><TD>Test Class</TD><TD>Score</TD><TD>Motif</TD><TD>Motif TF</TD></TR>")

	# Extract TF associated with motif 1 for each test group
	tfs <- c()
	no.motifs <- max (test.groups) + 1
	for (i in test.groups) {
		# Skip if no motifs were found
		if (!file.exists (sprintf ("%s/group_motifs/Group%d/homerResults.html", res.dir, i))) {
			no.motifs <- c(no.motifs, i)
			next
		}

		group.lines <- readLines (sprintf ("%s/group_motifs/Group%d/homerResults.html", res.dir, i))
		tfs <- c(tfs, gsub ('/.*', '', strsplit (group.lines[grep ('motif1', group.lines)], '<TD>')[[1]][8]))
	}
	test.groups <- setdiff (test.groups, no.motifs)

	test.class <- group.scores[test.groups, 'class']
	inds <- sort (group.scores[test.groups, 'score'], index.return=TRUE, decreasing=TRUE)$ix
	cols <- cbind ("<TR>", 
		sprintf ("<TD>&nbsp;<a href=\"group_motifs/Group%d/homerResults.html\">Group%d</a>&nbsp;</TD>", test.groups[inds], test.groups[inds]),
		sprintf ("<TD>%d</TD>", test.class[inds]),
		sprintf ("<TD>%.2f</TD>", group.scores[test.groups[inds], 'score']),
		sprintf ("<TD><IMG src=\"group_motifs/Group%d/homerResults/motif1.logo.png\" width=\"150\"/></TD>", test.groups[inds]),
		sprintf ("<TD>%s</TD>", tfs[inds]),
		"</TR>")
	lines <- c(lines, apply (cols, 1, paste, collapse=" "))
	lines <- c(lines, "</TABLE></BODY></HTML>")
	writeLines (lines, sprintf ("%s/group_motifs.html", res.dir))
    invisible (res)
}
