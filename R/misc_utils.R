
# Current time in seconds
get.time <- function () 
	as.numeric(format(Sys.time(), "%s"))



# Function to convert GRanges to bed
granges.to.bed <- function (regions, file) {

	# Build data frame
	df <- data.frame (chr=as.character (seqnames (regions)),
		start = start (regions),
		end = end (regions), name=as.character (regions$name))

	# Write to file
	write.table (df, file, row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
}
