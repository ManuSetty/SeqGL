\name{get.seqs}
\alias{get.seqs}
\title{Get sequences for genomic reanges}
\usage{
  get.seqs(org, regions, no.cores = 1)
}
\arguments{
  \item{org}{BSGenome object of the organism}

  \item{regions}{GenomicRanges object of the co ordinates}

  \item{no.cores}{Number of threads for parallel
  processing. Default 3}

  \item{metadat.field}{Metadata field name for naming
  sequences}

  \item{fasta.file}{Name of the fasta.file in which
  sequences need to written}
}
\value{
  DNAStringSet object of sequences
}
\description{
  Get sequences for genomic reanges
}

