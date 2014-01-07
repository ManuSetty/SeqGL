\name{run.seqGL.wrapper}
\alias{run.seqGL.wrapper}
\title{Wrapper for running Seqgl.}
\usage{
  run.seqGL.wrapper(pos.regions, neg.regions, org = "hg19",
    res.dir, dictionary.file, no.groups = 500,
    lambdas = c(0.01, 0.005, 0.001, 5e-04, 1e-04),
    no.cores = 1, motifs = FALSE, test.classes = c(1, -1),
    ...)
}
\arguments{
  \item{pos.regions}{GRanges object of positive regions}

  \item{neg.regions}{GRanges object of negative regions.}

  \item{org}{Organism from which data was derived}

  \item{res.dir}{Output directory}

  \item{dictionary.file}{The positional wildcard dictionary
  file built using ChIPKernel}

  \item{no.groups}{Number of groups to use for group lasso}

  \item{lambdas}{Regularization parameters for cross
  validation}

  \item{no.cores}{Number of cores for parallel processing}

  \item{motifs}{Logical indicating if motifs should be
  computed}

  \item{test.classes}{Which of the test classes to find
  motifs for}

  \item{...}{Additional arguments to bulid.train.test.data}
}
\value{
  GRanges object of enriched regions
}
\description{
  Wrapper for running Seqgl.
}
\details{
  The positive and negative regions will be split evenly to
  create training and test sets
}

