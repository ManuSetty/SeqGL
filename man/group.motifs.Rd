\name{group.motifs}
\alias{group.motifs}
\title{Function for finding motifs associate with each group}
\usage{
  group.motifs(res.dir, dictionary.file, no.cores = 1,
    min.examples = 25, test.classes = c(1, -1),
    org = "hg19")
}
\arguments{
  \item{res.dir}{Directory with group lasso results}

  \item{dictionary.file}{Dictionary file for the kernel
  build using \code{ChIPKernels}}

  \item{no.cores}{Number of cores for parallel processing}

  \item{min.examples}{Minimum examples for a group to be
  considered}

  \item{test.classes}{Group classes to be considered}

  \item{org}{Organism from which data was derived}
}
\value{
  GRanges object of enriched regions
}
\description{
  Function for finding motifs associate with each group
}
\details{
  The motifs for all groups will be found under folder
  group_motifs. The summary of all groups and motifs will
  be written to group_motifs.html
}

