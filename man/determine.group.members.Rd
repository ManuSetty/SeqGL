\name{determine.group.members}
\alias{determine.group.members}
\title{Determine group membership}
\usage{
  determine.group.members(labels, group.scores,
    group.error.changes, test.classes = c(1, -1),
    fdr.cutoff = 0.05, no.cores = 1)
}
\arguments{
  \item{labels}{Training labels (Should be +1/-1)}

  \item{group.scores}{Matrix of group scores}

  \item{group.error.changes}{Matrix of example group
  scores}

  \item{test.classes}{Classes to determine group members}

  \item{fdr.cutoff}{q-value from examples of other class}

  \item{no.cores}{No of cores for parallel processing}
}
\value{
  \item{group.members}{Binary matrix of size no. of
  sequences x no. of groups indicating group membership}
  \item{group.pvals}{FDR-corrected p-value matrix of size
  no. of sequences x no. of groups}
}
\description{
  Determine group membership
}
\seealso{
  \link{\code{run.group.lasso}}
}

