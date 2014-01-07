\name{determine.group.error.changes}
\alias{determine.group.error.changes}
\title{Group scores for each group and example}
\usage{
  determine.group.error.changes(train.features, labels,
    group.w, groups)
}
\arguments{
  \item{train.features}{Training feature matrix}

  \item{labels}{Training labels (Should be +1/-1)}

  \item{group.w}{Group lasso model}

  \item{groups}{Feature assignment to groups. Each feature
  should belong to exactly one group.}
}
\value{
  Matrix of size no. of sequences x no. of groups.
}
\description{
  Group scores for each group and example
}
\seealso{
  \link{\code{run.group.lasso}}
}

