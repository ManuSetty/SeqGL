\name{run.group.lasso}
\alias{run.group.lasso}
\title{Run group lasso for peaks}
\usage{
  run.group.lasso(train.features, train.labels,
    test.features, test.labels, groups, lambda1, lambda2,
    no.cores = 1)
}
\arguments{
  \item{train.features}{Training feature matrix}

  \item{train.labels}{Training labels (Should be +1/-1)}

  \item{test.features}{Test feature matrix}

  \item{test.labels}{Test labels (Should be +1/-1)}

  \item{groups}{Feature assignment to groups. Each feature
  should belong to exactly one cluster.}

  \item{lambad1}{L1 regularization parameter}

  \item{lambda2}{L2 regularization parameter}

  \item{no.cores}{Number of cores for parallel processing}
}
\value{
  List of results containing \code{w}, \code{test.preds},
  \code{aucs}
}
\description{
  Using the spams toolbox to run group lasso at
  regularization parameters selected using
  \link{\code{group.lasso.eval.parameters}}
}
\seealso{
  \link{\code{run.group.lasso}}
}

