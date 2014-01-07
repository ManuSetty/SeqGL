\name{run.clustering}
\alias{run.clustering}
\title{Hierarchical clustering of features}
\usage{
  run.clustering(features, no.groups = 200, no.cores = 1)
}
\arguments{
  \item{features}{Feature matrix for clustering. Rows are
  sequences, columns are kmers}

  \item{no.groups}{Number of desired groups}

  \item{no.cores}{Number of cores for parallel processing}
}
\value{
  List containing \code{groups} which show the feature
  assignments and \code{linkage}
}
\description{
  This function performs hierarchical clustering of
  features based on their correlation
}
\seealso{
  \link{\code{build.train.test.data}}
}

