\name{select.top.features}
\alias{select.top.features}
\title{Function to select top features}
\usage{
  select.top.features(features, labels,
    feature.count = 10000, min.example.fraction = 0.01)
}
\arguments{
  \item{features}{Feature matrix}

  \item{labels}{Feature labels (Should be +1/-1)}

  \item{feature.count}{Number of features to select.
  Default is 10000.}

  \item{min.example.fraction}{Minimal fraction of examples
  that a feature should defined in. Default is 1\%.}
}
\value{
  Vector of selected feature indexes.
}
\description{
  The differences are determined by weighted means. The
  number of examples in which a particular feature is
  defined is used to estimate the means in both classes.
}

