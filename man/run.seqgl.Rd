\name{run.seqGL}
\alias{run.seqGL}
\title{SeqGL pipeline}
\usage{
  run.seqGL(peaks.file, out.dir, data.type, org,
    span = 150, max.examples = 5000, no.cores = 1,
    no.groups = ifelse(data.type == "ChIP", 20, 100),
    dictionary.file = system.file("extdata/wildcard_dict_kmer8_mismatches2_alpha5_consecutive_mis.Rdata", package = "SeqGL"))
}
\arguments{
  \item{peaks.file}{Bed file containing the peaks. See
  Details.}

  \item{out.dir}{Path to the output directory}

  \item{data.type}{ChIP or DNase: indicating the experiment
  though which the peaks were derived}

  \item{org}{IUPAC organism code. hg19, hg18, mm10, mm9,
  mm8 are supported. Note that the corresponding BSGenome
  package has to be installed.}

  \item{span}{Width of the peaks used for analysis. Default
  is 150}

  \item{max.exmamples}{Maximum examples for training. Note
  that group membership will be determined for all
  examples}

  \item{no.cores}{Number of cores for parallel processing}

  \item{no.groups}{Number of groups to use for group lasso}

  \item{dictionary.file}{The positional wildcard dictionary
  file built using ChIPKernel. A dicionary will be built if
  it not specified.}
}
\description{
  SeqGL pipeline
}
\details{
  SeqGL results will be available in \code{out.dir}
}

