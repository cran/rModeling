\name{dataSplit}
\alias{dataSplit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A procedure to split whole dataset into multiple folds. 
}
\description{
  the whole dataset is split into multiple folds randomly (\code{batch=NULL}) or according to the batch information (\code{batch} is specified). The number of folds are defined by \code{nFold} in the former case. In the latter case, data belonging to each batch is used as one fold if \code{nBatch=0}, otherwise the dataset is split into \code{nBatch} folds according to the batch information (i.e., data from the same batch will be used exclusively in one fold). 
}
\usage{
 dataSplit(ixData, batch = NULL, 
           nBatch = 0, nFold = 10, 
           verbose = TRUE, seed = NULL)
}
\arguments{
  \item{ixData}{
    a vector of integers, demonstrating the indices of spectra.
}
  \item{batch}{
    a vector of sample identifications (e.g., batch/patient ID), must be the same length as \code{ixData}. Ideally, this should be the identification of the samples at the highest hierarchy (e.g., the patient ID rather than the spectral ID). If missing, the data is split randomly into \code{nFold} folds. 
}
  \item{nBatch}{
    an integer, the number of data folds in case of batch-wise cross-validaiton (if \code{nBatch=0}, each batch will be used as one fold). Ignored if \code{batch} is missing. 
}
  \item{nFold}{
    an integer, the number of data folds in case of normal k-fold cross-validaiton. Ignored if \code{batch} is given.
}
  \item{verbose}{
    a boolean value, if or not to print out the logging info.
}
  \item{seed}{
    an integer, if given, will be used as the random seed to split the data in case of k-fold cross-validation. Ignored if \code{batch} is given.
}
}

\value{
  a list, of which each element representing the indices of the sample belonging to one fold.
}
\references{
  S. Guo, T. Bocklitz, et al., Common mistakes in cross-validating classification models. Analytical methods 2017, 9 (30): 4410-4417.
}
\author{
  Shuxia Guo, Thomas Bocklitz, Juergen Popp
}