\name{predPcaLda}
\alias{predPcaLda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Predict new instances using the PCA-LDA model built from \code{\link{tunePcaLda}}).
}
\description{
  Predict new instances given in \code{newData} using the PCA-LDA model \code{objModel} built from \code{\link{tunePcaLda}}).
}

\usage{
 predPcaLda(objModel, newData)
}
\arguments{
  \item{objModel}{
    the classifier built from \code{\link{tunePcaLda}}).
}
  \item{newData}{
    data matrix composed of samples to be predicted.
}
}
\value{
  a vector of characters composed of the output of the prediction.
}
\author{
  Shuxia Guo, Thomas Bocklitz, Juergen Popp
}