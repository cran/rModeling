\name{predSummary}
\alias{predSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate the merit of the prediction
}
\description{
  produce the confusion matrix using function \code{\link[caret]{confusionMatrix}} from package caret if \code{classify=TRUE}, otherwise calculate the RMSE between the predicted and groundtruth values.
}
\usage{
  predSummary(reference, prediction, 
              lev = NULL, classify = TRUE)
}
\arguments{
  \item{reference}{
    groundtruth values.
}
  \item{prediction}{
    predicted values.
}
  \item{lev}{
    a vector of character, specifying the group names. Ignored if \code{classify=FALSE}.
}
  \item{classify}{
    a boolean value, telling whether a classification or regression task. 
}
}

\value{
 If \code{classify=TRUE}, a list, the output from \code{\link[caret]{confusionMatrix}}
 Otherwise a numeric value, giving the RMSE of the prediction.
}

\author{
  Shuxia Guo, Thomas Bocklitz, Juergen Popp
}

\seealso{
  \code{\link[caret]{confusionMatrix}}
}