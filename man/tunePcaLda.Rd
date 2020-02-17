\name{tunePcaLda}
\alias{tunePcaLda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Build a classifier with parameter tuning.
}
\description{
  optimize the number of principal component to be used in LDA based on a cross-validation procedure.
}
\usage{
 tunePcaLda(data, label, batch = NULL, nPC = 1:50, 
            optMerit = c("Accuracy", "Sensitivity")[2], 
            maximize = TRUE, 
            cv = c("CV", "BV")[2], 
            nPart = 10, ...)
}
\arguments{
  \item{data}{
    a data matrix, with samples saved in rows and features in columns.
}
  \item{label}{
    a vector of response variables (i.e., group/concentration info), must be the same length as the number of samples.
}
  \item{batch}{
    a vector of batch variables (i.e., batch/patient ID), must be given in case of \code{cv='BV'}. Ideally, this should be the identification of the samples at the highest hierarchy (e.g., the patient ID rather than the spectral ID). Ignored for \code{cv='CV'}.
}
  \item{nPC}{
    a vector of integers, the candidate numbers of principal components to be used for LDA, out of which an optimal value will be selected. 
}
  \item{optMerit}{
    a character value, the name of the merit to be optimized. The mean sensitivity will be optimized if \code{optMerit = "Sensitivity"}.
}
  \item{maximize}{
    a boolean value, if or not maximize the merit.
}
  \item{cv}{
    a character value, specifying the type of cross-validation. 
}
  \item{nPart}{
    an integer, the number of folds to be split for cross-validation. Equivelant to \code{nFold} of \code{\link{crossValidation}} for \code{cv='CV'} and to \code{nBatch} for \code{cv='BV'}. (NOTE: use \code{nPart=0} for leave-one-batch out cross-validaiton). 
}
  \item{\dots}{
    parameters for \code{\link{crossValidation}}
}
}
\details{
   build a classifier using each value in \code{nPC}, of which the performance is evaluated with a normal k-fold or batch-wise cross-validation. The optimal number is selected as the one giving the maximal (\code{maximize=TRUE}) or minimal (\code{maximize=FALSE}) merit. 

   A two-layer cross-validation can be performed by using \code{tunePcaLda} as the \code{method} in \code{\link{crossValidation}}.
}
\value{
  A list of elements:
  \item{PCA}{PCA model}
  \item{LDA}{LDA model built with the optimal number of principal components}
  \item{nPC}{the optimal number of principal components}
}
\references{
  S. Guo, T. Bocklitz, et al., Common mistakes in cross-validating classification models. Analytical methods 2017, 9 (30): 4410-4417.
}
\author{
  Shuxia Guo, Thomas Bocklitz, Juergen Popp
}

\seealso{
  \code{\link{crossValidation}}, \code{\link{tunePcaLda}}, \code{\link[MASS]{lda}}, \code{\link[stats]{prcomp}}
}
\examples{
  data(DATA)
  ### perform parameter tuning with a 3-fold cross-validaiton
  RES2 <- tunePcaLda(data=DATA$spec
                   ,label=DATA$labels
                   ,batch=DATA$batch
                   ,nPC=2:4
                   ,cv=c('CV', 'BV')[1]
                   ,nPart=3
                   ,optMerit=c('Accuracy', 'Sensitivity')[2]
                   ,center=TRUE
                   ,scale=FALSE)
}