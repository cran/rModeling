\name{fnPcaLda}
\alias{fnPcaLda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Build a classifier using PCA-LDA
}
\description{
  a classification function based on PCA following LDA. This function can be cooperated into \code{\link{crossValidation}} by setting parameter \code{method=fnPcaLda}
}
\usage{
  fnPcaLda(data, label, batch = NULL, nPC = 10, 
           cv = c("none", "CV", "BV")[1], 
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
    a vector of batch variables (i.e., batch/patient ID), must be given in case of \code{cv='BV'}. Ideally, this should be the identification of the samples at the highest hierarchy (e.g., the patient ID rather than the spectral ID). Ignored for \code{cv='None'} or \code{cv='CV'}.
}
  \item{nPC}{
    an integer, the number of principal components to be used in LDA.
}
  \item{cv}{
    a character value, specifying the type of cross-validation. 
}
  \item{nPart}{
    an integer, the number of folds to be split for cross-validation. Equivalent to \code{nFold} of \code{\link{crossValidation}} for \code{cv='CV'} and to \code{nBatch} for \code{cv='BV'}. (NOTE: use \code{nPart=0} for leave-one-batch out cross-validaiton). Ignored for \code{cv='None'}. 
}
  \item{\dots}{
    parameters for \code{\link[stats]{prcomp}} (\code{cv='None'}) or \code{\link{crossValidation}}.
}
}
\details{
    build a classifier based on the given data and return an object including the PCA and LDA models in case of \code{cv='none'}. Otherwise, a cross-validaiton is performed if \code{cv='CV'} or \code{cv='BV'}, corresponding to normal k-fold or batch-wise cross-validation, respectively. In the latter two cases, the function returns the results of the cross-validation (i.e., the output from \code{\link{crossValidation}}.
}
\value{
  For \code{cv='none'}, a list of elements:
  \item{PCA}{PCA model}
  \item{LDA}{LDA model}
  \item{nPC}{nPC used for modeling}

  For \code{cv='CV'} or \code{cv='BV'}, a list of elements:
  \item{Fold}{a list, each giving the sample indices of a fold}
  \item{True}{a vector of characters, groundtruth response variables, collected for each fold when it is used as testing data}
  \item{Pred}{a vector of characters, predicted results, collected for each fold when it is used as testing data}
  \item{Summ}{a list, the output of function \code{\link{predSummary}}. A confusion matrix (if \code{classify=TRUE}) from \code{\link[caret]{confusionMatrix}} or RMSE (if \code{classify=FALSE}) calculated from each fold being predicted.}
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
  ### perform classification with a 3-fold cross-validaiton
  RES1 <- fnPcaLda(data=DATA$spec
                   ,label=DATA$labels
                   ,batch=DATA$batch
                   ,nPC=3
                   ,cv=c('none', 'CV', 'BV')[2]
                   ,nPart=3
                   ,center=TRUE
                   ,scale=FALSE)
 }