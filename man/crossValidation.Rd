\name{crossValidation}
\alias{crossValidation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Conduct cross-validation 
}
\description{
  Conduct a cross-validation for a given classification/regression model and output the prediction results collected over the cross-validation loop. The cross-validation can be done in two ways: normal k-fold cross-validaiton (\code{batch=NULL}), or batch-wise cross-validation (\code{batch!=NULL}). The latter is particularly useful in the presence of significant intra-group heterogeneity. 
}
\usage{
crossValidation(data, label, batch = NULL, 
                method = lda, pred = predict, classify = TRUE, 
                folds = NULL, nBatch = 0, nFold = 10, 
                verbose = TRUE, seed = NULL, ...)
}
\arguments{
  \item{data}{
    a data matrix, with samples saved in rows and features in columns.
}
  \item{label}{
    a vector of response variables (i.e., group/concentration info), must be the same length as the number of samples.
}
  \item{batch}{
    a vector of sample identifications (e.g., batch/patient ID), must be the same length as the number of samples. Ideally, this should be the identification of the samples at the highest hierarchy (e.g., the patient ID rather than the spectral ID). If missing, a normal k-fold cross validaiton will be performed (i.e., the data is split randomly into k folds). Ignored if \code{folds} is given.
}
  \item{method}{
    the name of the function to be performed on training data (can be any model-based procedures, like classification/regression or even pre-processings). A user-defined function is possible, see \code{\link{fnPcaLda}} as an example.
}
  \item{pred}{
    the name of the function to be performed on testing data (eg. new substances) based on the model built by \code{method}. A user-defined function is possible, see \code{\link{predPcaLda}} as an example.
}
  \item{classify}{
    a boolean value, \code{classify=TRUE} means a classification task, otherwise a regression task. It is used in the function \code{\link{predSummary}}.
}
  \item{folds}{
    a list of indices specifying the sample index to be used in each fold, can be the output of function \code{\link{dataSplit}}. If missing, a data split will be done first before performing cross-validaiton
}
  \item{nBatch}{
    an integer, the number of data folds in case of batch-wise cross-validaiton (if \code{nBatch=0}, each batch will be used as one fold). Ignored if \code{folds} is given or if \code{batch} is missing.
}
  \item{nFold}{
    an integer, the value of k in case of normal k-fold cross-validaiton. Ignored if \code{folds} or \code{batch} is given.
}
  \item{verbose}{
    a boolean value, if or not to print out the logging info
}
  \item{seed}{
    an integer, if given, will be used as the random seed to split the data in case of k-fold cross-validation. Ignored if \code{batch} or \code{folds} is given.
}
  \item{\dots}{
    parameters to be passed to the \code{method}
}
}
\details{
  The cross-validaiton will be conducted based on the data partitions \code{folds}, each fold is predicted once using the model built on the rest folds. If \code{folds} is missing, a data split will be done first (see more in \code{\link{dataSplit}}). 

  The procedures to be performed within the cross-validation is given in the function \code{method}, for example, \code{\link{fnPcaLda}}. A user-defined function is possible, as long as the it follows the same structure as \code{\link{fnPcaLda}}. A two-layer cross-validation (see reference) can be done by using a tuning function as \code{method}, such as \code{\link{tunePcaLda}} (see examples). In this case, the parameters of a classifier are optimized using the training data within \code{\link{tunePcaLda}} and the optimal model is tested on the testing data. The parameters of pre-processing can be optimized in a similar way by involving the pre-processing steps into the function \code{method}.

 NOTE: It is recommended to specify the \code{seed} for a normal k-fold cross-validation in order to get the same results from repeated runnings. 
}
\value{
  A list with elements 
  \item{Fold}{a list, each giving the sample indices of a fold}
  \item{True}{a vector of characters, the groundtruth response variables, collected for each fold when it is used as testing data}
  \item{Pred}{a vector of characters, the results from prediction, collected for each fold when it is used as testing data}
  \item{Summ}{a list, the output of function \code{\link{predSummary}}. A confusion matrix (if \code{classify=TRUE}) from \code{\link[caret]{confusionMatrix}} or RMSE (if \code{classify=FALSE}) calculated from each fold being predicted.}
}
\references{
  S. Guo, T. Bocklitz, et al., Common mistakes in cross-validating classification models. Analytical methods 2017, 9 (30): 4410-4417.
}
\author{
  Shuxia Guo, Thomas Bocklitz, Juergen Popp
}

\seealso{
  \code{\link{dataSplit}}
}
\examples{
  data(DATA)
  ### perform batch-wise cross-validation using the function fnPcaLda
  RES3 <- crossValidation(data=DATA$spec
                          ,label=DATA$labels
                          ,batch=DATA$batch
                          ,method=fnPcaLda
                          ,pred=predPcaLda
                          ,folds=NULL 
                          ,nBatch=0
                          ,nFold=3
                          ,verbose=TRUE     
                          ,seed=NULL
                          
                          ### parameters to be passed to fnPcaLda
                          ,center=TRUE
                          ,scale=FALSE
   )


   ### perform a two-layer cross-validation using the function tunePcaLda,
   ### where the number of principal components used for LDA is optimized 
   ### (i.e., internal cross-validaiton).
   RES4 <- crossValidation(data=DATA$spec
                          ,label=DATA$labels				    
                          ,batch=DATA$batch	
                          ,method=tunePcaLda 
                          ,pred=predPcaLda     
                          ,folds=NULL      
                          ,nBatch=0			    
                          ,nFold=3					
                          ,verbose=TRUE     
                          ,seed=NULL
                          
                          ### parameters to be passed to tunePcaLda
                          ,nPC=2:4
                          ,cv=c('CV', 'BV')[2]
                          ,nPart=0
                          ,optMerit=c('Accuracy', 'Sensitivity')[2]
                          ,center=TRUE
                          ,scale=FALSE
  )
}