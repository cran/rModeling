
#' @export
crossValidation <- function(data     			    ## spectra
                            ,label				    ## vector of response variables
                            ,batch=NULL			  ## vector of batch information, used for batchwise validation
                            ,method=lda       ## method used for modeling
                            ,pred=predict     ## method used for prediction
                            ,classify=TRUE    ## if TRUE, classification, else regression
                            ,folds=NULL       ## if provided, cross-validation will be performed based on this partition
                            ,nBatch=0			    ## number of folds for batchwise validation (ignored if Batch=NULL)
                            ,nFold=10					## number of folds for normal cross validation (ignored if Batch supplied)
                            ,verbose=TRUE     ## flag about logging info
                            ,seed=NULL        ## seed for ramdomarize spectral indices, only used for normal CV
                            ,...              ## parameters used in method of modeling
                            )
{
  dimData <- dim(data)
  
  if(length(label)!=dimData[1])
  {
    stop('Error: numbers of spectra and labels do not match!')
  }
  
  indices <- 1:dimData[1]
  
  if(is.null(folds))    ## data partition is not provided
  {
    folds <- dataSplit(indices, batch, nBatch, nFold, verbose, seed)
  }
  
  nFold <- length(folds)
  
  allTrue <- rep(NA,length(label))
  allPred <- rep(NA,length(label))
  allSumm <- list()
  for (i in 1:nFold)
  {
    ixTrain <- indices[-folds[[i]]]
    
    model <- method(data[ixTrain,], label[ixTrain], batch=batch[ixTrain], ...)
    tmpPred <- pred(model, data[folds[[i]],,drop=F])
    
    allTrue[folds[[i]]] <- label[folds[[i]]]
    allPred[folds[[i]]] <- tmpPred
    
    if(classify)
    {
      allSumm[[i]] <- predSummary(reference=label[folds[[i]]] 
                                  ,prediction=tmpPred
                                  ,lev=unique(label)
                                  ,classify=TRUE
                                  )
    }
    else 
    {
      allSumm[[i]] <- predSummary(reference=label[folds[[i]]] 
                                  ,prediction=tmpPred
                                  ,lev=NULL
                                  ,classify=FALSE
                                  )
    }
  }
  
  return(list(Fold=folds, True=allTrue, Pred=allPred, Summ=allSumm))
}

#' @export
tunePcaLda <- function(data                                           ### Spectra
                       ,label                                         ### response variables
                       ,batch=NULL                                    ### batch information
                       ,nPC=1:50                                      ### number of PCs used for PCA+LDA
                       ,optMerit=c('Accuracy', 'Sensitivity')[2]      ### merit to be optimized
                       ,maximize=TRUE                                 ### maximize merit, otherwise minimize
                       ,cv=c('CV', 'BV')[2]                           ### parameter tuning by normal or batchwise cross-validaiton
                       ,nPart=10                                      ### number of folds to split for cross-validation
                       ,...                                           ### parameters used in crossValidaiton
                       )
{
  merits <- c()
  for(ncomp in nPC)
  {
    switch(cv,
           'CV'={
             RESULT <- crossValidation(data=data
                                       ,label=label
                                       ,batch=NULL
                                       ,method=fnPcaLda
                                       ,pred=predPcaLda
                                       ,nFold=nPart
                                       ,cv='none'
                                       ,nPC=ncomp
                                       ,...
             )
           },
           'BV'={
             
             if(is.null(batch))
             {
               stop('Error: Batch information needed for batch-wise cross-validation!')
             }
             
             RESULT <- crossValidation(data=data
                                       ,label=label
                                       ,batch=batch
                                       ,nBatch=nPart
                                       ,method=fnPcaLda
                                       ,pred=predPcaLda
                                       ,cv='none'
                                       ,nPC=ncomp
                                       ,...
             )
           }
    )
    switch(optMerit,
           'Accuracy' = {
             SUMM <- predSummary(reference=RESULT$True
                                 ,prediction=RESULT$Pred
                                 ,lev=unique(label)
                                 ,classify=TRUE
                                 )
             merits <- c(merits, SUMM$overall['Accuuracy'])
           },
           'Sensitivity' = {
             SUMM <- predSummary(reference=RESULT$True
                                 ,prediction=RESULT$Pred
                                 ,lev=unique(label)
                                 ,classify=TRUE
                                 )
             if(length(unique(label)>2)) 
             {
               merits <- c(merits, mean(SUMM$byClass[, 'Sensitivity']))
             }
             else
             {
                merits <- c(merits, 0.5*(SUMM$byClass['Sensitivity']+SUMM$byClass['Specificity']))
             }
           }
           )
  }
  
  if(maximize) optPC <- nPC[which.max(merits)]
  else optPC <- nPC[which.min(merits)]
  
  PCA <- prcomp(data, ...)
  LDA <- MASS::lda(PCA$x[, 1:optPC], label, ...)
  RESULT <- list(PCA=PCA, LDA=LDA, nPC=optPC)
  
  return(RESULT)
}

#' @export
fnPcaLda <- function(data                             ### spectra
                     ,label                           ### response variables
                     ,batch=NULL                      ### batch information
                     ,nPC=10                          ### number of PCs used in PCA-LDA
                     ,cv=c('none', 'CV', 'BV')[1]     ### type of cross-validaiton
                     ,nPart=10                        ### number of folds for cross-validaiton, ignored if cv='none'
                     ,...)
{
  switch(cv,
         'none'={
           PCA <- prcomp(data, ...)
           LDA <- MASS::lda(PCA$x[, 1:nPC], label, ...)
           RESULT <- list(PCA=PCA, LDA=LDA, nPC=nPC)
         },
         'CV'={
           RESULT <- crossValidation(data=data
                                     ,label=label
                                     ,batch=NULL
                                     ,method=fnPcaLda
                                     ,pred=predPcaLda
                                     ,nFold=nPart
                                     ,cv='none'
                                     ,nPC=nPC
                                     ,...
                                     )
         },
         'BV'={
           if(is.null(batch))
           {
             stop('Error: Batch information needed for batch-wise cross-validation!')
           }
           RESULT <- crossValidation(data=data
                                     ,label=label
                                     ,batch=batch
                                     ,nBatch=nPart
                                     ,method=fnPcaLda
                                     ,pred=predPcaLda
                                     ,cv='none'
                                     ,nPC=nPC
                                     ,...
                                     )
         }
         )
  
  return(RESULT)
}

#' @export
predPcaLda <- function(objModel, newData)
{
  pScores <- predict(objModel$PCA, newData)[, 1:objModel$nPC]
  pred <- as.character(predict(objModel$LDA, pScores)$class)
  return(pred)
}

### function of data partition
#' @export
dataSplit <- function(ixData     			  ## index of spectra
                      ,batch=NULL			  ## Batch vector for batchwise validation
                      ,nBatch=0			    ## Number of folds for batchwise validation (ignored if Batch=NULL)
                      ,nFold=10					## Number of folds for normal cross validation (ignored if Batch supplied)
                      ,verbose=TRUE     ## flag about logging info
                      ,seed=NULL        ## seed for ramdomarize spectral indices, only used for normal CV
                      )
{
  if (is.null(batch))       # there's no batch vector used for BV
  {
    if (verbose)
    {
      print(paste("CV", nFold))
    }
    folds <- vector("list", nFold)
    
    if(!is.null(seed)) set.seed(seed)
    
    indices <- sample(ixData, length(ixData))
    
    nRest <- length(ixData)%%nFold       # rest after partition
    nPart <- length(ixData)%/%nFold      # partition
    
    for (i in 1:nFold)               # distribute the rest spectra
    {
      folds[[i]] <- indices[(1:nPart)+(i-1)*nPart]
      
      if(nRest>0)
      {
        folds[[i]] <- c(folds[[i]], indices[nFold*nPart+nRest])
        nRest <- nRest-1
      }
    }
    if (verbose)
    {
      print('CV: data split finished')
    }
  }
  else	
  {
    if(length(batch)!=length(ixData))
    {
      stop('Error: numbers of spectra and batch do not match!')
    }
    
    uniBatch <- unique(batch)
    if(nBatch>1)
    {
      if (verbose)
      {
        print(paste("BV", nBatch))
      }
      levels <- vector("list", nBatch)
      bchIndex <- 1:length(uniBatch)
      for (i in 1:nBatch)
      {
        levels[[i]] <- uniBatch[which(bchIndex%%nBatch==(i-1))]  # to partion "batch" in nBatch parts
      }
      folds <- vector("list",nBatch)
      for (i in 1:nBatch)
      {
        folds[[i]] <- which(is.element(el=batch, set=levels[[i]]))
      }
    }
    else
    {
      nBatch <- length(uniBatch)
      
      if(verbose)
      {
        print(paste("BV", nBatch))
      }
      
      folds <- vector("list", nBatch)
      for (i in 1:nBatch)
      {
        folds[[i]] <- which(batch==uniBatch[i])
      }
    }
    if (verbose)
    {
      print('BV: data split finished')
    }
  }
  return(folds)
}

#' @export
predSummary <- function(reference, prediction, lev=NULL, classify=TRUE)
{
  if(length(reference)!=length(prediction))
  {
    stop('Error: length of reference and prediction must be the same!')
  }
  
  if(classify)
  {
    if(is.null(lev)) lev <- unique(reference)
    
    reference <- factor(reference, levels=lev)
    prediction <- factor(prediction, levels=lev)
    
    return(caret::confusionMatrix(prediction, reference))
  }
  else
  {
    return(sqrt(sum((reference-prediction)^2)))
  }
}


#' @import e1071
#' @importFrom caret confusionMatrix
#' @importFrom MASS lda
#' @import stats 