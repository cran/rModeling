\name{rModeling-package}
\alias{rModeling-package}
\alias{rModeling}
\docType{package}
\title{
 Cross-validation in spectroscopic data
}
\description{
 A cross-validation framework, allowing for model optimization and model evaluation based on batch-wise or normal k-fold cross-validation. It is built based on the ideas in \code{S. Guo, T. Bocklitz, et al., Analytical Methods 2017, 9 (30): 4410-4417}. In applications with significant intra-group heterogeneity, the batch-wise cross-validation ensures a robust and reliable statistical modeling and model evaluation.
}
\details{
\tabular{ll}{
Package: \tab rModeling\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2020-01-23\cr
License: \tab GPL-2\cr
Depends: \tab MASS \cr
         \tab caret \cr
}

The main function is \code{\link{crossValidation}}. It can be used as an independent function for model evaluation or as a wrapper of a user-defined function to optimize the parameters of a model.
}
\author{
 Shuxia Guo, Thomas Bocklitz, Juergen Popp

Maintainer: Shuxia Guo<shuxia.guo@uni-jena.de>, Thomas Bocklitz<thomas.bocklitz@uni-jena.de>, Juergen Popp<juergen.popp@ipht-jena.de>
}
\references{
 S. Guo, T. Bocklitz, et al., Common mistakes in cross-validating classification models. Analytical methods 2017, 9 (30): 4410-4417.
}
