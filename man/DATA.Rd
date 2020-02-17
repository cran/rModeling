\name{DATA}
\alias{DATA}
\docType{data}
\title{
 A Raman spectral data collected from cell lines.
}
\description{
 A Raman spectral data collected from cell lines composed of three cell types: MCF-7 ('m'), Leukocytes ('l') and Erythrocytes ('r').
}
\usage{data("DATA")}
\format{
  List of 3 elements:
  $ spec: 29 Raman spectra saved into a matrix, each row corresponding one spectrum.
  $ labels: a character vector of length 29, giving the cell type of each spectrum.
  $ batch : a character vector of length 29, giving the cultivation identification of each spectrum.
}
\references{
 U. Neugebauer, et al. Towards detection and identification of circulating tumour cells using Raman spectroscopy, Analyst 2010, 135.12: 3178-3182.
}
\examples{
 data(DATA)
}
\keyword{datasets}
