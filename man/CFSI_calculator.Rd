% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CFSI_calculator.R
\name{CFSI_calculator}
\alias{CFSI_calculator}
\title{CFSI_calculator.}
\usage{
CFSI_calculator(x, genename, genestart = 2, genesel = NULL, onecolumn = FALSE)
}
\arguments{
\item{x}{A matrix or data frame containing individual names and gene
information.}

\item{genename}{A vector or matrix containing only all gene names.}

\item{genestart}{A numeric representing the first column of gene information,
default is 2.}

\item{genesel}{A vector or matrix containing selected gene names, default is
\code{NULL}.}

\item{onecolumn}{A logical value, if \code{TRUE}, one column represents a
gene site, default is \code{FALSE}.}
}
\value{
A list contains CFSI matrix and CIBS matrix.
}
\description{
Calculate the CFSI and CIBS between individuals.
}
\details{
Calculate the CFSI(combined full sibling index) and CIBS(combined
  identity by state score) between every two individuals in a matrix or two
  individuals from two matrices.
}
