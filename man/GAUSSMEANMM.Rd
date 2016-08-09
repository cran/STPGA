\name{GAUSSMEANMM}
\alias{GAUSSMEANMM}

\title{GAUSSMEANMM}

\description{Kinv is a symmetric matrix (inverse of the relationship matrix K between individuals). \eqn{Z_{Test}} and \eqn{Z_{Train}} are the design matrices for the individuals in the test and training sets correspondingly. M is calculated as \eqn{1_{ntrain}(1_{ntrain}'1_{ntrain})^{-1}1_{ntrain}'.}
Calculate  \deqn{-mean(diag(Z_{Test}KZ_{Test}'-Z_{Test}KZ_{Train}'(Z_{Train}KZ_{Train}'+\lambda*I_{ntrain})^{-1}Z_{Train}KZ_{Test}').} C is taken as identity matrix if not specified. 
}
\usage{
GAUSSMEANMM(Train, Test, Kinv, K, lambda = 1e-05, C=NULL)
}
\arguments{
  \item{Train}{
vector of identifiers for individuals in the training set
}
  \item{Test}{
vector of identifiers for individuals in the test set
}
  \item{Kinv}{
\eqn{n \times n} inverse of the relationship matrix for the individuals with row and column names.
}
  \item{K}{
\eqn{n \times n} the relationship matrix for the individuals with row and column names.
}
  \item{lambda}{
scalar shrinkage parameter (\eqn{\lambda>0}).
}
 \item{C}{
Contrast Matrix.
}
}

\value{
value of the criterion.
}

\author{
Deniz Akdemir
}
