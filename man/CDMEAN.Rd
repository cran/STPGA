\name{CDMEAN}
\alias{CDMEAN}

\title{CDMEAN}

\description{
Calculate  \deqn{mean(diag(CP_{Test}(P'_{Train}P_{Train}+\lambda*I)^{-1}P'_{Test}C')/.diag(CP_{Test}P'_{Test}C')).} C is taken as identity matrix if not specified.  \eqn{/.} is used to denote elementwise division.
}
\usage{
CDMEAN(Train, Test, P, lambda = 1e-05, C=NULL)
}
\arguments{
  \item{Train}{
vector of identifiers for individuals in the training set
}
  \item{Test}{
vector of identifiers for individuals in the test set
}
  \item{P}{
\eqn{n \times k} matrix of the first PCs of the predictor variables. The matrix needs to have union of the identifiers of the candidate and test individuals as rownames.
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
