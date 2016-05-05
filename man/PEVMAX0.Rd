\name{PEVMAX0}
\alias{PEVMAX0}
\title{PEVMAX0}
\description{
Calculate  \deqn{max(diag(P_{Train}(P'_{Train}P_{Train}+lambda*I)^{-1}P'_{Train}))} 
}
\usage{
PEVMAX0(Train, Test, P, lambda = 1e-05)
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
}

\value{
value of the criterion.
}

\author{
Deniz Akdemir
}
