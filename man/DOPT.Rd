\name{DOPT}
\alias{DOPT}

\title{DOPT}
\description{
Calculate  \deqn{-det((P'_{Train}P_{Train}+lambda*I))} 
}
\usage{
DOPT(Train, Test, P, lambda = 1e-05)
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
