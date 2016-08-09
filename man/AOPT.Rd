\name{AOPT}
\alias{AOPT}

\title{AOPT}

\description{
Calculate  \deqn{trace(C(P'_{Train}P_{Train}+\lambda*I)^{-1}C').} C is taken as identity matrix if not specified. 
}
\usage{
AOPT(Train, Test, P, lambda = 1e-05, C=NULL)
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
