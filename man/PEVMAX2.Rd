\name{PEVMAX2}
\alias{PEVMAX2}


\title{PEVMAX2}

\description{
Calculate  \deqn{max(diag(CP_{Test}(P'_{Train}P_{Train}+lambda*I)^{-1}P'_{Train}P_{Train}(P'_{Train}P_{Train}+\lambda*I)^{-1}P'_{Test}C')} 
}
\usage{
PEVMAX2(Train, Test, P, lambda = 1e-05, C=NULL)
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
contrast matrix
}
}

\value{
value of the criterion.
}

\author{
Deniz Akdemir
}
