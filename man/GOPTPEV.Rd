\name{GOPTPEV}
\alias{GOPTPEV}

\title{GOPTPEV}

\description{
Calculate  \deqn{max(eigenval(CP_{Test}(P_{Train}'P_{Train}+\lambda*I_{ntrain})^{-1}P_{Test}'C')).} 
}
\usage{
GOPTPEV(Train, Test, P, lambda = 1e-05, C=NULL)
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
contrast matrix.
}
}

\value{
value of the criterion.
}

\author{
Deniz Akdemir
}
