\name{PEVMEAN}
\alias{PEVMEAN}

\title{Calculate  mean \eqn{PEV^{Ridge}(X_{Test})}.
}
\description{
Calculate  mean of the diagonal elements of \eqn{PEV^{Ridge}(X_{Test})} for a given training set.
}
\usage{
PEVMEAN(Train, Test, P, lambda = 1e-05)
}
\arguments{
  \item{Train}{
vector of identifiers for individuals in the training set
}
  \item{Test}{
vector of identifiers for individuals in the test set
}
  \item{P}{
\eqn{n \times k} matrix of the first PCs of the predictor variables. The matrix needs to have union of the identifiers of the candidate and test individuals as rownames. For low dimensional prediction problems there is no need to use the PCs, for those cases the matrix of predictor variables with identifiers as rownames can be used instead (the PEV measure will be exact).
}
  \item{lambda}{
scalar shrinkage parameter in \eqn{PEV^{Ridge}(X_{Test})} (\eqn{\lambda>0}).
}
}

\value{
trace of the \eqn{PEV^{Ridge}(X_{Test})} for a given set of training individuals.
}

\author{
Deniz Akdemir
}
