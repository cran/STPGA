\name{STPGA-package}
\alias{STPGA-package}
\alias{STPGA}
\docType{package}
\title{Selection of Training Populations by Genetic Algorithm}
\description{
Can be utilized to select a test data calibrated training population in high dimensional prediction problems and assumes that the explanatory variables are observed for all of the individuals. Once a "good" training set is identified the response variable can be obtained only for this set to build a model for predicting the response in the test set.
 }
\details{
\tabular{ll}{
Package: \tab STPGA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-05-29\cr
License: \tab GPL-3\cr
}

The package is  useful for high dimensional prediction problems where per individual cost of observing / analyzing the response variable is high and therefore a small number of training examples is sought or when the candidate set from which the training set must be chosen is not representative of the test data set. 

The main function "GenAlgForSubsetSelection" uses a simple genetic algorithm to identify a training set of a specified size from a larger set of candidates which minimizes the (approximate) mean prediction error variance (PEV mean) for a known test set.
Let  \eqn{X} be the \eqn{n\times m} matrix of explanatory variables partitioned as \deqn{X=\left[
\begin{array}{c}
X_{Candidate}\\ \hline
X_{Test}
\end{array}
\right]}
where \eqn{X_{Candidate}} is the matrix of explanatory variables for the individuals in the candidate set and \eqn{X_{Test}} is the matrix of explanatory variables for the individuals in the test set. Given a training set, PEV for the ridge regression for the test set is given by \deqn{PEV^{Ridge}(X_{Test})=(1, X_{Test})((1,X_{Train})'(1,X_{Train})+\lambda I)^{-1}(1,X_{Test})'} for a choice of \eqn{\lambda > 0.} 

A more suitable numerically efficient approximation to \eqn{PEV^{Ridge}(X_{Test})} can be obtained by using the first few principal components (PCs) of the marker matrix \eqn{X} instead of \eqn{X} itself. Let  \eqn{P} be the matrix of first \eqn{k \leq min(m, n)} PCs of \eqn{X} partitioned as \deqn{P=\left[
\begin{array}{c}
P_{Candidate}\\ \hline
P_{Test}
\end{array}
\right]}
where \eqn{P_{Candidate}} is the matrix of PCs for the individuals in the candidate set and \eqn{P_{Test}} is the matrix of PC's for the individuals in the test set. Now, \eqn{PEV^{Ridge}(X_{Test})} can be approximated by \deqn{PEV^{Ridge}(X_{Test})\approx (1, P_{Test})((1,P_{Train})'(1,P_{Train})+\lambda I)^{-1}(1,P_{Test})'.}
The algorithm tries to minimize the trace of the (approximate) \eqn{PEV^{Ridge}(X_{Test}).}
}
\author{
Maintainer: Deniz Akdemir <deniz.akdemir.work@gmail.com>
}
\references{
References: 
Akdemir, Deniz. "Training population selection for (breeding value) prediction." arXiv preprint arXiv:1401.7953 (2014).
}
