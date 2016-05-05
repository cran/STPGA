\name{STPGA-package}
\alias{STPGA-package}
\alias{STPGA}
\docType{package}
\title{Selection of Training Populations by Genetic Algorithm}
\description{
Can be utilized to select a (test data) calibrated training population in high dimensional prediction problems. Once a ''good'' training set is identified the response variable can be obtained only for this set to build a model for predicting the response (in the test set).
 }
\details{
\tabular{ll}{
Package: \tab STPGA\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2016-05-05\cr
License: \tab GPL-3\cr
}

The package is  useful for high dimensional prediction problems where per individual cost of observing / analyzing the response variable is high and therefore a small number of training examples is sought or when the candidate set from which the training set must be chosen (is not representative of the test data set). 

The function "GenAlgForSubsetSelection" uses a simple genetic algorithm to identify a training set of a specified size from a larger set of candidates which minimizes an optimization criterion (for a known test set). The function "GenAlgForSubsetSelectionNoTest" tries to identify a training set of a specified size from a larger set of candidates which minimizes an optimization criterion, no test set is specified.
Let  \eqn{P} be the \eqn{n\times m} matrix of explanatory variables (or their first few principal components)  partitioned as \deqn{P=\left[
\begin{array}{c}
P_{Candidate}\\ \hline
P_{Test}
\end{array}
\right]}
where \eqn{P_{Candidate}} is the matrix of explanatory variables for the individuals in the candidate set and \eqn{P_{Test}} is the matrix of explanatory variables for the individuals in the test set. \eqn{P_{Train}} is the set of individuals in the training set. 
}
\author{
Maintainer: Deniz Akdemir <deniz.akdemir.work@gmail.com>
}
\references{
References: 
Akdemir, Deniz. "Training population selection for (breeding value) prediction." arXiv preprint arXiv:1401.7953 (2014).
}
