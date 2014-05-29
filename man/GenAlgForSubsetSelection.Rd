\name{GenAlgForSubsetSelection}
\alias{GenAlgForSubsetSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Genetic algorithm for subset selection}
\description{
This is the main function in this package. It uses a genetic algorithm to select \eqn{n_{Train}} individuals from a larger set of \eqn{n_{candidate}} individuals so that the trace of the approximate \eqn{PEV^{Ridge}(X_{Test})} is minimum for a given set of \eqn{n_{Test}} individuals.
}

\usage{
GenAlgForSubsetSelection(P, Candidates, Test, ntoselect, 
npop, nelite, mutprob, niterations, lambda, plotiters=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
\eqn{n \times k} matrix of the first PCs of the predictor variables. The matrix needs to have union of the identifiers of the candidate and test individuals as rownames. For low dimensional prediction problems there is no need to use the PCs, for those cases the matrix of predictor variables with identifiers as rownames can be used instead (the PEV measure will be exact).
}
  \item{Candidates}{ 
  vector of identifiers for the individuals in the candidate set.
}
  \item{Test}{
  vector of identifiers for the individuals in the test set.
}
  \item{ntoselect}{
 \eqn{n_{Train}:} number of individuals to select in the training set.
}
  \item{npop}{
 genetic algorithm parameter, number of solutions at each iteration
}
  \item{nelite}{
 genetic algorithm parameter, number of solutions selected as elite parents which will  generate the next set of solutions. 
}
  \item{mutprob}{
 genetic algorithm parameter, probability of mutation for each generated solution.
}

  \item{niterations}{
 genetic algorithm parameter, number of iterations.  
}

  \item{lambda}{
 scalar shrinkage parameter in \eqn{PEV^{Ridge}(X_{Test})} (\eqn{\lambda>0}).
}
\item{plotiters}{
 plot the convergence: TRUE or FALSE. Default is TRUE.
}
}

\value{
A list of length nelite. The elements of the list are optimized training samples of size \eqn{n_{train}} and they are listed in increasing order of the optimization criterion.
}

\author{
Deniz Akdemir
}
\note{
The GA does not guarantee convergence to globally optimal solutions and it is highly recommended that the algorithm is replicated to obtain "good" training samples. The function PEVMEAN can be used to rank several solutions.
}


\examples{
data(iris)
#We will try to estimate petal width from 
#variables sepal length and width and petal length. 
y<-iris[,4]
X<-as.matrix(iris[,1:3])
names(y)<-rownames(X)<-paste(iris[,5], rep(1:50,3),sep='_')

#test data 25 iris plants selected at random from the virginica family, 
#candidates are the plants in the  setosa and versicolor families.
candidates<-rownames(X)[1:100]
test<-sample(setdiff(rownames(X),candidates), 25)
###test sample
ytest<-y[names(y)\%in\%test]
Xtest<-X[rownames(X)\%in\%test,]

#NOTE: Increase niterations and npop substantially for better convergence.

ListTrain<-GenAlgForSubsetSelection(P=X,Candidates=candidates,Test=test,ntoselect=25, 
npop=100, nelite=5, mutprob=.8, niterations=20, lambda=1e-5)


##predictions by optimized sample
ytrainopt<-y[names(y)\%in\% ListTrain[[1]]]
Xtrainopt<-X[rownames(X)\%in\%ListTrain[[1]],]

modelopt<-lm(ytrainopt~1+Xtrainopt)
predictopt<-cbind(rep(1, nrow(Xtest)),Xtest)\%*\%modelopt$coefficients

###predictions by a random sample of the same size
rs<-sample(candidates, 25)
ytrainrs<-y[names(y)\%in\%rs]
Xtrainrs<-X[rownames(X)\%in\%rs,]
modelrs<-lm(ytrainrs~1+Xtrainrs)
predictrs<-cbind(rep(1, nrow(Xtest)),Xtest)\%*\%modelrs$coefficients

#accuracies of the optimized sample and random sample. 
#(expect optimized sample to have better accuracies than a random sample)
cor(predictopt,ytest)
cor(predictrs, ytest)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

