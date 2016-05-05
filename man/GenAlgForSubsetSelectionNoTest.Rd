\name{GenAlgForSubsetSelectionNoTest}
\alias{GenAlgForSubsetSelectionNoTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Genetic algorithm for subset selection no given test}
\description{
It uses a genetic algorithm to select \eqn{n_{Train}} individuals so that optimality criterion is minimum.
}

\usage{
GenAlgForSubsetSelectionNoTest(P, ntoselect, 
npop, nelite, mutprob, niterations, lambda, plotiters=TRUE, errorstat="PEVMEAN")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
\eqn{n \times k} matrix of the first PCs of the predictor variables. The matrix needs identifiers of the individuals as rownames. 
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
 scalar shrinkage parameter (\eqn{\lambda>0}).
}
\item{plotiters}{
 plot the convergence: TRUE or FALSE. Default is TRUE.
}
\item{errorstat}{
 optimality criterion: One of the optimality criterion. Default is "PEVMEAN".
}
}

\value{
A list of length nelite. The elements of the list are optimized training samples of size \eqn{n_{train}} and they are listed in increasing order of the optimization criterion.
}

\author{
Deniz Akdemir
}
\note{
The GA does not guarantee convergence to globally optimal solutions and it is highly recommended that the algorithm is replicated to obtain ''good'' training samples.
}


\examples{
\dontrun{
data(iris)
#We will try to estimate petal width from 
#variables sepal length and width and petal length. 
y<-iris[,4]
X<-as.matrix(iris[,1:3])
names(y)<-rownames(X)<-paste(iris[,5], rep(1:50,3),sep='_')

#test data 25 iris plants selected at random from the virginica family, 
#NOTE: Increase niterations and npop substantially for better convergence.

ListTrain<-GenAlgForSubsetSelectionNoTest(P=X,ntoselect=25, 
npop=100, nelite=5, mutprob=.8, niterations=20, plotiters=FALSE, lambda=1e-5)

###test sample
ytest<-y[!(names(y)\%in\%ListTrain[[1]])]
Xtest<-X[!(rownames(X)\%in\%ListTrain[[1]]),]



##predictions by optimized sample
ytrainopt<-y[names(y)\%in\% ListTrain[[1]]]
Xtrainopt<-X[rownames(X)\%in\%ListTrain[[1]],]

modelopt<-lm(ytrainopt~1+Xtrainopt)
predictopt<-cbind(rep(1, nrow(Xtest)),Xtest)\%*\%modelopt$coefficients

###predictions by a random sample of the same size
rs<-sample(names(y), 25)
ytrainrs<-y[names(y)\%in\%rs]
Xtrainrs<-X[rownames(X)\%in\%rs,]
modelrs<-lm(ytrainrs~1+Xtrainrs)
ytestrs<-y[!(names(y)\%in\%rs)]
Xtestrs<-X[!(rownames(X)\%in\%rs),]
predictrs<-cbind(rep(1, nrow(Xtestrs)),Xtestrs)\%*\%modelrs$coefficients

#accuracies of the optimized sample and random sample. 
#(expect optimized sample to have better accuracies than a random sample.)
cor(predictopt,ytest)
cor(predictrs, ytestrs)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

