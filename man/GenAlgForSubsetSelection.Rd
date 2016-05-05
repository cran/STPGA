\name{GenAlgForSubsetSelection}
\alias{GenAlgForSubsetSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Genetic algorithm for subset selection}
\description{
It uses a genetic algorithm to select \eqn{n_{Train}} individuals so that optimality criterion is minimum.
}

\usage{
GenAlgForSubsetSelection(P, Candidates, Test, ntoselect, 
npop, nelite, mutprob, niterations, lambda, plotiters=TRUE,errorstat="PEVMEAN")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
\eqn{n \times k} matrix of the first PCs of the predictor variables. The matrix needs to have union of the identifiers of the candidate and test individuals as rownames. 
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
 scalar shrinkage parameter  (\eqn{\lambda>0}).
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
npop=100, nelite=5, mutprob=.8, niterations=20, lambda=1e-5,
plotiters=FALSE,errorstat="PEVMEAN")


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

 \dontrun{
data(iris)
#We will try to estimate petal width from
#variables sepal length and width and petal length.
y<-iris[,4]
X<-as.matrix(iris[,1:3])
X<-cbind(rep(1,nrow(X)),X)
names(y)<-rownames(X)<-paste(iris[,5], rep(1:50,3),sep="_" )
#test data 25 iris plants selected at random from the virginica family,
#candidates are the plants in the  setosa and versicolor families.
candidates<-rownames(X)[1:100]
test<-sample(setdiff(rownames(X),candidates), 25)
###test sample
ytest<-y[names(y)\%in\%test]
Xtest<-X[rownames(X)\%in\%test,]
#NOTE: Increase niterations and npop substantially for better convergence.
ListTrain2<-GenAlgForSubsetSelection(P=X,Candidates=candidates,
Test=test,ntoselect=25,npop=100,
nelite=5, mutprob=.8, niterations=200,
lambda=1e-5, errorstat="PEVMEAN2")

ListTrain<-GenAlgForSubsetSelection(P=X,Candidates=candidates,
Test=test,ntoselect=25,npop=100,
nelite=5, mutprob=.8, niterations=200, lambda=1e-5, errorstat="PEVMEAN")

ListTrain3<-GenAlgForSubsetSelection(P=X,Candidates=candidates,
Test=test,ntoselect=25,npop=100, 
nelite=5, mutprob=.8, niterations=200, lambda=1e-5, errorstat="CDMEAN")
ListTrain4<-GenAlgForSubsetSelection(P=X,Candidates=candidates,
Test=test,ntoselect=25,npop=100,
nelite=5, mutprob=.8, niterations=200, lambda=1e-5, errorstat="CDMEAN2")

ListTrain5<-GenAlgForSubsetSelection(P=X,Candidates=candidates,
Test=test,ntoselect=25,npop=100,
nelite=5, mutprob=.8, niterations=200, lambda=1e-5, errorstat="PEVMEAN0")


ListTrain6<-GenAlgForSubsetSelection(P=X,Candidates=candidates,
Test=test,ntoselect=25,npop=100,
nelite=5, mutprob=.8, niterations=200, lambda=1e-5, errorstat="CDMEAN0")


ListTrain7<-GenAlgForSubsetSelection(P=X,Candidates=candidates,
Test=test,ntoselect=25, npop=100,
nelite=5, mutprob=.8, niterations=200, lambda=1e-5, errorstat="DOPT")


ListTrain8<-GenAlgForSubsetSelection(P=X,Candidates=candidates,
Test=test,ntoselect=25, npop=100,
nelite=5, mutprob=.8, niterations=200, lambda=1e-5, errorstat="AOPT")


ListTrain9<-GenAlgForSubsetSelection(P=X,Candidates=candidates,
Test=test,ntoselect=25, npop=100,
nelite=5, mutprob=.8, niterations=200, lambda=1e-5, errorstat="PEVMAX")


ListTrain10<-GenAlgForSubsetSelection(P=X,Candidates=candidates,
Test=test,ntoselect=25, npop=100,
nelite=5, mutprob=.8, niterations=200, lambda=1e-5, errorstat="CDMAX2")

#
#library(kernlab)
#kpc <- kpca(~.,data=data.frame(X),kernel="rbfdot",
 #           kpar=list(sigma=0.2),features=2)
#print the principal component vectors

#plot(rotated(kpc),col=as.integer(iris[,5]), 
#xlab="1st Principal Component",ylab="2nd #Principal Component")

#ListTrain11<-GenAlgForSubsetSelection(P=rotated(kpc),Candidates=candidates,
#Test=test,ntoselect=25,npop=100,
#nelite=5, mutprob=.8, niterations=200, #lambda=1e-5, errorstat="PEVMAX")

#


##predictions by optimized sample
ytrainopt<-y[names(y)\%in\% ListTrain[[1]]]
Xtrainopt<-X[rownames(X)\%in\%ListTrain[[1]],]
modelopt<-lm(ytrainopt~-1+Xtrainopt)
predictopt<-Xtest\%*\%modelopt$coefficients

ytrainopt2<-y[names(y)\%in\% ListTrain2[[1]]]
Xtrainopt2<-X[rownames(X)\%in\%ListTrain2[[1]],]
modelopt2<-lm(ytrainopt2~-1+Xtrainopt2)
predictopt2<-Xtest\%*\%modelopt2$coefficients


ytrainopt3<-y[names(y)\%in\% ListTrain3[[1]]]
Xtrainopt3<-X[rownames(X)\%in\%ListTrain3[[1]],]
modelopt3<-lm(ytrainopt3~-1+Xtrainopt3)
predictopt3<-Xtest\%*\%modelopt3$coefficients

ytrainopt4<-y[names(y)\%in\% ListTrain4[[1]]]
Xtrainopt4<-X[rownames(X)\%in\%ListTrain4[[1]],]
modelopt4<-lm(ytrainopt4~-1+Xtrainopt4)
predictopt4<-Xtest\%*\%modelopt4$coefficients


ytrainopt5<-y[names(y)\%in\% ListTrain5[[1]]]
Xtrainopt5<-X[rownames(X)\%in\%ListTrain5[[1]],]
modelopt5<-lm(ytrainopt5~-1+Xtrainopt5)
predictopt5<-Xtest\%*\%modelopt5$coefficients


ytrainopt6<-y[names(y)\%in\% ListTrain6[[1]]]
Xtrainopt6<-X[rownames(X)\%in\%ListTrain6[[1]],]
modelopt6<-lm(ytrainopt6~-1+Xtrainopt6)
predictopt6<-Xtest\%*\%modelopt6$coefficients


ytrainopt7<-y[names(y)\%in\% ListTrain7[[1]]]
Xtrainopt7<-X[rownames(X)\%in\%ListTrain7[[1]],]
modelopt7<-lm(ytrainopt7~-1+Xtrainopt7)
predictopt7<-Xtest\%*\%modelopt7$coefficients


ytrainopt8<-y[names(y)\%in\% ListTrain8[[1]]]
Xtrainopt8<-X[rownames(X)\%in\%ListTrain8[[1]],]
modelopt8<-lm(ytrainopt8~-1+Xtrainopt8)
predictopt8<-Xtest\%*\%modelopt8$coefficients


ytrainopt9<-y[names(y)\%in\% ListTrain9[[1]]]
Xtrainopt9<-X[rownames(X)\%in\%ListTrain9[[1]],]
modelopt9<-lm(ytrainopt9~-1+Xtrainopt9)
predictopt9<-Xtest\%*\%modelopt9$coefficients


ytrainopt10<-y[names(y)\%in\% ListTrain10[[1]]]
Xtrainopt10<-X[rownames(X)\%in\%ListTrain10[[1]],]
modelopt10<-lm(ytrainopt10~-1+Xtrainopt10)
predictopt10<-Xtest\%*\%modelopt10$coefficients

#ytrainopt11<-y[names(y)\%in\% ListTrain11[[1]]]
#Xtrainopt11<-X[rownames(X)\%in\%ListTrain11[[1]],]
#modelopt11<-lm(ytrainopt11~-1+Xtrainopt11)
#predictopt11<-Xtest\%*\%modelopt11$coefficients


###predictions by a random sample of the same size
rs<-sample(candidates, 25)
ytrainrs<-y[names(y)\%in\%rs]
Xtrainrs<-X[rownames(X)\%in\%rs,]
modelrs<-lm(ytrainrs~-1+Xtrainrs)
predictrs<-Xtest\%*\%modelrs$coefficients
#accuracies of the optimized sample and random sample.


#(expect optimized sample to have better accuracies than a random sample)
cor(predictopt,ytest)
cor(predictopt2,ytest)
cor(predictopt3,ytest)
cor(predictopt4,ytest)
cor(predictopt5,ytest)
cor(predictopt6,ytest)
cor(predictopt7,ytest)
cor(predictopt8,ytest)
cor(predictopt9,ytest)
cor(predictopt10,ytest)
#cor(predictopt11,ytest)

cor(predictrs, ytest)


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

