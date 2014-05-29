GenAlgForSubsetSelection <-
function(P,Candidates,Test,ntoselect, npop, nelite, mutprob, niterations, lambda, plotiters=TRUE){
	InitPop<-lapply(1:npop, function(x){return(sample(Candidates, ntoselect))})
	InitPopFuncValues<-as.numeric(unlist(lapply(InitPop, FUN=function(x){PEVMEAN(Train=x, Test=Test,P=P, lambda=lambda)})))
	orderofInitPop<-order(InitPopFuncValues, decreasing=FALSE)
	ElitePop<-lapply(orderofInitPop[1:nelite], FUN=function(x){return(InitPop[[x]])})
	ElitePopFuncValues<-InitPopFuncValues[orderofInitPop[1:nelite]]
	meanvec<-c()
	for (iters in 1:niterations){
	CurrentPop<-GenerateCrossesfromElites(Elites=ElitePop,Candidates=Candidates, npop=npop, mutprob=mutprob)	
	
	CurrentPopFuncValues<-as.numeric(unlist(lapply(CurrentPop, FUN=function(x){PEVMEAN(Train=x, Test=Test,P=P, lambda=lambda)})))
	
	orderofCurrentPop<-order(CurrentPopFuncValues, decreasing=FALSE)
	ElitePop<-lapply(orderofCurrentPop[1:nelite], FUN=function(x){return(CurrentPop[[x]])})
	ElitePopFuncValues<-CurrentPopFuncValues[orderofCurrentPop[1:nelite]]
	meanvec<-c(meanvec,min(ElitePopFuncValues))
	##print(sort(ElitePop[[1]]))
	if(plotiters){plot(meanvec)}
	}
	return(ElitePop)
}
