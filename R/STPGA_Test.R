STPGA_Test<-function(P,Candidates,Test, ntoselect, npop, nelite, mutprob, mutintensity=1,
                       niterations, preppoolrep=200, preppoolsize=5,lambda=1e-5, plotiters = TRUE, errorstat="PEVMEAN",C=NULL, mc.cores=1, tolconv=1e-5, Vg=NULL, Ve=NULL){
  

  listofgoodsols2<-NULL
  newlengthinrep<- 0
  for (reppool in 1:preppoolrep){
    if (((newlengthinrep<=ntoselect+.02*ntoselect)&&(newlengthinrep>=ntoselect-.02*ntoselect))) {break;}
    listofgoodsolstemp<-parallel::mclapply(1:preppoolsize,function(x){return(GenAlgForSubsetSelection(P=P,Candidates=Candidates,Test=Test,ntoselect=ntoselect,npop=npop, nelite=nelite, mutprob=mutprob, mutintensity=mutintensity,niterations=min(c(10,ceiling(niterations/100))), lambda=lambda,
                                                                plotiters=F,errorstat=errorstat, C=NULL,mc.cores=1, InitPop=listofgoodsols2, tolconv=tolconv, Vg=Vg, Ve=Ve))},mc.cores=min(c(mc.cores,preppoolsize)), mc.preschedule=F)
    
    
    listofgoodsols<-vector(mode="list")
    eijk=1
    for (eij in 1:preppoolsize){
    for (nij in 1:nelite){
        listofgoodsols[[eijk]]<-listofgoodsolstemp[[eij]][[nij]]
        eijk=eijk+1
      }
    }
    
    tableofgoodsols<-table(factor(unlist(listofgoodsols), levels=rownames(P)))
    
    
    newlengthinrep<- sum(tableofgoodsols>0)
    itemfreq<-tableofgoodsols
    probs<-itemfreq/sum(itemfreq)
    listofgoodsols2<-c(listofgoodsols,lapply(1:100, function(x){sample(names(itemfreq),ntoselect, prob=probs, replace=F)}))
    
  }
  ListTrain1<-GenAlgForSubsetSelection(P=P,Candidates=Candidates,Test=Test,ntoselect=ntoselect, 
                                             npop=npop, nelite=nelite, mutprob=mutprob, mutintensity=mutintensity,niterations=niterations, lambda=lambda,
                                             plotiters=plotiters,errorstat=errorstat, C=NULL,mc.cores=mc.cores, InitPop=listofgoodsols2, tolconv=tolconv, Vg=Vg, Ve=Ve)


return(ListTrain1)
}