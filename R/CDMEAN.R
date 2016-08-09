
CDMEAN <-
  function(Train,Test, P, lambda=1e-5, C=NULL){
    PTrain<-P[rownames(P)%in%Train,]
    PTest<-P[rownames(P)%in%Test,]
    if (!is.null(C)){ PTest<-C%*%PTest}
    CDmean<-mean(diag(PTest%*%solve(crossprod(PTrain)+lambda*diag(ncol(PTrain)),t(PTest)))/diag(tcrossprod(PTest)))
    return(CDmean)
  }
