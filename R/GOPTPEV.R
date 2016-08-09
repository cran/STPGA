GOPTPEV <-
  function(Train,Test, P, lambda=1e-5,C=NULL){
    PTrain<-P[rownames(P)%in%Train,]
    PTest<-P[rownames(P)%in%Test,]
    if (!is.null(C)){ PTest<-C%*%PTest}
    svdD<-svd(PTest%*%solve(crossprod(PTrain)+lambda*diag(ncol(PTrain)),t(PTest)), nu=0,nv=0)
    return(svdD$d[1])
  }



