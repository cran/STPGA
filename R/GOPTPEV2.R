
GOPTPEV2 <-
  function(Train,Test, P, lambda=1e-5, C=NULL){
    PTrain<-P[rownames(P)%in%Train,]
    PTest<-P[rownames(P)%in%Test,]
    if (!is.null(C)){PTest<-C%*%PTest}
    PEV<-PTest%*%solve(crossprod(PTrain)+lambda*diag(ncol(PTrain)),t(PTrain))
    svdD<-svd(tcrossprod(PEV),nu=0,nv=0)
    return(svdD$d[1])
  }

