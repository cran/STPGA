GenAlgForSubsetSelectionNoTest<-function (P, ntoselect, npop, nelite, mutprob, mutintensity=2,
          niterations, lambda, plotiters = TRUE, errorstat="PEVMEAN",C=NULL, mc.cores=1,  InitPop=NULL, tolconv=1e-5, Vg=NULL, Ve=NULL) {

  if (errorstat=="CDMEANMM"){
    K=solve(P)
  }
  if (errorstat=="PEVMEANMM"){
    K=solve(P)
  }
  if (errorstat=="GAUSSMEANMM"){
    K=solve(P)
  }
  
  linenames<-rownames(P)
  
  completeton<-function(x){if (length(x)<ntoselect){
    x<-c(x, sample(setdiff(unique(unlist(InitPop)), x), ntoselect-length(x)))
  }
    return(x)
  }
  if (is.null(InitPop) ){InitPop <- lapply(1:npop, function(x) {
    return(sample(linenames, ntoselect))
  })}
  else{
    InitPop<-lapply(InitPop, completeton)
  }
  
  
  if (errorstat=="PEVMEAN"){
   InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
    PEVMEAN(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
  },mc.cores=mc.cores,mc.preschedule=T)))
  }
  
  else if (errorstat=="PEVMEANMM"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      PEVMEANMM(Train = x, Test = setdiff(linenames,x), Kinv = P, lambda = lambda, C=C, Vg=Vg, Ve=Ve)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  
  else if (errorstat=="CDMEANMM"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      CDMEANMM(Train = x, Test = setdiff(linenames,x), Kinv = P,K=K, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  else if (errorstat=="GAUSSMEANMM"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      GAUSSMEANMM(Train = x, Test = setdiff(linenames,x), Kinv = P,K=K, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  else if (errorstat=="GOPTPEV"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      GOPTPEV(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  else if (errorstat=="GOPTPEV2"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      GOPTPEV2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  else if (errorstat=="PEVMEAN0"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      PEVMEAN0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  else if (errorstat=="PEVMEAN2"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      PEVMEAN2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  
  else if (errorstat=="PEVMAX"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      PEVMAX(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  
  else if (errorstat=="PEVMAX0"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      PEVMAX0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  
  else if (errorstat=="PEVMAX2"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      PEVMAX2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  else if (errorstat=="CDMEAN"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      CDMEAN(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  else if (errorstat=="CDMEAN0"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      CDMEAN0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  else if (errorstat=="CDMEAN2"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      CDMEAN2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  
  else if (errorstat=="CDMAX"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      CDMAX(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  else if (errorstat=="CDMAX0"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      CDMAX0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  
  else if (errorstat=="CDMAX2"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      CDMAX2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  
  else if (errorstat=="DOPT"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      DOPT(Train = x, Test = NULL, P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  
  else if (errorstat=="AOPT"){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      AOPT(Train = x, Test = NULL, P = P, lambda = lambda, C=C)
    },mc.cores=mc.cores,mc.preschedule=T)))
  }
  
  else{
    if (!is.null(errorstat)){
    InitPopFuncValues <- as.numeric(unlist(parallel::mclapply(InitPop, FUN = function(x) {
      do.call(errorstat, list(x,setdiff(linenames,x),P,lambda, C))
    },mc.cores=mc.cores,mc.preschedule=T)))
    }
    }
  
  orderofInitPop <- order(InitPopFuncValues, decreasing = FALSE)
  ElitePop <- parallel::mclapply(orderofInitPop[1:nelite], FUN = function(x) {
    return(InitPop[[x]])
  },mc.cores=mc.cores,mc.preschedule=T)
  ElitePopFuncValues <- InitPopFuncValues[orderofInitPop[1:nelite]]
  meanvec <- c()
  for (iters in 1:niterations) {
    if (iters>200){
      maxmeans<-max(meanvec[(length(meanvec)-200):length(meanvec)])
      minmeans<-min(meanvec[(length(meanvec)-200):length(meanvec)])
      meansdiff<-maxmeans-minmeans
       if(meansdiff<tolconv){break('No change in last 200 iterations')
     }
  }
    CurrentPop <- GenerateCrossesfromElites(Elites = ElitePop, 
                                            Candidates = linenames, npop = npop, mutprob = mutprob,mc.cores=mc.cores, mutintensity=mutintensity)
    
    CurrentPop<-c(CurrentPop, ElitePop[1])
    
     if (errorstat=="PEVMEAN"){
    CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                     FUN = function(x) {
                                                       PEVMEAN(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                     },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
   
    else if (errorstat=="PEVMEANMM"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMEANMM(Train = x, Test = setdiff(linenames,x), Kinv = P,  lambda = lambda, C=C, Vg=Vg, Ve=Ve)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    else if (errorstat=="CDMEANMM"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMEANMM(Train = x, Test = setdiff(linenames,x), Kinv = P, K=K,lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    else if (errorstat=="GAUSSMEANMM"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         GAUSSMEANMM(Train = x, Test = setdiff(linenames,x), Kinv = P, K=K,lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    else if (errorstat=="GOPTPEV"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         GOPTPEV(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    else if (errorstat=="GOPTPEV2"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         GOPTPEV2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    else if (errorstat=="PEVMEAN0"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMEAN0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    else if (errorstat=="PEVMEAN2"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMEAN2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
    else if (errorstat=="PEVMAX"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMAX(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    else if (errorstat=="PEVMAX0"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMAX0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
    else if (errorstat=="PEVMAX2"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         PEVMAX2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
    else if (errorstat=="CDMEAN0"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMEAN0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
    else if (errorstat=="CDMEAN"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMEAN(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
    else if (errorstat=="CDMEAN2"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMEAN2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
    else if (errorstat=="CDMAX"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMAX(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
    else if (errorstat=="CDMAX0"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMAX0(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    else if (errorstat=="CDMAX2"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         CDMAX2(Train = x, Test = setdiff(linenames,x), P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
    else if (errorstat=="DOPT"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         DOPT(Train = x, Test = NULL, P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
    else if (errorstat=="AOPT"){
      CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                       FUN = function(x) {
                                                         AOPT(Train = x, Test = NULL, P = P, lambda = lambda, C=C)
                                                       },mc.cores=mc.cores,mc.preschedule=T)))
    }
    
    else {
      if (!is.null(errorstat)){
        CurrentPopFuncValues <- as.numeric(unlist(parallel::mclapply(CurrentPop, 
                                                         FUN = function(x) {
                                                          do.call(errorstat,list(x, setdiff(linenames,x), P, lambda, C))
                                                         },mc.cores=mc.cores,mc.preschedule=T)))
      }
    }
    orderofCurrentPop <- order(CurrentPopFuncValues, decreasing = FALSE)
    ElitePop <- lapply(orderofCurrentPop[1:nelite], FUN = function(x) {
      return(CurrentPop[[x]])
    })
    ElitePopFuncValues <- CurrentPopFuncValues[orderofCurrentPop[1:nelite]]
    meanvec <- c(meanvec, min(ElitePopFuncValues))
    if (plotiters) {
      plot(meanvec)
    }
  }
  ElitePop[[nelite+1]]<-meanvec
  
  return(ElitePop)
}


