Ares.SNV<-function(V, R, purity, CNv, CNr, K.cluster, replicates=10, cores=1){
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.repara, V=V, R=R, purity=purity, CNv=CNv, CNr=CNr, K.cluster=K.cluster, mc.cores=cores)
    
    for(rep in 1:replicates){
      if(is.na(res.mixed.model[[rep]]$obj)){
        obj.tmp[rep]<-1e20
      }else{
        obj.tmp[rep]<-res.mixed.model[[rep]]$obj
        BIC.tmp[rep]<-res.mixed.model[[rep]]$BIC
        AIC.tmp[rep]<-res.mixed.model[[rep]]$AIC
      }
    }
    Best.rep<-which.min(obj.tmp)
    res<-res.mixed.model[[Best.rep]]
    
    ### reorder cluster
    New.cluster.order<-order(res$Omega)
    Omega.new<-res$Omega[New.cluster.order]
    Mu.new<-res$Mu[New.cluster.order]
    Pi.prior.new<-res$Pi.prior[New.cluster.order]
    Pi.new<-Posterior.Cluster.Finite.Binomial.repara(K.cluster, Pi.prior.new, V, R, purity, CNv, CNr, Omega.new, Mu.new)
    cluster.new<-apply(Pi.new,1,which.max)
    res$Omega<-Omega.new
    res$Mu<-Mu.new
    res$Pi<-Pi.new
    res$Pi.prior<-Pi.prior.new
    res$cluster<-cluster.new
    res$obj.seq<-obj.tmp
    res$BIC.seq<-BIC.tmp
    res$AIC.seq<-AIC.tmp
  }else{
    res <- tryCatch(Finite.Binomial.repara(X=1, V, R, purity, CNv, CNr, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  return(res)
}

Ares.SNP<-function(B, R, purity, CNA, CNB, K.cluster, replicates=10, cores=1){
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.SNP.C, B=B, R=R, purity=purity, CNA=CNA, CNB=CNB, K.cluster=K.cluster, mc.cores=cores)
    
    for(rep in 1:replicates){
      if(is.na(res.mixed.model[[rep]]$obj)){
        obj.tmp[rep]<-1e20
      }else{
        obj.tmp[rep]<-res.mixed.model[[rep]]$obj
        BIC.tmp[rep]<-res.mixed.model[[rep]]$BIC
        AIC.tmp[rep]<-res.mixed.model[[rep]]$AIC
      }
    }
    Best.rep<-which.min(obj.tmp)
    res<-res.mixed.model[[Best.rep]]
    
    ### reorder cluster
    New.cluster.order<-order(res$Omega)
    Alpha.new<-res$Alpha[New.cluster.order]
    Omega.new<-res$Omega[New.cluster.order]
    Mu.new<-res$Mu[New.cluster.order]
    Pi.prior.new<-res$Pi.prior[New.cluster.order]
    Pi.new<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.prior.new, B, R, purity, CNA, CNB, Alpha.new, Omega.new, Mu.new)
    cluster.new<-apply(Pi.new,1,which.max)
    res$Alpha<-Alpha.new
    res$Omega<-Omega.new
    res$Mu<-Mu.new
    res$Pi<-Pi.new
    res$Pi.prior<-Pi.prior.new
    res$cluster<-cluster.new
    res$obj.seq<-obj.tmp
    res$BIC.seq<-BIC.tmp
    res$AIC.seq<-AIC.tmp
  }else{
    res <- tryCatch(Finite.Binomial.Ares.mix.SNP.C(X=1, B, R, purity, CNA, CNB, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  return(res)
}
