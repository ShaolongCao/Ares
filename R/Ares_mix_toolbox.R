library(parallel)

### reparameterize
Binomial.Likelihood.repara<-function(V, R, purity, CNv, CNr, omega, mu){
  prob <- purity*CNv*omega/(2*(1-purity)/mu+purity*(CNv*omega+CNr))
  theta <- V/R
  return(exp(-R*(theta-prob)^2/(2*theta*(1.0-theta)))/sqrt(V*(1-theta))+1e-250)
}

Posterior.Cluster.Finite.Binomial.repara<-function(K.cluster, Pi.prior,V, R, purity, CNv, CNr, Omega, Mu){
  n<-length(V)
  Pi.all<-matrix(0,n,K.cluster)
  for(k in 1:K.cluster){
    Pi.all[,k]<-Pi.prior[k]*Binomial.Likelihood.repara(V, R, purity, CNv, CNr, Omega[k], Mu[k])
  }
  Pi.all<-Pi.all/rowSums(Pi.all)
  return(Pi.all)
}

Finite.Binomial.obj.repara<-function(K.cluster, Pi.prior, V, R, purity, CNv, CNr, Omega, Mu){
  n<-length(V)
  OBJ<-matrix(0,n,K.cluster)
  for(k in 1:K.cluster){
    OBJ[,k]<-Pi.prior[k]*Binomial.Likelihood.repara(V, R, purity, CNv, CNr, omega=Omega[k], mu=Mu[k])
  }
  res<-(-sum(log(rowSums(OBJ))-0.5*log(2*pi)))
  return(res)
}

Finite.Binomial.repara<-function(X, V, R, purity, CNv, CNr, K.cluster){
  n<-length(V)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Omega and Eta
    Omega <- exp(runif(K.cluster, min=-log(20), max = log(30)))
    Mu <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Finite.Binomial.obj.repara(K.cluster, Pi.prior, V, R, purity, CNv, CNr, Omega, Mu)
    while(err>1e-4 & iter<200 & obj>0 & nondecrease){
      Pi_old<-Pi
      Mu_old<-Mu
      Omega_old<-Omega
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi<-Posterior.Cluster.Finite.Binomial.repara(K.cluster, Pi.prior, V, R, purity, CNv, CNr, Omega, Mu)
      ### M step
      for(k in 1:K.cluster){
        Pi.nk<-Pi[,k]
        tmp.res<-Solve_mixed_Binomial_repara(Pi.nk, V, R, purity, CNv, CNr)
        Mu[k]<-tmp.res$mu
        Omega[k]<-tmp.res$omega
      }
      obj<-Finite.Binomial.obj.repara(K.cluster, Pi.prior, V, R, purity, CNv, CNr, Omega, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Mu<-Mu_old
      Omega<-Omega_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    vaf.fit<-rep(0,n)
    for(i in 1:n){
      vaf.fit[i]<-purity[i]*CNv[i]*Omega[cluster[i]]/(2*(1-purity[i])/Mu[cluster[i]]+purity[i]*(CNv[i]*Omega[cluster[i]]+CNr[i]))
    }
    obj<-Finite.Binomial.obj.repara(K.cluster, Pi.prior, V, R, purity, CNv, CNr, Omega, Mu)
    BIC<-log(n)*2*K.cluster+2*obj
    AIC<-4*K.cluster+2*obj
  }else{
    tmp.res<-Solve_mixed_Binomial_repara(rep(1,length(V)), V, R, purity, CNv, CNr)
    Mu<-tmp.res$mu
    Omega<-tmp.res$omega
    obj<-Finite.Binomial.obj.repara(K.cluster, 1, V, R, purity, CNv, CNr, Omega, Mu)
    
    cluster<-rep(1,n)
    vaf.fit<-purity*CNv*Omega/(2*(1-purity)/Mu+purity*(CNv*Omega+CNr))
    BIC<-log(n)*2*K.cluster+2*obj
    AIC<-4*K.cluster+2*obj
  }
  return(list(Omega=Omega, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.repara.rep<-function(V, R, purity, CNv, CNr, K.cluster, replicates=10, cores=1){
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

### True obj
Binomial.Likelihood.true<-function(V, R, purity, CNv, CNr, omega, mu){
  theta <- purity*CNv*omega/(2*(1-purity)/mu+purity*(CNv*omega+CNr))
  return(dbinom(x=V, size=R, prob=theta))
}

Posterior.Cluster.Finite.Binomial.true<-function(K.cluster, Pi.prior,V, R, purity, CNv, CNr, Omega, Mu){
  n<-length(V)
  Pi.all<-matrix(0,n,K.cluster)
  for(k in 1:K.cluster){
    Pi.all[,k]<-Pi.prior[k]*Binomial.Likelihood.true(V, R, purity, CNv, CNr, Omega[k], Mu[k])
  }
  Pi.all<-Pi.all+1e-250
  Pi.all<-Pi.all/rowSums(Pi.all)
  return(Pi.all)
}

Finite.Binomial.obj.true<-function(K.cluster, Pi.prior, V, R, purity, CNv, CNr, Omega, Mu){
  n<-length(V)
  OBJ<-matrix(0,n,K.cluster)
  for(k in 1:K.cluster){
    OBJ[,k]<-Pi.prior[k]*Binomial.Likelihood.true(V, R, purity, CNv, CNr, omega=Omega[k], mu=Mu[k])
  }
  res<-(-sum(log(rowSums(OBJ)+1e-250)))
  return(res)
}

Finite.Binomial.Goldensection<-function(X, V, R, purity, CNv, CNr, K.cluster){
  n<-length(V)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Omega and Eta
    Omega <- exp(runif(K.cluster, min=-log(20), max = log(30)))
    Mu <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Finite.Binomial.obj.true(K.cluster, Pi.prior, V, R, purity, CNv, CNr, Omega, Mu)
    while(err>1e-4 & iter<200 & obj>0 & nondecrease){
      Pi_old<-Pi
      Mu_old<-Mu
      Omega_old<-Omega
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi<-Posterior.Cluster.Finite.Binomial.true(K.cluster, Pi.prior, V, R, purity, CNv, CNr, Omega, Mu)
      ### M step
      for(k in 1:K.cluster){
        Pi.nk<-Pi[,k]
        tmp.res<-Solve_mixed_Binomial_Goldensection(Pi.nk, V, R, purity, CNv, CNr)
        Mu[k]<-tmp.res$mu
        Omega[k]<-tmp.res$omega
      }
      obj<-Finite.Binomial.obj.true(K.cluster, Pi.prior, V, R, purity, CNv, CNr, Omega, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Mu<-Mu_old
      Omega<-Omega_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    vaf.fit<-rep(0,n)
    for(i in 1:n){
      vaf.fit[i]<-purity[i]*CNv[i]*Omega[cluster[i]]/(2*(1-purity[i])/Mu[cluster[i]]+purity[i]*(CNv[i]*Omega[cluster[i]]+CNr[i]))
    }
    obj<-Finite.Binomial.obj.true(K.cluster, Pi.prior, V, R, purity, CNv, CNr, Omega, Mu)
    BIC<-log(n)*2*K.cluster+2*obj
    AIC<-4*K.cluster+2*obj
  }else{
    tmp.res<-Solve_mixed_Binomial_Goldensection(rep(1,length(V)), V, R, purity, CNv, CNr)
    Mu<-tmp.res$mu
    Omega<-tmp.res$omega
    obj<-Finite.Binomial.obj.true(K.cluster, 1, V, R, purity, CNv, CNr, Omega, Mu)
    
    cluster<-rep(1,n)
    vaf.fit<-purity*CNv*Omega/(2*(1-purity)/Mu+purity*(CNv*Omega+CNr))
    BIC<-log(n)*2*K.cluster+2*obj
    AIC<-4*K.cluster+2*obj
  }
  return(list(Omega=Omega, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Goldensection.rep<-function(V, R, purity, CNv, CNr, K.cluster, replicates=10, cores=1){
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Goldensection, V=V, R=R, purity=purity, CNv=CNv, CNr=CNr, K.cluster=K.cluster, mc.cores=cores)
    
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
    Pi.new<-Posterior.Cluster.Finite.Binomial.true(K.cluster, Pi.prior.new, V, R, purity, CNv, CNr, Omega.new, Mu.new)
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
    res <- tryCatch(Finite.Binomial.Goldensection(X=1, V, R, purity, CNv, CNr, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  return(res)
}

### SNP method
Binomial.Likelihood.SNP<-function(B, R, purity, CNA, CNB, alpha, omega, mu){
  theta <- 1/ (1 + ( 1-purity + purity*CNA*mu)/((1-purity)*alpha + purity*CNB*omega*mu))
  return(dbinom(x=B, size=R, prob=theta))
}

Posterior.Cluster.Ares.mix.SNP<-function(K.cluster, Pi.prior,B, R, purity, CNA, CNB, Alpha, Omega, Mu){
  n<-length(R)
  Pi.all<-matrix(0,n,K.cluster)
  for(k in 1:K.cluster){
    Pi.all[,k]<-Pi.prior[k]*Binomial.Likelihood.SNP(B, R, purity, CNA, CNB, Alpha[k], Omega[k], Mu[k])
  }
  Pi.all<-Pi.all+1e-250
  Pi.all<-Pi.all/rowSums(Pi.all)
  return(Pi.all)
}

Ares.mix.SNP.obj<-function(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu){
  n<-length(R)
  OBJ<-matrix(0,n,K.cluster)
  for(k in 1:K.cluster){
    OBJ[,k]<-Pi.prior[k]*Binomial.Likelihood.SNP(B, R, purity, CNA, CNB, Alpha[k], omega=Omega[k], mu=Mu[k])
  }
  res<-(-sum(log(rowSums(OBJ)+1e-250)))
  return(res)
}

Finite.Binomial.Ares.mix.SNP<-function(X, B, R, purity, CNA, CNB, K.cluster){
  n<-length(R)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Alpha <- exp(runif(K.cluster, min=-log(5), max = log(5)))
    Omega <- exp(runif(K.cluster, min=-log(20), max = log(20)))
    Mu <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    # Alpha <- exp(seq(-log(5), log(5), length.out = K.cluster))
    # Omega <- exp(seq(-log(20), log(20), length.out = K.cluster))
    # Mu <- exp(seq(-log(20), log(20), length.out = K.cluster))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    while(err>1e-4 & iter<100 & obj>0 & nondecrease){
      Pi_old<-Pi
      Alpha_old<-Alpha
      Mu_old<-Mu
      Omega_old<-Omega
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      ### M step
      for(k in 1:K.cluster){
        Pi.nk<-Pi[,k]
        tmp.res<-Solve_Ares_mix_SNP(Pi.nk, B, R, purity, CNA, CNB)
        Alpha[k]<-tmp.res$alpha
        Mu[k]<-tmp.res$mu
        Omega[k]<-tmp.res$omega
      }
      obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Alpha<-Alpha_old
      Mu<-Mu_old
      Omega<-Omega_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    vaf.fit<-rep(0,n)
    for(i in 1:n){
      vaf.fit[i] <- 1/ (1 + ( 1-purity[i] + purity[i]*CNA[i]*Mu[cluster[i]])/((1-purity[i])*Alpha[cluster[i]] + purity[i]*CNB[i]*Omega[cluster[i]]*Mu[cluster[i]]))
    }
    obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }else{
    tmp.res<-Solve_Ares_mix_SNP(rep(1,length(B)), B, R, purity, CNA, CNB)
    Alpha<-tmp.res$alpha
    Mu<-tmp.res$mu
    Omega<-tmp.res$omega
    obj<-Ares.mix.SNP.obj(K.cluster, 1, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    
    cluster<-rep(1,n)
    vaf.fit<-1/ (1 + ( 1-purity + purity*CNA*Mu)/((1-purity)*Alpha + purity*CNB*Omega*Mu))
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega=Omega, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.SNP.rep<-function(B, R, purity, CNA, CNB, K.cluster, replicates=10, cores=1){
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.SNP, B=B, R=R, purity=purity, CNA=CNA, CNB=CNB, K.cluster=K.cluster, mc.cores=cores)
    
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
    res <- tryCatch(Finite.Binomial.Ares.mix.SNP(X=1, B, R, purity, CNA, CNB, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  return(res)
}

Finite.Binomial.Ares.mix.SNP.C<-function(X, B, R, purity, CNA, CNB, K.cluster){
  n<-length(R)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Alpha <- exp(runif(K.cluster, min=-log(5), max = log(5)))
    Omega <- exp(runif(K.cluster, min=-log(20), max = log(20)))
    Mu <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    # Alpha <- exp(seq(-log(5), log(5), length.out = K.cluster))
    # Omega <- exp(seq(-log(20), log(20), length.out = K.cluster))
    # Mu <- exp(seq(-log(20), log(20), length.out = K.cluster))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    while(err>1e-4 & iter<50 & obj>0 & nondecrease){
      Pi_old<-Pi
      Alpha_old<-Alpha
      Mu_old<-Mu
      Omega_old<-Omega
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      ### M step
      tmp.res<-Ares_mix_SNP_Mstep(Pi, B, R, purity, CNA, CNB, K.cluster)
      Alpha<-tmp.res$Alpha
      Mu<-tmp.res$Mu
      Omega<-tmp.res$Omega
      obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Alpha<-Alpha_old
      Mu<-Mu_old
      Omega<-Omega_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    vaf.fit<-rep(0,n)
    for(i in 1:n){
      vaf.fit[i] <- 1/ (1 + ( 1-purity[i] + purity[i]*CNA[i]*Mu[cluster[i]])/((1-purity[i])*Alpha[cluster[i]] + purity[i]*CNB[i]*Omega[cluster[i]]*Mu[cluster[i]]))
    }
    obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }else{
    tmp.res<-Solve_Ares_mix_SNP(matrix(1,length(B),1), B, R, purity, CNA, CNB)
    Alpha<-tmp.res$alpha
    Mu<-tmp.res$mu
    Omega<-tmp.res$omega
    obj<-Ares.mix.SNP.obj(K.cluster, 1, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    
    cluster<-rep(1,n)
    vaf.fit<-1/ (1 + ( 1-purity + purity*CNA*Mu)/((1-purity)*Alpha + purity*CNB*Omega*Mu))
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega=Omega, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.SNP.C.rep<-function(B, R, purity, CNA, CNB, K.cluster, replicates=10, cores=1){
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

Finite.Binomial.Ares.mix.SNP.Newton<-function(X, B, R, purity, CNA, CNB, K.cluster){
  n<-length(R)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Alpha <- exp(runif(K.cluster, min=-log(5), max = log(5)))
    Omega <- exp(runif(K.cluster, min=-log(20), max = log(20)))
    Mu <- exp(runif(K.cluster, min=-log(20), max = log(20)))
    # Alpha <- exp(seq(-log(5), log(5), length.out = K.cluster))
    # Omega <- exp(seq(-log(20), log(20), length.out = K.cluster))
    # Mu <- exp(seq(-log(20), log(20), length.out = K.cluster))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    while(err>1e-4 & iter<100 & obj>0 & nondecrease){
      Pi_old<-Pi
      Alpha_old<-Alpha
      Mu_old<-Mu
      Omega_old<-Omega
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      ### M step
      tmp.res<-Ares_mix_SNP_Mstep_Newton(Pi, B, R, purity, CNA, CNB, K.cluster)
      Alpha<-tmp.res$Alpha
      Mu<-tmp.res$Mu
      Omega<-tmp.res$Omega
      obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Alpha<-Alpha_old
      Mu<-Mu_old
      Omega<-Omega_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    vaf.fit<-rep(0,n)
    for(i in 1:n){
      vaf.fit[i] <- 1/ (1 + ( 1-purity[i] + purity[i]*CNA[i]*Mu[cluster[i]])/((1-purity[i])*Alpha[cluster[i]] + purity[i]*CNB[i]*Omega[cluster[i]]*Mu[cluster[i]]))
    }
    obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }else{
    tmp.res<-Ares_mix_SNP_Mstep_Newton(matrix(1,length(B),1), B, R, purity, CNA, CNB, K.cluster)
    Alpha<-tmp.res$Alpha
    Mu<-tmp.res$Mu
    Omega<-tmp.res$Omega
    obj<-Ares.mix.SNP.obj(K.cluster, 1, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    
    cluster<-rep(1,n)
    vaf.fit<-1/ (1 + ( 1-purity + purity*CNA*Mu)/((1-purity)*Alpha + purity*CNB*Omega*Mu))
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega=Omega, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.SNP.Newton.rep<-function(B, R, purity, CNA, CNB, K.cluster, replicates=10, cores=1){
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.SNP.Newton, B=B, R=R, purity=purity, CNA=CNA, CNB=CNB, K.cluster=K.cluster, mc.cores=cores)
    
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
    res <- tryCatch(Finite.Binomial.Ares.mix.SNP.Newton(X=1, B, R, purity, CNA, CNB, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  return(res)
}

Finite.Binomial.Ares.mix.SNP.Gradient<-function(X, B, R, purity, CNA, CNB, K.cluster){
  n<-length(R)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Alpha <- exp(runif(K.cluster, min=-log(5), max = log(5)))
    Omega <- exp(runif(K.cluster, min=-log(20), max = log(20)))
    Mu <- exp(runif(K.cluster, min=-log(20), max = log(20)))
    # Alpha <- exp(seq(-log(5), log(5), length.out = K.cluster))
    # Omega <- exp(seq(-log(20), log(20), length.out = K.cluster))
    # Mu <- exp(seq(-log(20), log(20), length.out = K.cluster))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    while(err>1e-4 & iter<100 & obj>0 & nondecrease){
      Pi_old<-Pi
      Alpha_old<-Alpha
      Mu_old<-Mu
      Omega_old<-Omega
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      ### M step
      tmp.res<-Ares_mix_SNP_Gradient(Pi, B, R, purity, CNA, CNB, K.cluster)
      Alpha<-tmp.res$Alpha
      Mu<-tmp.res$Mu
      Omega<-tmp.res$Omega
      obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Alpha<-Alpha_old
      Mu<-Mu_old
      Omega<-Omega_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    vaf.fit<- 1/ (1 + ( 1-purity + purity*CNA*Mu[cluster])/((1-purity)*Alpha[cluster] + purity*CNB*Omega[cluster]*Mu[cluster]))
    obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }else{
    tmp.res<-Ares_mix_SNP_Mstep_Newton(matrix(1,length(B),1), B, R, purity, CNA, CNB, K.cluster)
    Alpha<-tmp.res$Alpha
    Mu<-tmp.res$Mu
    Omega<-tmp.res$Omega
    obj<-Ares.mix.SNP.obj(K.cluster, 1, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    
    cluster<-rep(1,n)
    vaf.fit<-1/ (1 + ( 1-purity + purity*CNA*Mu)/((1-purity)*Alpha + purity*CNB*Omega*Mu))
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega=Omega, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.SNP.Gradient.rep<-function(B, R, purity, CNA, CNB, K.cluster, replicates=10, cores=1){
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.SNP.Gradient, B=B, R=R, purity=purity, CNA=CNA, CNB=CNB, K.cluster=K.cluster, mc.cores=cores)
    
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
    res <- tryCatch(Finite.Binomial.Ares.mix.SNP.Gradient(X=1, B, R, purity, CNA, CNB, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  return(res)
}



Finite.Binomial.Ares.mix.SNP.C.Stochastic<-function(X, B, R, purity, CNA, CNB, K.cluster, srate){
  n<-length(R)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Alpha <- exp(runif(K.cluster, min=-log(5), max = log(5)))
    Omega <- exp(runif(K.cluster, min=-log(20), max = log(20)))
    Mu <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    # Alpha <- exp(seq(-log(5), log(5), length.out = K.cluster))
    # Omega <- exp(seq(-log(20), log(20), length.out = K.cluster))
    # Mu <- exp(seq(-log(20), log(20), length.out = K.cluster))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    
    obj_old<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    while(err>1e-4 & iter<50 & obj>0 & nondecrease){
      Pi_old<-Pi
      Alpha_old<-Alpha
      Mu_old<-Mu
      Omega_old<-Omega
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      ### random subset
      random.id<-sample(x=c(1:n), size=floor(srate*n))
      B.random<-B[random.id]
      R.random<-R[random.id]
      purity.random<-purity[random.id]
      CNA.random<-CNA[random.id]
      CNB.random<-CNB[random.id]
      ### M step
      tmp.res<-Ares_mix_SNP_Mstep(Pi[random.id,], B.random, R.random, purity.random, CNA.random, CNB.random, K.cluster)
      Alpha<-tmp.res$Alpha
      Mu<-tmp.res$Mu
      Omega<-tmp.res$Omega
      obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Alpha<-Alpha_old
      Mu<-Mu_old
      Omega<-Omega_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    vaf.fit<-rep(0,n)
    for(i in 1:n){
      vaf.fit[i] <- 1/ (1 + ( 1-purity[i] + purity[i]*CNA[i]*Mu[cluster[i]])/((1-purity[i])*Alpha[cluster[i]] + purity[i]*CNB[i]*Omega[cluster[i]]*Mu[cluster[i]]))
    }
    obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }else{
    tmp.res<-Solve_Ares_mix_SNP(matrix(1,length(B),1), B, R, purity, CNA, CNB)
    Alpha<-tmp.res$alpha
    Mu<-tmp.res$mu
    Omega<-tmp.res$omega
    obj<-Ares.mix.SNP.obj(K.cluster, 1, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    
    cluster<-rep(1,n)
    vaf.fit<-1/ (1 + ( 1-purity + purity*CNA*Mu)/((1-purity)*Alpha + purity*CNB*Omega*Mu))
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega=Omega, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.SNP.C.Stochastic.rep<-function(B, R, purity, CNA, CNB, K.cluster, srate, replicates=10, cores=1){
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.SNP.C.Stochastic, B=B, R=R, purity=purity, CNA=CNA, CNB=CNB, K.cluster=K.cluster, srate=srate, mc.cores=cores)
    
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
    res <- tryCatch(Finite.Binomial.Ares.mix.SNP.C.Stochastic(X=1, B, R, purity, CNA, CNB, K.cluster, srate), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  return(res)
}


Finite.Binomial.Ares.mix.SNP.Newton.Stochastic<-function(X, B, R, purity, CNA, CNB, K.cluster, srate){
  n<-length(R)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Alpha <- exp(runif(K.cluster, min=-log(5), max = log(5)))
    Omega <- exp(runif(K.cluster, min=-log(20), max = log(20)))
    Mu <- exp(runif(K.cluster, min=-log(20), max = log(20)))
    # Alpha <- exp(seq(-log(5), log(5), length.out = K.cluster))
    # Omega <- exp(seq(-log(20), log(20), length.out = K.cluster))
    # Mu <- exp(seq(-log(20), log(20), length.out = K.cluster))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    
    obj_old<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    while(err>1e-4 & iter<100 & obj>0 & nondecrease){
      Pi_old<-Pi
      Alpha_old<-Alpha
      Mu_old<-Mu
      Omega_old<-Omega
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      ### random subset
      random.id<-sample(x=c(1:n), size=floor(srate*n))
      B.random<-B[random.id]
      R.random<-R[random.id]
      purity.random<-purity[random.id]
      CNA.random<-CNA[random.id]
      CNB.random<-CNB[random.id]
      ### M step
      tmp.res<-Ares_mix_SNP_Mstep_Newton(Pi[random.id,], B.random, R.random, purity.random, CNA.random, CNB.random, K.cluster)
      Alpha<-tmp.res$Alpha
      Mu<-tmp.res$Mu
      Omega<-tmp.res$Omega
      obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Alpha<-Alpha_old
      Mu<-Mu_old
      Omega<-Omega_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    vaf.fit<-rep(0,n)
    for(i in 1:n){
      vaf.fit[i] <- 1/ (1 + ( 1-purity[i] + purity[i]*CNA[i]*Mu[cluster[i]])/((1-purity[i])*Alpha[cluster[i]] + purity[i]*CNB[i]*Omega[cluster[i]]*Mu[cluster[i]]))
    }
    obj<-Ares.mix.SNP.obj(K.cluster, Pi.prior, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }else{
    tmp.res<-Ares_mix_SNP_Mstep_Newton(matrix(1,length(B),1), B, R, purity, CNA, CNB, K.cluster)
    Alpha<-tmp.res$Alpha
    Mu<-tmp.res$Mu
    Omega<-tmp.res$Omega
    obj<-Ares.mix.SNP.obj(K.cluster, 1, B, R, purity, CNA, CNB, Alpha, Omega, Mu)
    
    cluster<-rep(1,n)
    vaf.fit<-1/ (1 + ( 1-purity + purity*CNA*Mu)/((1-purity)*Alpha + purity*CNB*Omega*Mu))
    BIC<-log(n)*3*K.cluster+2*obj
    AIC<-6*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega=Omega, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.SNP.Newton.Stochastic.rep<-function(B, R, purity, CNA, CNB, K.cluster, srate, replicates=10, cores=1){
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.SNP.Newton.Stochastic, B=B, R=R, purity=purity, CNA=CNA, CNB=CNB, K.cluster=K.cluster, srate=srate, mc.cores=cores)
    
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
    res <- tryCatch(Finite.Binomial.Ares.mix.SNP.Newton.Stochastic(X=1, B, R, purity, CNA, CNB, K.cluster, srate), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  return(res)
}


### Combined model
Binomial.Likelihood.Combined<-function(B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, alpha, omega.SNP, omega.SNV, mu){
  theta.SNP <- 1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*mu)/((1-purity.SNP)*alpha + purity.SNP*CNB*omega.SNP*mu))
  theta.SNV = purity.SNV*CNv*omega.SNV/(2.0*(1.0-purity.SNV)/mu+purity.SNV*(CNv*omega.SNV+CNr))
  res<-c(dbinom(x=B, size=R.SNP, prob=theta.SNP), dbinom(x=V, size=R.SNV, prob=theta.SNV))
  return(res)
}

Ares.mix.Combined.obj<-function(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu){
  n<-length(R.SNP)+length(R.SNV)
  OBJ<-matrix(0,n,K.cluster)
  for(k in 1:K.cluster){
    OBJ[,k]<-Pi.prior[k]*Binomial.Likelihood.Combined(B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha[k], omega.SNP=Omega.SNP[k], omega.SNV=Omega.SNV[k], mu=Mu[k])
  }
  res<-(-sum(log(rowSums(OBJ)+1e-250)))
  return(res)
}

Finite.Binomial.Ares.mix.Combined.Mu<-function(X, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster){
  n<-length(R.SNP)+length(R.SNV)
  n.SNP<-length(R.SNP)
  n.SNV<-length(R.SNV)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega.SNP<-rep(1,K.cluster)
  Omega.SNV<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Mu <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    # Alpha <- exp(seq(-log(5), log(5), length.out = K.cluster))
    # Omega <- exp(seq(-log(20), log(20), length.out = K.cluster))
    # Mu <- exp(seq(-log(20), log(20), length.out = K.cluster))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
    while(err>1e-4 & iter<50 & obj>0 & nondecrease){
      Pi_old<-Pi
      Alpha_old<-Alpha
      Mu_old<-Mu
      Omega.SNP_old<-Omega.SNP
      Omega.SNV_old<-Omega.SNV
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi.SNP<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu)
      Pi.SNV<-Posterior.Cluster.Finite.Binomial.true(K.cluster, Pi.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu)
      Pi<-rbind(Pi.SNP, Pi.SNV)
      ### M step
      tmp.res<-Ares_mix_Combined_Mstep_Mu(Pi.SNP, B, R.SNP, purity.SNP, CNA, CNB, Pi.SNV, V, R.SNV, purity.SNV, CNv, CNr, K.cluster)
      
      Alpha<-tmp.res$Alpha
      Mu<-tmp.res$Mu
      Omega.SNP<-tmp.res$Omega_SNP
      Omega.SNV<-tmp.res$Omega_SNV
      obj<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Alpha<-Alpha_old
      Mu<-Mu_old
      Omega.SNP<-Omega.SNP_old
      Omega.SNV<-Omega.SNV_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    
    cluster1<-cluster[1:n.SNP]
    cluster2<-cluster[(n.SNP+1):n]
    vaf.fit1 <- 1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*Mu[cluster1])/((1-purity.SNP)*Alpha[cluster1] + purity.SNP*CNB*Omega.SNP[cluster1]*Mu[cluster1]))
    vaf.fit2 <- purity.SNV*CNv*Omega.SNV[cluster2]/(2*(1-purity.SNV)/Mu[cluster2]+purity.SNV*(CNv*Omega.SNV[cluster2]+CNr))
    vaf.fit<-c(vaf.fit1, vaf.fit2)
    
    # vaf.fit<-rep(0,n)
    # for(i in 1:n.SNP){
    #   vaf.fit[i] <- 1/ (1 + ( 1-purity.SNP[i] + purity.SNP[i]*CNA[i]*Mu[cluster[i]])/((1-purity.SNP[i])*Alpha[cluster[i]] + purity.SNP[i]*CNB[i]*Omega.SNP[cluster[i]]*Mu[cluster[i]]))
    # }
    # for(i in 1:n.SNV){
    #   vaf.fit[i+n.SNP] <- purity.SNV[i]*CNv[i]*Omega.SNV[cluster[i+n.SNP]]/(2*(1-purity.SNV[i])/Mu[cluster[i+n.SNP]]+purity.SNV[i]*(CNv[i]*Omega.SNV[cluster[i+n.SNP]]+CNr[i]))
    # }
    obj<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
    BIC<-log(n)*K.cluster+2*obj
    AIC<-2*K.cluster+2*obj
  }else{
    tmp.res<-Solve_Ares_mix_Combined_Mu(matrix(1,length(B),1), B, R.SNP, purity.SNP, CNA, CNB, matrix(1,length(V),1), V, R.SNV, purity.SNV, CNv, CNr)
    Alpha<-tmp.res$Alpha
    Mu<-tmp.res$Mu
    Omega.SNP<-tmp.res$Omega_SNP
    Omega.SNV<-tmp.res$Omega_SNV
    obj<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
    
    cluster<-rep(1,n)
    vaf.fit.SNP<-1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*Mu)/((1-purity.SNP)*Alpha + purity.SNP*CNB*Omega.SNP*Mu))
    vaf.fit.SNV<-purity.SNV*CNv*Omega.SNV/(2*(1-purity.SNV)/Mu+purity.SNV*(CNv*Omega.SNV+CNr))
    vaf.fit<-c(vaf.fit.SNP, vaf.fit.SNV)
    BIC<-log(n)*K.cluster+2*obj
    AIC<-2*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega.SNP=Omega.SNP, Omega.SNV=Omega.SNV, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.Combined.Mu.rep<-function(B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster, replicates=10, cores=1){
  Mutation<-c(rep("SNP", length(B)), rep("SNV", length(V)))
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.Combined.Mu, B=B, R.SNP=R.SNP, purity.SNP=purity.SNP, CNA=CNA, CNB=CNB, V=V, R.SNV=R.SNV, 
                              purity.SNV=purity.SNV, CNv=CNv, CNr=CNr, K.cluster=K.cluster, mc.cores=cores)
    
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
    New.cluster.order<-order(res$Mu)
    Alpha.new<-res$Alpha[New.cluster.order]
    Omega.SNP.new<-res$Omega.SNP[New.cluster.order]
    Omega.SNV.new<-res$Omega.SNV[New.cluster.order]
    Mu.new<-res$Mu[New.cluster.order]
    Pi.prior.new<-res$Pi.prior[New.cluster.order]
    New.cluster.order.match<-order(New.cluster.order)
    cluster.new<-New.cluster.order.match[res$cluster]
    res$Alpha<-Alpha.new
    res$Omega.SNP<-Omega.SNP.new
    res$Omega.SNV<-Omega.SNV.new
    res$Mu<-Mu.new
    res$Pi.prior<-Pi.prior.new
    res$cluster<-cluster.new
    res$obj.seq<-obj.tmp
    res$BIC.seq<-BIC.tmp
    res$AIC.seq<-AIC.tmp
  }else{
    res <- tryCatch(Finite.Binomial.Ares.mix.Combined.Mu(X=1, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  res$Mutation<-Mutation
  return(res)
}


Finite.Binomial.Ares.mix.Combined<-function(X, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster){
  n<-length(R.SNP)+length(R.SNV)
  n.SNP<-length(R.SNP)
  n.SNV<-length(R.SNV)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega.SNP<-rep(1,K.cluster)
  Omega.SNV<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Mu <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    Alpha <- exp(runif(K.cluster, min=-log(5), max = log(5)))
    Omega.SNP <- exp(runif(K.cluster, min=-log(10), max = log(10)))
    Omega.SNV <- exp(runif(K.cluster, min=-log(10), max = log(10)))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
    while(err>1e-4 & iter<50 & obj>0 & nondecrease){
      Pi_old<-Pi
      Alpha_old<-Alpha
      Mu_old<-Mu
      Omega.SNP_old<-Omega.SNP
      Omega.SNV_old<-Omega.SNV
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi.SNP<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu)
      Pi.SNV<-Posterior.Cluster.Finite.Binomial.true(K.cluster, Pi.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu)
      Pi<-rbind(Pi.SNP, Pi.SNV)
      ### M step
      tmp.res<-Ares_mix_Combined_Mstep(Pi.SNP, B, R.SNP, purity.SNP, CNA, CNB, Pi.SNV, V, R.SNV, purity.SNV, CNv, CNr, K.cluster)
      
      Alpha<-tmp.res$Alpha
      Mu<-tmp.res$Mu
      Omega.SNP<-tmp.res$Omega_SNP
      Omega.SNV<-tmp.res$Omega_SNV
      obj<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Alpha<-Alpha_old
      Mu<-Mu_old
      Omega.SNP<-Omega.SNP_old
      Omega.SNV<-Omega.SNV_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    
    cluster1<-cluster[1:n.SNP]
    cluster2<-cluster[(n.SNP+1):n]
    vaf.fit1 <- 1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*Mu[cluster1])/((1-purity.SNP)*Alpha[cluster1] + purity.SNP*CNB*Omega.SNP[cluster1]*Mu[cluster1]))
    vaf.fit2 <- purity.SNV*CNv*Omega.SNV[cluster2]/(2*(1-purity.SNV)/Mu[cluster2]+purity.SNV*(CNv*Omega.SNV[cluster2]+CNr))
    vaf.fit<-c(vaf.fit1, vaf.fit2)
    
    # vaf.fit<-rep(0,n)
    # for(i in 1:n.SNP){
    #   vaf.fit[i] <- 1/ (1 + ( 1-purity.SNP[i] + purity.SNP[i]*CNA[i]*Mu[cluster[i]])/((1-purity.SNP[i])*Alpha[cluster[i]] + purity.SNP[i]*CNB[i]*Omega.SNP[cluster[i]]*Mu[cluster[i]]))
    # }
    # for(i in 1:n.SNV){
    #   vaf.fit[i+n.SNP] <- purity.SNV[i]*CNv[i]*Omega.SNV[cluster[i+n.SNP]]/(2*(1-purity.SNV[i])/Mu[cluster[i+n.SNP]]+purity.SNV[i]*(CNv[i]*Omega.SNV[cluster[i+n.SNP]]+CNr[i]))
    # }
    obj<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
    BIC<-log(n)+log(n.SNV)*2*(K.cluster-1) +log(n.SNP)*3*(K.cluster-1) +2*obj
    AIC<-6*K.cluster+2*obj
  }else{
    tmp.res<-Solve_Ares_mix_Combined(matrix(1,length(B),1), B, R.SNP, purity.SNP, CNA, CNB, matrix(1,length(V),1), V, R.SNV, purity.SNV, CNv, CNr)
    Alpha<-tmp.res$Alpha
    Mu<-tmp.res$Mu
    Omega.SNP<-tmp.res$Omega_SNP
    Omega.SNV<-tmp.res$Omega_SNV
    obj<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
    
    cluster<-rep(1,n)
    vaf.fit.SNP<-1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*Mu)/((1-purity.SNP)*Alpha + purity.SNP*CNB*Omega.SNP*Mu))
    vaf.fit.SNV<-purity.SNV*CNv*Omega.SNV/(2*(1-purity.SNV)/Mu+purity.SNV*(CNv*Omega.SNV+CNr))
    vaf.fit<-c(vaf.fit.SNP, vaf.fit.SNV)
    BIC<-log(n)+2*obj
    AIC<-6*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega.SNP=Omega.SNP, Omega.SNV=Omega.SNV, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.Combined.rep<-function(B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster, replicates=10, cores=1){
  Mutation<-c(rep("SNP", length(B)), rep("SNV", length(V)))
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.Combined, B=B, R.SNP=R.SNP, purity.SNP=purity.SNP, CNA=CNA, CNB=CNB, V=V, R.SNV=R.SNV, 
                              purity.SNV=purity.SNV, CNv=CNv, CNr=CNr, K.cluster=K.cluster, mc.cores=cores)
    
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
    New.cluster.order<-order(res$Mu)
    Alpha.new<-res$Alpha[New.cluster.order]
    Omega.SNP.new<-res$Omega.SNP[New.cluster.order]
    Omega.SNV.new<-res$Omega.SNV[New.cluster.order]
    Mu.new<-res$Mu[New.cluster.order]
    Pi.prior.new<-res$Pi.prior[New.cluster.order]
    New.cluster.order.match<-order(New.cluster.order)
    cluster.new<-New.cluster.order.match[res$cluster]
    res$Alpha<-Alpha.new
    res$Omega.SNP<-Omega.SNP.new
    res$Omega.SNV<-Omega.SNV.new
    res$Mu<-Mu.new
    res$Pi.prior<-Pi.prior.new
    res$cluster<-cluster.new
    res$obj.seq<-obj.tmp
    res$BIC.seq<-BIC.tmp
    res$AIC.seq<-AIC.tmp
  }else{
    res <- tryCatch(Finite.Binomial.Ares.mix.Combined(X=1, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  res$Mutation<-Mutation
  return(res)
}

Finite.Binomial.Ares.mix.Combined2<-function(X, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster){
  n<-length(R.SNP)+length(R.SNV)
  n.SNP<-length(R.SNP)
  n.SNV<-length(R.SNV)
  Pi.prior<-rep(1/K.cluster,K.cluster)
  Pi<-matrix(1/K.cluster,n,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu<-rep(1,K.cluster)
  Omega.SNP<-rep(1,K.cluster)
  Omega.SNV<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Mu <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    Alpha <- exp(runif(K.cluster, min=-log(5), max = log(5)))
    Omega.SNP <- exp(runif(K.cluster, min=-log(10), max = log(10)))
    Omega.SNV <- exp(runif(K.cluster, min=-log(10), max = log(10)))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
    while(err>1e-4 & iter<50 & obj>0 & nondecrease){
      Pi_old<-Pi
      Alpha_old<-Alpha
      Mu_old<-Mu
      Omega.SNP_old<-Omega.SNP
      Omega.SNV_old<-Omega.SNV
      Pi.prior_old<-Pi.prior
      ### E step
      Pi.prior<-colMeans(Pi)
      Pi.SNP<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu)
      Pi.SNV<-Posterior.Cluster.Finite.Binomial.true(K.cluster, Pi.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu)
      Pi<-rbind(Pi.SNP, Pi.SNV)
      ### M step
      tmp.res<-Ares_mix_Combined_Mstep2(Pi.SNP, B, R.SNP, purity.SNP, CNA, CNB, Pi.SNV, V, R.SNV, purity.SNV, CNv, CNr, K.cluster)
      
      Alpha<-tmp.res$Alpha
      Mu<-tmp.res$Mu
      Omega.SNP<-tmp.res$Omega_SNP
      Omega.SNV<-tmp.res$Omega_SNV
      obj<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi<-Pi_old
      Alpha<-Alpha_old
      Mu<-Mu_old
      Omega.SNP<-Omega.SNP_old
      Omega.SNV<-Omega.SNV_old
      Pi.prior<-Pi.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    cluster<-apply(Pi,1,which.max)
    
    cluster1<-cluster[1:n.SNP]
    cluster2<-cluster[(n.SNP+1):n]
    vaf.fit1 <- 1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*Mu[cluster1])/((1-purity.SNP)*Alpha[cluster1] + purity.SNP*CNB*Omega.SNP[cluster1]*Mu[cluster1]))
    vaf.fit2 <- purity.SNV*CNv*Omega.SNV[cluster2]/(2*(1-purity.SNV)/Mu[cluster2]+purity.SNV*(CNv*Omega.SNV[cluster2]+CNr))
    vaf.fit<-c(vaf.fit1, vaf.fit2)
    
    # vaf.fit<-rep(0,n)
    # for(i in 1:n.SNP){
    #   vaf.fit[i] <- 1/ (1 + ( 1-purity.SNP[i] + purity.SNP[i]*CNA[i]*Mu[cluster[i]])/((1-purity.SNP[i])*Alpha[cluster[i]] + purity.SNP[i]*CNB[i]*Omega.SNP[cluster[i]]*Mu[cluster[i]]))
    # }
    # for(i in 1:n.SNV){
    #   vaf.fit[i+n.SNP] <- purity.SNV[i]*CNv[i]*Omega.SNV[cluster[i+n.SNP]]/(2*(1-purity.SNV[i])/Mu[cluster[i+n.SNP]]+purity.SNV[i]*(CNv[i]*Omega.SNV[cluster[i+n.SNP]]+CNr[i]))
    # }
    obj<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
    BIC<-log(n)+log(n.SNV)*2*(K.cluster-1) +log(n.SNP)*3*(K.cluster-1) +2*obj
    AIC<-6*K.cluster+2*obj
  }else{
    tmp.res<-Solve_Ares_mix_Combined(matrix(1,length(B),1), B, R.SNP, purity.SNP, CNA, CNB, matrix(1,length(V),1), V, R.SNV, purity.SNV, CNv, CNr)
    Alpha<-tmp.res$Alpha
    Mu<-tmp.res$Mu
    Omega.SNP<-tmp.res$Omega_SNP
    Omega.SNV<-tmp.res$Omega_SNV
    obj<-Ares.mix.Combined.obj(K.cluster, Pi.prior, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, Alpha, Omega.SNP, Omega.SNV, Mu)
    
    cluster<-rep(1,n)
    vaf.fit.SNP<-1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*Mu)/((1-purity.SNP)*Alpha + purity.SNP*CNB*Omega.SNP*Mu))
    vaf.fit.SNV<-purity.SNV*CNv*Omega.SNV/(2*(1-purity.SNV)/Mu+purity.SNV*(CNv*Omega.SNV+CNr))
    vaf.fit<-c(vaf.fit.SNP, vaf.fit.SNV)
    BIC<-log(n)+2*obj
    AIC<-6*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega.SNP=Omega.SNP, Omega.SNV=Omega.SNV, Mu=Mu, obj=obj, Pi.prior=Pi.prior, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.Combined2.rep<-function(B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster, replicates=10, cores=1){
  Mutation<-c(rep("SNP", length(B)), rep("SNV", length(V)))
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.Combined2, B=B, R.SNP=R.SNP, purity.SNP=purity.SNP, CNA=CNA, CNB=CNB, V=V, R.SNV=R.SNV, 
                              purity.SNV=purity.SNV, CNv=CNv, CNr=CNr, K.cluster=K.cluster, mc.cores=cores)
    
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
    New.cluster.order<-order(res$Mu)
    Alpha.new<-res$Alpha[New.cluster.order]
    Omega.SNP.new<-res$Omega.SNP[New.cluster.order]
    Omega.SNV.new<-res$Omega.SNV[New.cluster.order]
    Mu.new<-res$Mu[New.cluster.order]
    Pi.prior.new<-res$Pi.prior[New.cluster.order]
    New.cluster.order.match<-order(New.cluster.order)
    cluster.new<-New.cluster.order.match[res$cluster]
    res$Alpha<-Alpha.new
    res$Omega.SNP<-Omega.SNP.new
    res$Omega.SNV<-Omega.SNV.new
    res$Mu<-Mu.new
    res$Pi.prior<-Pi.prior.new
    res$cluster<-cluster.new
    res$obj.seq<-obj.tmp
    res$BIC.seq<-BIC.tmp
    res$AIC.seq<-AIC.tmp
  }else{
    res <- tryCatch(Finite.Binomial.Ares.mix.Combined2(X=1, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  res$Mutation<-Mutation
  return(res)
}

Finite.Binomial.Ares.mix.Combined.Independent<-function(X, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster){
  n<-length(R.SNP)+length(R.SNV)
  n.SNP<-length(R.SNP)
  n.SNV<-length(R.SNV)
  Pi.SNP.prior<-rep(1/K.cluster,K.cluster)
  Pi.SNV.prior<-rep(1/K.cluster,K.cluster)
  Pi.SNP<-matrix(1/K.cluster,n.SNP,K.cluster)
  Pi.SNV<-matrix(1/K.cluster,n.SNV,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu.SNP<-rep(1,K.cluster)
  Mu.SNV<-rep(1,K.cluster)
  Omega.SNP<-rep(1,K.cluster)
  Omega.SNV<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Mu.SNP <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    Mu.SNV <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    Alpha <- exp(runif(K.cluster, min=-log(5), max = log(5)))
    Omega.SNP <- exp(runif(K.cluster, min=-log(10), max = log(10)))
    Omega.SNV <- exp(runif(K.cluster, min=-log(10), max = log(10)))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Ares.mix.SNP.obj(K.cluster, Pi.SNP.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu.SNP) + Finite.Binomial.obj.repara(K.cluster, Pi.SNV.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu.SNV)
    while(err>1e-4 & iter<50 & obj>0 & nondecrease){
      Pi.SNP_old<-Pi.SNP
      Pi.SNV_old<-Pi.SNV
      Alpha_old<-Alpha
      Mu.SNP_old<-Mu.SNP
      Mu.SNV_old<-Mu.SNV
      Omega.SNP_old<-Omega.SNP
      Omega.SNV_old<-Omega.SNV
      Pi.SNP.prior_old<-Pi.SNP.prior
      Pi.SNV.prior_old<-Pi.SNV.prior
      ### E step
      Pi.SNP.prior<-colMeans(Pi.SNP)
      Pi.SNV.prior<-colMeans(Pi.SNV)
      Pi.SNP<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.SNP.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu.SNP)
      Pi.SNV<-Posterior.Cluster.Finite.Binomial.true(K.cluster, Pi.SNV.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu.SNV)
      ### M step
      tmp.res<-Ares_mix_Combined_Independent(Pi.SNP, B, R.SNP, purity.SNP, CNA, CNB, Pi.SNV, V, R.SNV, purity.SNV, CNv, CNr, K.cluster)
      
      Alpha<-tmp.res$Alpha
      Mu.SNP<-tmp.res$Mu_SNP
      Mu.SNV<-tmp.res$Mu_SNV
      Omega.SNP<-tmp.res$Omega_SNP
      Omega.SNV<-tmp.res$Omega_SNV
      obj<-Ares.mix.SNP.obj(K.cluster, Pi.SNP.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu.SNP) + Finite.Binomial.obj.repara(K.cluster, Pi.SNV.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu.SNV)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi.SNP<-Pi.SNP_old
      Pi.SNV<-Pi.SNV_old
      Alpha<-Alpha_old
      Mu.SNP<-Mu.SNP_old
      Mu.SNV<-Mu.SNV_old
      Omega.SNP<-Omega.SNP_old
      Omega.SNV<-Omega.SNV_old
      Pi.SNP.prior<-Pi.SNP.prior_old
      Pi.SNV.prior<-Pi.SNV.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    
    Pi<-rbind(Pi.SNP, Pi.SNV)
    cluster<-apply(Pi,1,which.max)
    cluster1<-cluster[1:n.SNP]
    cluster2<-cluster[(n.SNP+1):n]
    vaf.fit1 <- 1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*Mu.SNP[cluster1])/((1-purity.SNP)*Alpha[cluster1] + purity.SNP*CNB*Omega.SNP[cluster1]*Mu.SNP[cluster1]))
    vaf.fit2 <- purity.SNV*CNv*Omega.SNV[cluster2]/(2*(1-purity.SNV)/Mu.SNV[cluster2]+purity.SNV*(CNv*Omega.SNV[cluster2]+CNr))
    vaf.fit<-c(vaf.fit1, vaf.fit2)
    
    obj<-Ares.mix.SNP.obj(K.cluster, Pi.SNP.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu.SNP) + Finite.Binomial.obj.repara(K.cluster, Pi.SNV.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu.SNV)
    BIC<-log(n.SNV)+log(n.SNP)+log(n.SNV)*2*(K.cluster-1) +log(n.SNP)*3*(K.cluster-1) +2*obj
    AIC<-6*K.cluster+2*obj
  }else{
    tmp.res<-Solve_Ares_mix_Combined_Independent(matrix(1,length(B),1), B, R.SNP, purity.SNP, CNA, CNB, matrix(1,length(V),1), V, R.SNV, purity.SNV, CNv, CNr)
    Alpha<-tmp.res$Alpha
    Mu.SNP<-tmp.res$Mu_SNP
    Mu.SNV<-tmp.res$Mu_SNV
    Omega.SNP<-tmp.res$Omega_SNP
    Omega.SNV<-tmp.res$Omega_SNV
    obj<-Ares.mix.SNP.obj(K.cluster, Pi.SNP.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu.SNP) + Finite.Binomial.obj.repara(K.cluster, Pi.SNV.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu.SNV)
    cluster<-rep(1,n)
    vaf.fit.SNP<-1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*Mu.SNP)/((1-purity.SNP)*Alpha + purity.SNP*CNB*Omega.SNP*Mu.SNP))
    vaf.fit.SNV<-purity.SNV*CNv*Omega.SNV/(2*(1-purity.SNV)/Mu.SNV+purity.SNV*(CNv*Omega.SNV+CNr))
    vaf.fit<-c(vaf.fit.SNP, vaf.fit.SNV)
    BIC<-log(n.SNV)+log(n.SNP)+2*obj
    AIC<-6*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega.SNP=Omega.SNP, Omega.SNV=Omega.SNV, Mu.SNP=Mu.SNP, Mu.SNV=Mu.SNV, obj=obj, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.Combined.Independent.rep<-function(B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster, replicates=10, cores=1){
  Mutation<-c(rep("SNP", length(B)), rep("SNV", length(V)))
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.Combined.Independent, B=B, R.SNP=R.SNP, purity.SNP=purity.SNP, CNA=CNA, CNB=CNB, V=V, R.SNV=R.SNV, 
                              purity.SNV=purity.SNV, CNv=CNv, CNr=CNr, K.cluster=K.cluster, mc.cores=cores)
    
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
    New.cluster.order<-order(res$Mu.SNV)
    Alpha.new<-res$Alpha[New.cluster.order]
    Omega.SNP.new<-res$Omega.SNP[New.cluster.order]
    Omega.SNV.new<-res$Omega.SNV[New.cluster.order]
    Mu.SNP.new<-res$Mu.SNP[New.cluster.order]
    Mu.SNV.new<-res$Mu.SNV[New.cluster.order]
    Pi.prior.new<-res$Pi.prior[New.cluster.order]
    New.cluster.order.match<-order(New.cluster.order)
    cluster.new<-New.cluster.order.match[res$cluster]
    res$Alpha<-Alpha.new
    res$Omega.SNP<-Omega.SNP.new
    res$Omega.SNV<-Omega.SNV.new
    res$Mu.SNP<-Mu.SNP.new
    res$Mu.SNV<-Mu.SNV.new
    res$cluster<-cluster.new
    res$obj.seq<-obj.tmp
    res$BIC.seq<-BIC.tmp
    res$AIC.seq<-AIC.tmp
  }else{
    res <- tryCatch(Finite.Binomial.Ares.mix.Combined.Independent(X=1, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  res$Mutation<-Mutation
  return(res)
}

Finite.Binomial.Ares.mix.Joint<-function(X, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster){
  n<-length(R.SNP)+length(R.SNV)
  n.SNP<-length(R.SNP)
  n.SNV<-length(R.SNV)
  Pi.SNP.prior<-rep(1/K.cluster,K.cluster)
  Pi.SNV.prior<-rep(1/K.cluster,K.cluster)
  Pi.SNP<-matrix(1/K.cluster,n.SNP,K.cluster)
  Pi.SNV<-matrix(1/K.cluster,n.SNV,K.cluster)
  Alpha<-rep(1,K.cluster)
  Mu.SNP<-rep(1,K.cluster)
  Mu.SNV<-rep(1,K.cluster)
  Omega.SNP<-rep(1,K.cluster)
  Omega.SNV<-rep(1,K.cluster)
  if(K.cluster>1){
    ### random initial Alpha, Omega and Eta
    Mu.SNP <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    Mu.SNV <- exp(runif(K.cluster, min=-log(10), max = log(40)))
    Alpha <- exp(runif(K.cluster, min=-log(5), max = log(5)))
    Omega.SNP <- exp(runif(K.cluster, min=-log(10), max = log(10)))
    Omega.SNV <- exp(runif(K.cluster, min=-log(10), max = log(10)))
    
    err<-1
    iter<-1
    obj<-1e10
    nondecrease<-TRUE
    obj_old<-Ares.mix.SNP.obj(K.cluster, Pi.SNP.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu.SNP) + Finite.Binomial.obj.repara(K.cluster, Pi.SNV.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu.SNV)
    while(err>1e-4 & iter<50 & obj>0 & nondecrease){
      Pi.SNP_old<-Pi.SNP
      Pi.SNV_old<-Pi.SNV
      Alpha_old<-Alpha
      Mu.SNP_old<-Mu.SNP
      Mu.SNV_old<-Mu.SNV
      Omega.SNP_old<-Omega.SNP
      Omega.SNV_old<-Omega.SNV
      Pi.SNP.prior_old<-Pi.SNP.prior
      Pi.SNV.prior_old<-Pi.SNV.prior
      ### E step
      Pi.SNP.prior<-colMeans(Pi.SNP)
      Pi.SNV.prior<-colMeans(Pi.SNV)
      Pi.SNP<-Posterior.Cluster.Ares.mix.SNP(K.cluster, Pi.SNP.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu.SNP)
      Pi.SNV<-Posterior.Cluster.Finite.Binomial.true(K.cluster, Pi.SNV.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu.SNV)
      ### M step
      tmp.res<-Ares_mix_Joint(Pi.SNP, B, R.SNP, purity.SNP, CNA, CNB, Pi.SNV, V, R.SNV, purity.SNV, CNv, CNr, K.cluster)
      
      Alpha<-tmp.res$Alpha
      Mu.SNP<-tmp.res$Mu_SNP
      Mu.SNV<-tmp.res$Mu_SNV
      Omega.SNP<-tmp.res$Omega_SNP
      Omega.SNV<-tmp.res$Omega_SNV
      obj<-Ares.mix.SNP.obj(K.cluster, Pi.SNP.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu.SNP) + Finite.Binomial.obj.repara(K.cluster, Pi.SNV.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu.SNV)
      err<-abs((obj_old-obj)/obj_old)
      #print(obj)
      obj_old<-obj
      if(is.na((obj_old-obj)/obj_old) | (obj_old-obj)/obj_old<(-0.1)){
        nondecrease=FALSE
      }
      iter=iter+1
    }
    if(obj<0 | !nondecrease){
      obj<-obj_old
      Pi.SNP<-Pi.SNP_old
      Pi.SNV<-Pi.SNV_old
      Alpha<-Alpha_old
      Mu.SNP<-Mu.SNP_old
      Mu.SNV<-Mu.SNV_old
      Omega.SNP<-Omega.SNP_old
      Omega.SNV<-Omega.SNV_old
      Pi.SNP.prior<-Pi.SNP.prior_old
      Pi.SNV.prior<-Pi.SNV.prior_old
    }
    
    if(!nondecrease){
      Converge=FALSE
    }else{
      Converge=TRUE
    }
    Reach.iter=FALSE
    if(iter>199){
      Reach.iter=TRUE
    }
    
    Pi<-rbind(Pi.SNP, Pi.SNV)
    cluster<-apply(Pi,1,which.max)
    cluster1<-cluster[1:n.SNP]
    cluster2<-cluster[(n.SNP+1):n]
    vaf.fit1 <- 1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*Mu.SNP[cluster1])/((1-purity.SNP)*Alpha[cluster1] + purity.SNP*CNB*Omega.SNP[cluster1]*Mu.SNP[cluster1]))
    vaf.fit2 <- purity.SNV*CNv*Omega.SNV[cluster2]/(2*(1-purity.SNV)/Mu.SNV[cluster2]+purity.SNV*(CNv*Omega.SNV[cluster2]+CNr))
    vaf.fit<-c(vaf.fit1, vaf.fit2)
    
    obj<-Ares.mix.SNP.obj(K.cluster, Pi.SNP.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu.SNP) + Finite.Binomial.obj.repara(K.cluster, Pi.SNV.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu.SNV)
    BIC<-log(n.SNV)+log(n.SNP)+log(n.SNV)*2*K.cluster +log(n.SNP)*3*K.cluster +2*obj
    AIC<-6*K.cluster+2*obj
  }else{
    tmp.res<-Solve_Ares_mix_Combined_Independent(matrix(1,length(B),1), B, R.SNP, purity.SNP, CNA, CNB, matrix(1,length(V),1), V, R.SNV, purity.SNV, CNv, CNr)
    Alpha<-tmp.res$Alpha
    Mu.SNP<-tmp.res$Mu_SNP
    Mu.SNV<-tmp.res$Mu_SNV
    Omega.SNP<-tmp.res$Omega_SNP
    Omega.SNV<-tmp.res$Omega_SNV
    obj<-Ares.mix.SNP.obj(K.cluster, Pi.SNP.prior, B, R.SNP, purity.SNP, CNA, CNB, Alpha, Omega.SNP, Mu.SNP) + Finite.Binomial.obj.repara(K.cluster, Pi.SNV.prior, V, R.SNV, purity.SNV, CNv, CNr, Omega.SNV, Mu.SNV)
    cluster<-rep(1,n)
    vaf.fit.SNP<-1/ (1 + ( 1-purity.SNP + purity.SNP*CNA*Mu.SNP)/((1-purity.SNP)*Alpha + purity.SNP*CNB*Omega.SNP*Mu.SNP))
    vaf.fit.SNV<-purity.SNV*CNv*Omega.SNV/(2*(1-purity.SNV)/Mu.SNV+purity.SNV*(CNv*Omega.SNV+CNr))
    vaf.fit<-c(vaf.fit.SNP, vaf.fit.SNV)
    BIC<-log(n.SNV)+log(n.SNP)+2*obj
    AIC<-6*K.cluster+2*obj
  }
  return(list(Alpha=Alpha, Omega.SNP=Omega.SNP, Omega.SNV=Omega.SNV, Mu.SNP=Mu.SNP, Mu.SNV=Mu.SNV, obj=obj, cluster=cluster, vaf.fit=vaf.fit, BIC=BIC, AIC=AIC))
}

Finite.Binomial.Ares.mix.Joint.rep<-function(B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster, replicates=10, cores=1){
  Mutation<-c(rep("SNP", length(B)), rep("SNV", length(V)))
  if(K.cluster>1){
    obj.tmp<-rep(1e20,replicates)
    BIC.tmp<-rep(1e20,replicates)
    AIC.tmp<-rep(1e20,replicates)
    
    res.mixed.model<-mclapply(X=1:replicates, FUN=Finite.Binomial.Ares.mix.Joint, B=B, R.SNP=R.SNP, purity.SNP=purity.SNP, CNA=CNA, CNB=CNB, V=V, R.SNV=R.SNV, 
                              purity.SNV=purity.SNV, CNv=CNv, CNr=CNr, K.cluster=K.cluster, mc.cores=cores)
    
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
    New.cluster.order<-order(res$Mu.SNV)
    Alpha.new<-res$Alpha[New.cluster.order]
    Omega.SNP.new<-res$Omega.SNP[New.cluster.order]
    Omega.SNV.new<-res$Omega.SNV[New.cluster.order]
    Mu.SNP.new<-res$Mu.SNP[New.cluster.order]
    Mu.SNV.new<-res$Mu.SNV[New.cluster.order]
    Pi.prior.new<-res$Pi.prior[New.cluster.order]
    New.cluster.order.match<-order(New.cluster.order)
    cluster.new<-New.cluster.order.match[res$cluster]
    res$Alpha<-Alpha.new
    res$Omega.SNP<-Omega.SNP.new
    res$Omega.SNV<-Omega.SNV.new
    res$Mu.SNP<-Mu.SNP.new
    res$Mu.SNV<-Mu.SNV.new
    res$cluster<-cluster.new
    res$obj.seq<-obj.tmp
    res$BIC.seq<-BIC.tmp
    res$AIC.seq<-AIC.tmp
  }else{
    res <- tryCatch(Finite.Binomial.Ares.mix.Joint(X=1, B, R.SNP, purity.SNP, CNA, CNB, V, R.SNV, purity.SNV, CNv, CNr, K.cluster), error = function(e) {list(0)})
    res$obj.seq<-rep(res$obj, replicates)
    res$BIC.seq<-rep(res$BIC, replicates)
    res$AIC.seq<-rep(res$AIC, replicates)
  }
  res$Mutation<-Mutation
  return(res)
}
