
# Weighted Counting process
Nt.G<-function(Y, delta, case, Gx) # case = 1, case = 2,...
{
  n<-length(Y)
  obs.mat<-matrix(Y,n,n,byrow=T)
  delta.mat<-matrix(delta,n,n,byrow=T)
  Nt.mat<-t(ifelse(obs.mat<=Y & delta.mat == case , 1, 0))  #(ixj) n x n  I(Vi<=Vj, casei=1)

  return(Nt.mat/Gx)
}

# Weighted Counting process
Nt.G.Qt<-function(Y, delta, case, Gx, Qseq, Rtime) # case = 1, case = 2,...
{
  n<-length(Y)
  obs.mat<-matrix(Y,n,n,byrow=T)
  delta.mat<-matrix(delta,n,n,byrow=T)
  Nt.mat<-t(ifelse(obs.mat<=Y & delta.mat == case , 1, 0))  #(ixj) n x n  I(Vi<=Vj, casei=1)
  
  # time-varying weight 
  Qseq.mat<-matrix(Qseq, n, n,byrow=T)
  Rt.mat<-matrix(Rtime, n, n,byrow=T)
  Qt.mat<-t(ifelse(Rt.mat > Y, 1, Qseq.mat))  # Qi(Vj)

  return(list(NtG.m=Nt.mat/Gx, Qt.mat=Qt.mat))
}

# Score function
Ub.sc<-function(beta, X, NtG, Qseq)
{
  np<-length(beta)
  e.xb<-exp(beta%*%t(X))
  P<-1-exp(-e.xb)
  D.eta<-X*c((1-P)*e.xb)
  
  #Un.sc<-colSums(D.eta*c(NtG-P)*Qseq)
  Un.sc<-colSums(D.eta*c(NtG-P)*Qseq*length(Qseq)/sum(Qseq))
  if(sum(is.na(Un.sc))>0){
    Un.sc<-rep(1e-20,np)
  }
  return(Un.sc)
}

estb.Vn.sc<-function(beta.m, X, time, event, weight, Gt, sub.NtG.m, uni.time.e1)
{
  ## for Bn ##
  pts<-time
  n<-dim(X)[1] # = length(time) = n of subj in this regimen = dim(X)[1]
  np<-dim(X)[2]
  npts<-length(pts)
  obs.mat<-matrix(time, npts, n,byrow=T)
  delta.mat<-matrix(event,npts,n,byrow=T)
  nt<-length(uni.time.e1)
  
  Ntt.mat<-t(ifelse(obs.mat>=pts, 1, 0))  #(ixj) n x n  I(Vi>=Vj) not order
  wpi<-colMeans(Ntt.mat*weight) # n-vector
  
  wdNc.mat<-t(ifelse(obs.mat==time & delta.mat == 0, 1, 0))*weight  #(ixj) n x n I(Vi=Vj, casei=0)
  wqi<-colMeans(wdNc.mat) # n-vector
  
  demo.s<-(wpi-wqi)
  demo.s[(demo.s==0)]<-1e-20
  
  wYc.mat<-Ntt.mat*weight 
  whc<-Gt$hx
  wdAc<-t(t(wYc.mat)*whc)
  dMc<-wdNc.mat-wdAc # (ixj) nxn
  
  Gm<-t(t(dMc)/demo.s)
  Ik.mat<-Ntt.mat #t(ifelse(obs.mat>=pts, 1, 0))  #(kxi) n x n  I(Vk>=Vi) 
  
  Ikk<-matrix(rep(t(Ik.mat),n), , npts, byrow=T) 
  dim(Ikk)<-c(n,n,npts) 
  
  Gmm<-array(rep(Gm,each=n),c(n,n,npts))
  Gm.mat<-Gmm*Ikk
  Mc.v<-apply(Gm.mat,1,sum)  # for each Ti
  
  #Indx<-which(event==1)
  
  sd.m<-c()
  V.mm<-c()
  dtau.W.mm<-array(matrix(0, n, nt), c(n, np, nt))
  inv.In.mm<-array(matrix(0, np, np),c(np,  np, nt))
  WI.mm<-array(matrix(0, n, nt),c(n, np, nt))
  
  ind.time<-order(uni.time.e1)
  rnk.time<-rank(uni.time.e1)
  ## sort beta by time
  Beta<-cbind(uni.time.e1,beta.m)[ind.time, ]
  dtau<-c(diff(c(sort(uni.time.e1),max(time))))
  rang<-diff(range(c(uni.time.e1,max(time))))
  avg.beta<-colSums(Beta[,-1]*dtau)/rang
  avg.beta.m<-matrix(avg.beta, nt, np, byrow = T)
  T2.obs<-apply(abs(beta.m-avg.beta.m),2,max)
  
  T3.obs<-colSums((beta.m-avg.beta.m)^2*dtau[rnk.time])
  ## change back
  r.dtau<-dtau[rnk.time]/rang
  
  for (indx in 1:nt)
  {
    beta<-beta.m[indx, ]
    np<-length(beta)
    e.xb<-exp(beta%*%t(X))
    P<-1-exp(-e.xb)
    D.eta<-X*c((1-P)*e.xb)
    D.eta[(e.xb==Inf),]<-0
    
    In.eta<-t(D.eta)%*%(D.eta*weight)
    
    An<-D.eta*c(sub.NtG.m[,indx]-P)*weight
    Bn<--D.eta*c(sub.NtG.m[,indx]*Mc.v)*weight
    
    W1<-An+Bn
    Wn1<-t(W1)%*%W1
    diag.1<-try(solve(In.eta))
    if(class(diag.1)=="try-error"){
      Vn.t<-matrix(0,np,np)
      dtau.W.mm[,,indx]<-matrix(0,n,np)*r.dtau[indx]
      WI.mm[,,indx]<-t(matrix(0,np,np)%*%t(W1))
      #inv.In.mm[,,indx]<-matrix(0, np, np)
      #dtau.m[,,indx]<-matrix(r.dtau[indx],n,np)
      #sd.est<-rep(10000,length(beta))
    }else{
      Vn.t<-solve(In.eta)%*%Wn1%*%solve(In.eta)
      dtau.W.mm[,,indx]<-t(solve(In.eta)%*%t(W1))*r.dtau[indx] #W1*r.dtau[indx]
      WI.mm[,,indx]<-t(solve(In.eta)%*%t(W1))
      #inv.In.mm[,,indx]<-solve(In.eta)
      #dtau.m[,,indx]<-matrix(r.dtau[indx],n,np)
      #sd.est<-sqrt(diag(solve(In.eta)%*%Wn1%*%solve(In.eta)))
    }
    #sd.est<-sqrt(Vn.t)
    sd.m<-rbind(sd.m, sqrt(diag(Vn.t)))
    V.mm<-rbind(V.mm, c(Vn.t))
  }
  
  T1.obs<-apply(abs(beta.m/sd.m),2,max)
  avg.dtau.w<-apply(dtau.W.mm,c(1,2),sum)    
  W.T2.m<-WI.mm-array(avg.dtau.w, c(n, np, nt))
  
  set.seed(15213)
  B<-1000
  sup.T1<-sup.T2<-T3.w<-c()
  for(b in 1:B)
  {
    gi.m<-matrix(rnorm(n, 0, 1), n, np)
    Gi.mm<-array(gi.m, c(n, np, nt))
    
    # FOR T1 #
    WG<-WI.mm*Gi.mm
    T1.w<-abs(t(apply(WG, c(2,3),sum))/sd.m)
    T1.w[is.na(T1.w[,1]),]<-0
    sup.T1<-rbind(sup.T1, apply(T1.w,2,max))
    
    # FOR T2 #
    WG.T2<-W.T2.m*Gi.mm
    T2.w<-abs(t(apply(WG.T2,c(2,3),sum)))
    T2.w[is.na(T2.w[,1]),]<-0
    sup.T2<-rbind(sup.T2, apply(T2.w,2,max))
    
    # FOR T3 #
    WG.T3<-t(apply(WG.T2, c(2,3),sum))^2*r.dtau*rang#*matrix(r.dtau*rang, nt, np)#aperm(array(matrix(r.dtau*rang, nt, np), c(nt, np, n)), c(3,2,1))
    WG.T3[is.na(WG.T3[,1]),]<-0
    T3.w<-rbind(T3.w, apply(WG.T3, 2, sum))
  }
  
  tt1<-ifelse(sup.T1-matrix(rep(T1.obs,B),B,np,byrow = T)>0,1,0)
  PvalueT1<-colSums(tt1)/B
  T1.crit.po<-(apply(sup.T1,2,sort))[round(0.95*B),]
  
  tt2<-ifelse(sup.T2-matrix(rep(T2.obs,B),B,np,byrow = T)>0,1,0)
  PvalueT2<-colSums(tt2)/B
  T2.crit.po<-(apply(sup.T2,2,sort))[round(0.95*B),]
  
  tt3<-ifelse(T3.w-matrix(rep(T3.obs,B), B, np, byrow = T)>0,1,0)
  PvalueT3<-colSums(tt3)/B
  T3.crit.po<-(apply(T3.w,2,sort))[round(0.95*B),]   
  
  out.Vmm<-cbind(uni.time.e1, V.mm)
  out.Vmm<-out.Vmm[order(out.Vmm[,1]), ]
  
  out.estb<-cbind(uni.time.e1, beta.m)
  out.estb<-out.estb[order(out.estb[,1]), ]
  
  out.sd<-cbind(uni.time.e1,  sd.m)
  out.sd<-out.sd[order(out.sd[,1]), ]
  
  colnames(out.estb)<-c("time","(Intercept)",colnames(X)[-1])
  colnames(out.sd)<-c("time","(Intercept)",colnames(X)[-1])
  
  return(list(estb=out.estb, sd=out.sd, 
              est.V.mm=out.Vmm, est.WI.mm=WI.mm[,,ind.time],
              pvalue.testeq0=PvalueT1, crit.testeq0=T1.crit.po, obs.testeq0=T1.obs,
              pvalue.KS=PvalueT2, crit.KS=T2.crit.po, obs.KS=T2.obs,
              pvalue.CvM=PvalueT3, crit.CvM=T3.crit.po, obs.CvM=T3.obs))
}

estb.Vn.twsc<-function(beta.m, X, time, event, Qt.mat, Gt, sub.NtG.m, Indx)
{
  ## for Bn ##
  n<-dim(X)[1] # = length(time) = n of subj in this regimen = dim(X)[1]
  np<-dim(X)[2]
  obs.mat<-matrix(time, n, n,byrow=T)
  delta.mat<-matrix(event, n, n,byrow=T)
  
  uni.time.e1<-time[Indx] # use for the estimation
  nt<-length(uni.time.e1)
  
  Ntt.mat<-t(ifelse(obs.mat>=time, 1, 0))  #(ixj) n x n  I(Vi>=Vj) not order
  wpi<-colMeans(Ntt.mat*Qt.mat)
  
  wdNc.mat<-t(ifelse(obs.mat==time & delta.mat == 0, 1, 0))*Qt.mat  #(ixj) n x n I(Vi=Vj, casei=0)
  wqi<-colMeans(wdNc.mat) # n-vector
  
  demo.s<-(wpi-wqi)
  demo.s[(demo.s==0)]<-1e-20

  wYc.mat<-Ntt.mat*Qt.mat
  whc<-Gt$hx
  wdAc<-t(t(wYc.mat)*whc)
  dMc<-wdNc.mat-wdAc # (ixj) nxn
  
  Gm<-t(t(dMc)/demo.s)
  Ik.mat<-Ntt.mat #t(ifelse(obs.mat>=pts, 1, 0))  #(kxi) n x n  I(Vk>=Vi) 
  
  Ikk<-matrix(rep(t(Ik.mat),n), , n, byrow=T) 
  dim(Ikk)<-c(n,n,n) 
  
  Gmm<-array(rep(Gm,each=n),c(n,n,n))
  Gm.mat<-Gmm*Ikk
  Mc.v<-apply(Gm.mat,1,sum)  # for each Ti
  
  sd.m<-c()
  V.mm<-c()
  dtau.W.mm<-array(matrix(0, n, nt), c(n, np, nt))
  inv.In.mm<-array(matrix(0, np, np),c(np,  np, nt))
  WI.mm<-array(matrix(0, n, nt),c(n, np, nt))
  
  ind.time<-order(uni.time.e1)
  rnk.time<-rank(uni.time.e1)
  ## sort beta by time
  Beta<-cbind(uni.time.e1,beta.m)[ind.time, ]
  dtau<-c(diff(c(sort(uni.time.e1),max(time))))
  rang<-diff(range(c(uni.time.e1,max(time))))
  avg.beta<-colSums(Beta[,-1]*dtau)/rang
  avg.beta.m<-matrix(avg.beta, nt, np, byrow = T)
  T2.obs<-apply(abs(beta.m-avg.beta.m),2,max)
  
  T3.obs<-colSums((beta.m-avg.beta.m)^2*dtau[rnk.time])
  ## change back
  r.dtau<-dtau[rnk.time]/rang
  
  sub.Qt<-Qt.mat[,Indx]
  
  for (indx in 1:nt)
  {
    beta<-beta.m[indx, ]
    np<-length(beta)
    e.xb<-exp(beta%*%t(X))
    P<-1-exp(-e.xb)
    D.eta<-X*c((1-P)*e.xb)
    D.eta[(e.xb==Inf),]<-0
    
    In.eta<-t(D.eta)%*%(D.eta*sub.Qt[,indx])
    
    An<-D.eta*c(sub.NtG.m[,indx]-P)*sub.Qt[,indx]
    Bn<--D.eta*c(sub.NtG.m[,indx]*Mc.v)*sub.Qt[,indx]
    
    W1<-An+Bn
    Wn1<-t(W1)%*%W1
    diag.1<-try(solve(In.eta))
    if(class(diag.1)=="try-error"){
      Vn.t<-matrix(0,np,np)
      dtau.W.mm[,,indx]<-matrix(0,n,np)*r.dtau[indx]
      WI.mm[,,indx]<-t(matrix(0,np,np)%*%t(W1))
      #inv.In.mm[,,indx]<-matrix(0, np, np)
      #dtau.m[,,indx]<-matrix(r.dtau[indx],n,np)
      #sd.est<-rep(10000,length(beta))
    }else{
      Vn.t<-solve(In.eta)%*%Wn1%*%solve(In.eta)
      dtau.W.mm[,,indx]<-t(solve(In.eta)%*%t(W1))*r.dtau[indx] #W1*r.dtau[indx]
      WI.mm[,,indx]<-t(solve(In.eta)%*%t(W1))
      #inv.In.mm[,,indx]<-solve(In.eta)
      #dtau.m[,,indx]<-matrix(r.dtau[indx],n,np)
      #sd.est<-sqrt(diag(solve(In.eta)%*%Wn1%*%solve(In.eta)))
    }
    #sd.est<-sqrt(Vn.t)
    sd.m<-rbind(sd.m, sqrt(diag(Vn.t)))
    V.mm<-rbind(V.mm, c(Vn.t))
  }
  
  T1.obs<-apply(abs(beta.m/sd.m),2,max)
  avg.dtau.w<-apply(dtau.W.mm,c(1,2),sum)    
  W.T2.m<-WI.mm-array(avg.dtau.w, c(n, np, nt))
  
  set.seed(15213)
  B<-1000
  sup.T1<-sup.T2<-T3.w<-c()
  for(b in 1:B)
  {
    gi.m<-matrix(rnorm(n, 0, 1), n, np)
    Gi.mm<-array(gi.m, c(n, np, nt))
    
    # FOR T1 #
    WG<-WI.mm*Gi.mm
    T1.w<-abs(t(apply(WG, c(2,3),sum))/sd.m)
    T1.w[is.na(T1.w[,1]),]<-0
    sup.T1<-rbind(sup.T1, apply(T1.w,2,max))
    
    # FOR T2 #
    WG.T2<-W.T2.m*Gi.mm
    T2.w<-abs(t(apply(WG.T2,c(2,3),sum)))
    T2.w[is.na(T2.w[,1]),]<-0
    sup.T2<-rbind(sup.T2, apply(T2.w,2,max))
    
    # FOR T3 #
    WG.T3<-t(apply(WG.T2, c(2,3),sum))^2*r.dtau*rang#*matrix(r.dtau*rang, nt, np)#aperm(array(matrix(r.dtau*rang, nt, np), c(nt, np, n)), c(3,2,1))
    WG.T3[is.na(WG.T3[,1]),]<-0
    T3.w<-rbind(T3.w, apply(WG.T3, 2, sum))
  }
  
  tt1<-ifelse(sup.T1-matrix(rep(T1.obs,B),B,np,byrow = T)>0,1,0)
  PvalueT1<-colSums(tt1)/B
  T1.crit.po<-(apply(sup.T1,2,sort))[round(0.95*B),]
  
  tt2<-ifelse(sup.T2-matrix(rep(T2.obs,B),B,np,byrow = T)>0,1,0)
  PvalueT2<-colSums(tt2)/B
  T2.crit.po<-(apply(sup.T2,2,sort))[round(0.95*B),]
  
  tt3<-ifelse(T3.w-matrix(rep(T3.obs,B), B, np, byrow = T)>0,1,0)
  PvalueT3<-colSums(tt3)/B
  T3.crit.po<-(apply(T3.w,2,sort))[round(0.95*B),]   
  
  out.Vmm<-cbind(uni.time.e1, V.mm)
  out.Vmm<-out.Vmm[order(out.Vmm[,1]), ]
  
  out.estb<-cbind(uni.time.e1, beta.m)
  out.estb<-out.estb[order(out.estb[,1]), ]
  
  out.sd<-cbind(uni.time.e1,  sd.m)
  out.sd<-out.sd[order(out.sd[,1]), ]
  
  colnames(out.estb)<-c("time","(Intercept)",colnames(X)[-1])
  colnames(out.sd)<-c("time","(Intercept)",colnames(X)[-1])
  
  return(list(estb=out.estb, sd=out.sd, 
              est.V.mm=out.Vmm, est.WI.mm=WI.mm[,,ind.time],
              pvalue.testeq0=PvalueT1, crit.testeq0=T1.crit.po, obs.testeq0=T1.obs,
              pvalue.KS=PvalueT2, crit.KS=T2.crit.po, obs.KS=T2.obs,
              pvalue.CvM=PvalueT3, crit.CvM=T3.crit.po, obs.CvM=T3.obs))
}

### CIF ###
CIF.sc.w<-function(est.Vn.out, predx, predt)
{
  n<-dim(est.Vn.out$est.WI.mm)[1]
  np<-dim(est.Vn.out$estb)[2]-1
  nt0<-length(predt)
  x.use<-matrix(predx,,np)
  nx<-dim(x.use)[1]
  nt<-length(time<-est.Vn.out$estb[,1])
  
  est.p1.m<-est.p1.sd.m<- dp.m<-ind.tp<-c()
  
  for (k in 1:nt0){
    ind.tp<-c(ind.tp,order(ifelse((predt[k]-time)<0,100000,predt[k]-time))[1])
  }
  
  for (indx in 1:nt ){
    beta.use<-est.Vn.out$estb[indx, -1]
    e.xb<-exp(beta.use%*%t(x.use))
    est.p1<-1-exp(-e.xb)
    dp<-c((1-est.p1)*e.xb)
    dp.m<-rbind(dp.m, dp)
    
    est.p1.sd<-sqrt((dp)^2*diag(x.use%*%matrix(est.Vn.out$est.V.mm[indx, -1],np,np)%*%t(x.use)))
    est.p1.m<-rbind(est.p1.m, est.p1)
    est.p1.sd.m<-rbind(est.p1.sd.m, est.p1.sd)
    
  }  
  
  set.seed(15213)
  B<-1000
  sup.CB<-c()
  for(b in 1:B)
  {
    gi.m<-matrix(rnorm(n, 0, 1), n, np)
    Gi.mm<-array(gi.m, c(n, np, nt))
    
    W.CB<-est.Vn.out$est.WI.mm*Gi.mm
    
    Bb<-abs(dp.m*t((x.use)%*%apply(W.CB, c(2,3),sum))/est.p1.sd.m)
    Bb[is.na(Bb[,1]),]<-0
    sup.CB<-rbind(sup.CB, apply(Bb,2,max))
  }
  cb.crit.po<-(apply(sup.CB,2,sort))[round(0.95*B),]
  
  out.p1.all<-cbind(time, est.p1.m)
  out.p1.se.all<-cbind(time, est.p1.sd.m)
  out.p1.t0<-cbind(predt, est.p1.m[ind.tp,])
  out.p1.se.t0<-cbind(predt, est.p1.sd.m[ind.tp,])
  colnames(out.p1.t0)<-colnames(out.p1.all)<-c("predt",paste0("p1.sbj",1:nx))
  colnames(out.p1.se.t0)<-colnames(out.p1.se.all)<-c("predt",paste0("p1.se.sbj",1:nx))
  
  return(list(p1.all=out.p1.all, p1.se.all=out.p1.se.all,
              p1.t0=out.p1.t0, p1.se.t0=out.p1.se.t0,
              Cb=cb.crit.po))
}

est.wSC<-function(subdata, predx, predt, int.time){
  
  time<-subdata$time
  event<-subdata$event
  Qseq<-subdata$Qseq
  
  nsub<-length(time)
  X<-as.matrix(subdata$X, nsub, ) 
  np<-dim(X)[2]
  
  Indx<-which(event==1)
  Indx2<-(((time[Indx]%in%int.time)*!duplicated(time[Indx]))==1)
  
  Gt<-wSurvKM.CHazNA(time, event, 0, Qseq)
  NtG.m<-Nt.G(time, event, case=1, Gt$Sx) 
  
  ## select event=1 and no duplicated time
  sub.NtG.m<-NtG.m[ ,Indx[Indx2]]
  uni.time.e1<-time[Indx[Indx2]]
  
  estb.sc.w.v<-c()
  for(inx in 1:sum(Indx2)){
    estb.sc.w.v<-rbind(estb.sc.w.v, multiroot(Ub.sc, start=rep(0, np), X=X, NtG=sub.NtG.m[,inx], Qseq=Qseq, maxiter = 1000)$root)
  }
  out.wsc<-estb.Vn.sc(estb.sc.w.v, X, time, event, Qseq, Gt, sub.NtG.m, uni.time.e1)
  p1.wsc<-CIF.sc.w(out.wsc, predx, predt)
  
  return(list(estb=out.wsc$estb, estb.se=out.wsc$sd, 
              est.p1.all=p1.wsc$p1.all, est.p1.se.all=p1.wsc$p1.se.all, 
              est.p1.t0=p1.wsc$p1.t0, est.p1.se.t0=p1.wsc$p1.se.t0, 
              pvalue.testeq0=out.wsc$pvalue.testeq0, crit.testeq0=out.wsc$crit.testeq0, obs.testeq0=out.wsc$obs.testeq0,
              pvalue.KS=out.wsc$pvalue.KS, crit.KS=out.wsc$crit.KS, obs.KS=out.wsc$obs.KS,
              pvalue.CvM=out.wsc$pvalue.CvM, crit.CvM=out.wsc$crit.CvM, obs.CvM=out.wsc$obs.CvM,
              CB.p1=p1.wsc$Cb
  ))
}


sc.CIF.t0<-function(sc, time.seq)
{
  time.s<-sc$time[!is.na(sc$P1)]
  P1<-sc$P1[!is.na(sc$P1)]
  P1.se<-sc$se.P1[!is.na(sc$P1)]
  
  npt<-length(time.s)
  f1<-f1.se<-c()
  for (tt in time.seq)
  {
    tpt<-max(ifelse(time.s<=tt,1,0)*c(1:npt))
    f1<-c(f1, ifelse(tpt>0, P1[tpt],P1[1]))
    f1.se<-c(f1.se, ifelse(tpt>0, P1.se[tpt],P1.se[1]))
  }
  return(list(p1=f1, p1.se=f1.se))
}


est.twSC<-function(Data.all, Qseq, predx, predt, int.time){
  
  time<-Data.all$time
  event<-Data.all$event
  Rtime<-Data.all$Rtime
    
  nsub<-length(time)
  X<-as.matrix(Data.all$X, nsub, ) 
  np<-dim(X)[2]
  
  Indx<-((event==1)*(Qseq!=0)*(time%in%int.time)==1)
  
  Gt<-twSurvKM.CHazNA(time, event, 0, Qseq, Rtime)
  NtG.Qt.m<-Nt.G.Qt(time, event, case=1, Gt$Sx, Qseq, Rtime) 
  
  ## select event=1 and in the intersect time points
  sub.NtG.m<-NtG.Qt.m$NtG.m[ ,Indx]
  Q.mat<-NtG.Qt.m$Qt.mat[ ,Indx]
  
  estb.sc.w.v<-c()
  for(inx in 1:sum(Indx)){
    estb.sc.w.v<-rbind(estb.sc.w.v, multiroot(Ub.sc, start=rep(0, np), X=X, NtG=sub.NtG.m[,inx], Qseq=Q.mat[,inx], maxiter = 1000)$root)
  }
  
  out.twsc<-estb.Vn.twsc(estb.sc.w.v, X, time, event, NtG.Qt.m$Qt.mat, Gt, sub.NtG.m, Indx)
  p1.twsc<-CIF.sc.w(out.twsc, predx, predt)
  
  return(list(estb=out.twsc$estb, estb.se=out.twsc$sd, 
              est.p1.all=p1.twsc$p1.all, est.p1.se.all=p1.twsc$p1.se.all, 
              est.p1.t0=p1.twsc$p1.t0, est.p1.se.t0=p1.twsc$p1.se.t0, 
              pvalue.testeq0=out.twsc$pvalue.testeq0, crit.testeq0=out.twsc$crit.testeq0, obs.testeq0=out.twsc$obs.testeq0,
              pvalue.KS=out.twsc$pvalue.KS, crit.KS=out.twsc$crit.KS, obs.KS=out.twsc$obs.KS,
              pvalue.CvM=out.twsc$pvalue.CvM, crit.CvM=out.twsc$crit.CvM, obs.CvM=out.twsc$obs.CvM,
              CB.p1=p1.twsc$Cb
  ))
}


sc.CIF.t0<-function(sc, time.seq)
{
  time.s<-sc$time[!is.na(sc$P1)]
  P1<-sc$P1[!is.na(sc$P1)]
  P1.se<-sc$se.P1[!is.na(sc$P1)]
  
  npt<-length(time.s)
  f1<-f1.se<-c()
  for (tt in time.seq)
  {
    tpt<-max(ifelse(time.s<=tt,1,0)*c(1:npt))
    f1<-c(f1, ifelse(tpt>0, P1[tpt],P1[1]))
    f1.se<-c(f1.se, ifelse(tpt>0, P1.se[tpt],P1.se[1]))
  }
  return(list(p1=f1, p1.se=f1.se))
}

