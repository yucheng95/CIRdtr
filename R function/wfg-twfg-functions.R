
## weighted Survival Function / Cumulative Hazard Function
wSurvKM.CHazNA<-function(Y, delta, case, Qseq) # case =0 : censored; case = 1, case = 2,...
{
  pts<-sort(unique(Y))
  n<-length(Y)
  npts<-length(pts)
  Qseq.mat<-matrix(Qseq,npts,n,byrow=T)
  obs.mat<-matrix(Y,npts,n,byrow=T)
  delta.mat<-matrix(delta,npts,n,byrow=T)
  wDt.mat<-ifelse(obs.mat==pts & delta.mat == case , Qseq.mat, 0)  # npts x n
  Yt.mat<-ifelse(obs.mat >= pts,1,0)  # npts x n
  wYt.mat<-t(t(ifelse(obs.mat >= pts,1,0))*Qseq)  # npts x n
  wlambda<-rowSums(wDt.mat)/rowSums(wYt.mat) # npts-vector
  Sxx<-cumprod((1-wlambda))  # Survival for sort.time
  Sx<-pmax(Sxx[colSums(Yt.mat)],1e-20)
  wLambda<-cumsum(wlambda)   # Cumu. Harazd for sort.time
  wHx<-wLambda[colSums(Yt.mat)] 
  whx<-wlambda[colSums(Yt.mat)]
  
  return(list(Y=Y, delta=delta, Sx=Sx, Hx=wHx, hx=whx))
}

## time-varying weighted Survival Function / Cumulative Hazard Function
twSurvKM.CHazNA<-function(Y, delta, case, Qseq, Rtime) # case =0 : censored; case = 1, case = 2,...
{
  pts<-sort(unique(Y))
  n<-length(Y)
  npts<-length(pts)
  Qseq.mat<-matrix(Qseq, npts,n,byrow=T)
  obs.mat<-matrix(Y,npts,n,byrow=T)
  Rt.mat<-matrix(Rtime, npts, n,byrow=T)
  delta.mat<-matrix(delta, npts,n,byrow=T)
  Yt.mat<-ifelse(obs.mat >= pts,1,0)  # npts x n

  # time-vary weight#
  Qt.mat<-ifelse(Rt.mat > pts, 1, Qseq.mat)  # npts x n
  twDt.mat<-ifelse(obs.mat==pts & delta.mat == case, Qt.mat, 0)  # npts x n
  twYt.mat<-ifelse(obs.mat >= pts, Qt.mat, 0)  # npts x n
  twlambda<-rowSums(twDt.mat)/rowSums(twYt.mat) # npts-vector
  twlambda<-ifelse(twlambda=="NaN",0,twlambda)
  twSxx<-cumprod((1-twlambda))  # Survival for sort.time
  twSx<-pmax(twSxx[colSums(Yt.mat)],1e-20)
  twLambda<-cumsum(twlambda)  # Cumu. Harazd  for sort.time
  twHx<-twLambda[colSums(Yt.mat)]
  twhx<-twlambda[colSums(Yt.mat)]
  
  return(list(Y=Y, delta=delta, Sx=twSx, Hx=twHx, hx=twhx))
}


# Weighted Counting process
wNt.wYt<-function(Y, delta, case, Qseq, Gx) # case = 1, case = 2,...
{
  n<-length(Y)
  obs.mat<-matrix(Y,n,n,byrow=T)
  delta.mat<-matrix(delta,n,n,byrow=T)
  Nt.mat<-t(ifelse(obs.mat<=Y & delta.mat == case , 1, 0))  #(ixj) n x n  I(Vi<=Vj, casei=1)
  dNt.mat<-t(ifelse(obs.mat==Y & delta.mat == case , 1, 0))  #(ixj) n x n I(Vi=Vj, casei=1)
  Yt.mat<-t(1-ifelse(obs.mat < Y & delta.mat == case, 1, 0))  # (ixj) nxn 1-Ni(Vj-)
  wt.mat<-t(ifelse(obs.mat >= Y , 1, ifelse(obs.mat < Y & (delta.mat !=0), outer(Gx,1/Gx), 0))) 
  return(list(wNt=Nt.mat*Qseq, dwNt=dNt.mat*Qseq, wYt=Yt.mat*Qseq, wt=wt.mat))
}

# time-varying Weighted Counting process
twNt.twYt<-function(Y, delta, case, Qseq, twGx, Rtime) # case = 1, case = 2,...
{
  pts<-Y
  n<-length(Y)
  npts<-length(pts)
  obs.mat<-matrix(Y, npts, n, byrow=T)
  delta.mat<-matrix(delta, npts, n, byrow=T)
  Nt.mat<-t(ifelse(obs.mat<=pts & delta.mat == case , 1, 0))  #(ixj) n x n  I(Vi<=Vj, casei=1)
  dNt.mat<-t(ifelse(obs.mat==pts & delta.mat == case , 1, 0))  #(ixj) n x n I(Vi=Vj, casei=1)
  Yt.mat<-t(1-ifelse(obs.mat < pts & delta.mat == case, 1, 0))  # (ixj) nxn 1-Ni(Vj-)
  
  # Time-varying weight #
  Qseq.mat<-matrix(Qseq, npts,n,byrow=T)
  Rt.mat<-matrix(Rtime, npts, n,byrow=T)
  Qt.mat<-t(ifelse(Rt.mat > pts, 1, Qseq.mat))  # Qi(Vj)
  twt.mat<-t(ifelse(obs.mat >= pts , 1, ifelse(obs.mat < pts & (delta.mat !=0), outer(twGx,1/twGx), 0))) 
  
  return(list(twNt=Nt.mat*Qt.mat, dtwNt=dNt.mat*Qt.mat, twYt=Yt.mat*Qt.mat, twt=twt.mat))
}


# Score function
Ub.fg<-function(beta, X, wt, Yt, dNt)
{
  np<-length(beta)
  wY.mat<-wt*Yt
  s0<-colMeans(wY.mat*c(exp(beta%*%t(X))))
  s1<-c()
  for(pp in 1:np)
  {
    s1<-cbind(s1, colMeans(wY.mat*c(X[,pp]*exp(beta%*%t(X)))) )
  }
  Un.fg<-colSums((X-s1/s0)*diag(wt*dNt))
  #Un.fg<-c(sum((c(X[,1]-(s1/s0)[,1])*wt*dNt)[,unik]) , sum((X[,2]-(s1/s0)[,2])*(wt*dNt)[,unik]) )
  return(Un.fg)
}

Omega.fg<-function(beta, X, wt, Yt, dNt)
{
  np<-length(beta)
  wY.mat<-wt*Yt
  s0<-colMeans(wY.mat*c(exp(beta%*%t(X))))
  s1<-c()
  for(pp in 1:np)
  {
    s1<-cbind(s1, colMeans(wY.mat*c(X[,pp]*exp(beta%*%t(X)))) )
  }
  xbar<-s1/s0
  
  Omega<-matrix(0,np,np)
  for(i in 1:np)
  { 
    for(j in i:np){
      omg<-mean((colMeans(wY.mat*c(X[,i]*X[,j]*exp(beta%*%t(X))))/s0-xbar[,i]*xbar[,j])*diag(wt*dNt))
      Omega[i,j]<-Omega[j,i]<-omg
    }
  }
  return(Omega)
}


Sigma.fg<-function(time, event, case=1, beta, X, wNY.mat, weight, Gt.Lambda)
{
  np<-length(beta)
  wY.mat<-wNY.mat$wt*wNY.mat$wYt
  pts<-time
  n<-length(time)
  npts<-length(pts)
  obs.mat<-matrix(time,npts,n,byrow=T)
  delta.mat<-matrix(event,npts,n,byrow=T)
  Yt.mat<-t(1-ifelse(obs.mat < pts & delta.mat == case, 1, 0))  # (ixj) nxn 1-Ni(Vj-)
  
  s0<-colMeans(wY.mat*c(exp(beta%*%t(X))))
  s1<-c()
  for(pp in 1:np){
    s1<-cbind(s1, colMeans(wY.mat*c(X[,pp]*exp(beta%*%t(X)))) )
  }
  xbar<-s1/s0
  h0t<-colMeans(wNY.mat$wt*wNY.mat$dwNt/s0)
  
  unik<-!duplicated(time) 
  X.arr<-array(rep(t(X),each=n), c(dim(X),n)) # x:rep, y:# variable, z: value of variable
  Xbar.arr<-array(xbar, c(dim(X),n)) # x:avg of value at t, y:#variable, z: rep
  dXxb<-X.arr-Xbar.arr # (jxi)  
  
  dA<-t(t(wNY.mat$wYt*c(exp(beta%*%t(X))))*h0t)
  dM<-wNY.mat$dwNt -dA  # dMij
  weta<-c()
  for(k in 1:np){
    weta<-cbind(weta, colSums((dXxb[,k,]*t(wNY.mat$wt*dM))[unik,]))
  }
  
  Ntt.mat<-t(ifelse(obs.mat>=pts, 1, 0))  #(ixj) n x n  I(Vi>=Vj) not order
  wpi<-colMeans(Ntt.mat*weight) # n-vector
  
  wdNc.mat<-t(ifelse(obs.mat==time & delta.mat == 0, 1, 0))*weight  #(ixj) n x n I(Vi=Vj, casei=0)
  #  wdNc.mat2<-diag(ifelse(event==0,weight,0)) #nxn
  wYc.mat<-Ntt.mat*weight 
  whc<-Gt.Lambda$hx
  #  whc<-wCHar.NA(time,event,0,weight)$hx
  wdAc<-t(t(wYc.mat)*whc)
  dMc<-wdNc.mat-wdAc # (ixj) nxn
  
  Ij.mat<-t(ifelse(obs.mat<pts,  1, 0))  #(jxi) n x n  I(Vj<Vi)
  Ik.mat<-Ntt.mat #t(ifelse(obs.mat>=pts, 1, 0))  #(kxi) n x n  I(Vk>=Vi) 
  
  Ikk<-matrix(rep(t(Ik.mat),n), , npts, byrow=T) 
  dim(Ikk)<-c(n,n,npts) 
  Ijj<-array(rep(Ij.mat,each=n),c(n,n,npts))
  Ikj.mat<-Ikk*Ijj
  rm(Ik.mat,Ij.mat,Ikk,Ijj)
  
  wphi<-c()
  for(k in 1:np){
    dXx.mat<-Ikj.mat*array(dXxb[,k,]*t(wNY.mat$wt*dM),c(n,n,npts))
    qx<-(-1)*apply(dXx.mat[unik,,],3,sum)/n
    rm(dXx.mat)
    wphi<-cbind(wphi, colSums((t(dMc)*qx/wpi)[unik,]))
  }
  rm(dXxb,Ikj.mat)
  Sigma<-cov(weta+wphi)*(n-1)/n
  
  return(list(Sigma=Sigma, wUn=weta+wphi, dM=dM, dMc=dMc, wpi=wpi, xbar=xbar))
}


Sigma.fg.tw<-function(time, event, case=1, beta, X, NY.mat.all, weight, Gt.Lambda.all, Rtime)
{
  x1.sub<-X[,1]
  x2.sub<-X[,2]
  wY.mat<-NY.mat.all$twt*NY.mat.all$twYt
  pts<-time
  n<-length(time)
  npts<-length(pts)
  obs.mat<-matrix(time,npts,n,byrow=T)
  delta.mat<-matrix(event,npts,n,byrow=T)
  Yt.mat<-t(1-ifelse(obs.mat < pts & delta.mat == case, 1, 0))  # (ixj) nxn 1-Ni(Vj-)
  
  
  s0<-colMeans(wY.mat*c(exp(beta%*%t(X))))
  s1<-cbind(colMeans(wY.mat*c(x1.sub*exp(beta%*%t(X)))),colMeans(wY.mat*c(x2.sub*exp(beta%*%t(X)))))
  xbar<-s1/s0
  #h0t<-c(1/s0*diag(NY.mat.all$twt*NY.mat.all$dtwNt)/npts)
  h0t<-colMeans(NY.mat.all$twt*NY.mat.all$dtwNt/s0)
  
  unik<-!duplicated(time)  
  X.arr<-array(rep(t(X),each=n),c(dim(X),n)) 
  Xbar.arr<-array(xbar,c(dim(X),n))
  dXxb<-X.arr-Xbar.arr # (jxi) 
  
  dA<-t(t(NY.mat.all$twYt*c(exp(beta%*%t(X))))*h0t)
  dM<-NY.mat.all$dtwNt -dA  # dMij
  deta1<-colSums((dXxb[,1,]*t(NY.mat.all$twt*dM))[unik,]) 
  deta2<-colSums((dXxb[,2,]*t(NY.mat.all$twt*dM))[unik,])
  #c(mean(deta1), mean(deta2))
  weta<-cbind(deta1,deta2)
  rm(deta1,deta2)
  
  Ntt.mat<-t(ifelse(obs.mat>=pts, 1, 0))  #(ixj) n x n  I(Vi>=Vj) not order
  
  ###
  Qseq.mat<-matrix(weight, npts,n,byrow=T)
  Rt.mat<-matrix(Rtime, npts, n,byrow=T)
  Qt.mat<-t(ifelse(Rt.mat > pts, 1, Qseq.mat))  # Qi(Vj)
  
  wpi<-colMeans(Ntt.mat*Qt.mat)
  
  wdNc.mat<-t(ifelse(obs.mat==time & delta.mat == 0, 1, 0))*Qt.mat  #(ixj) n x n I(Vi=Vj, casei=0)
  
  wYc.mat<-Ntt.mat*Qt.mat 
  whc<-Gt.Lambda.all$hx
  #  whc<-wCHar.NA(time,event,0,weight)$hx
  wdAc<-t(t(wYc.mat)*whc)
  dMc<-wdNc.mat-wdAc # (ixj)
  
  Ij.mat<-t(ifelse(obs.mat<pts,  1, 0))  #(jxi) n x n  I(Vj<Vi)
  Ik.mat<-Ntt.mat #t(ifelse(obs.mat>=pts, 1, 0))  #(kxi) n x n  I(Vk>=Vi) 
  
  Ikk<-matrix(rep(t(Ik.mat),n), , npts, byrow=T) 
  dim(Ikk)<-c(n,n,npts) 
  Ijj<-array(rep(Ij.mat,each=n),c(n,n,npts))
  Ikj.mat<-Ikk*Ijj
  rm(Ik.mat,Ij.mat,Ikk,Ijj)
  
  dX1.mat<-Ikj.mat*array(dXxb[,1,]*t(NY.mat.all$twt*dM),c(n,n,npts))
  dX2.mat<-Ikj.mat*array(dXxb[,2,]*t(NY.mat.all$twt*dM),c(n,n,npts))
  rm(dXxb,Ikj.mat)
  q1<-(-1)*apply(dX1.mat[unik,,],3,sum)/n #colSums(dXxb[,1,]*t(wNY.mat$wt*dM)*)
  q2<-(-1)*apply(dX2.mat[unik,,],3,sum)/n #colSums(dXxb[,2,]*t(wNY.mat$wt*dM))
  rm(dX1.mat,dX2.mat)
  wphi1<-colSums((ifelse(t(dMc)*q1==0,0,t(dMc)*q1/wpi))[unik,])
  wphi2<-colSums((ifelse(t(dMc)*q2==0,0,t(dMc)*q2/wpi))[unik,])
  wphi<-cbind(wphi1,wphi2)
  Sigma<-cov(weta+wphi)*(n-1)/n
  rm(wphi1,wphi2)
  return(list(Sigma=Sigma, wUn=weta+wphi, dM=dM, dMc=dMc, wpi=wpi, xbar=xbar))
}

### CIF ###
CIF.fg.w<-function(time, event, X, x, beta, wNY.mat, Sgm, omg, predt)
{
  wY.mat<-wNY.mat$wt*wNY.mat$wYt
  pts<-time
  n<-length(time)
  npts<-length(pts)
  unik<-!duplicated(time) 
  obs.mat<-matrix(time,npts,n,byrow=T)
  Ntt2.mat<-t(ifelse(obs.mat<=pts, 1, 0))  #(ixj) n x n  I(Vi<=Vj)
  s0.p<-exp(x%*%beta)
  s0<-colMeans(wY.mat*c(exp(beta%*%t(X))))
  h0t<-colMeans(wNY.mat$wt*wNY.mat$dwNt/s0)
  Ht<-colSums((Ntt2.mat*h0t)[unik,])*s0.p
  Ft<-1-exp(-Ht)
  
  dM<-Sgm$dM 
  wpi<-Sgm$wpi 
  dMc<-Sgm$dMc 
  
  j11<-t(t(wNY.mat$wt*dM)/s0)*c(s0.p) # j is time points (rep)
  
  tpp<-cpt<-c()
  for (k in 1:length(predt) ){
    tpp<-c(tpp, order(ifelse((predt[k]-time)<0,100000,predt[k]-time))[1])
  }  
  
  J11<-array(j11,c(n, n, n)) # z : rep j11
  J12<-matrix(rep(Ntt2.mat, each=n),,npts,byrow=F) # x: rep Ntt2.mat
  dim(J12)<-c(n,n,npts) 
  J1.1<-apply((J11*J12)[,unik,tpp],c(1,3),sum)    
  
  rm(j11, J11, J12)
  
  np<-length(x)
  hx0.1<-(matrix(x,n,np,byrow=T)-Sgm$xbar)*h0t
  hx0<-c()
  for(k in 1:np){
    hx0<- cbind(hx0, colSums((hx0.1[,k]*Ntt2.mat)[unik,tpp])*c(s0.p))
  }
  J2.1<-hx0%*%solve(omg)
  
  J3.1<-c()
  for (tt in predt)
  {
    Ij.mat<-t(ifelse(obs.mat<pts,  1, 0))  # (ixj) n x n  I(Vi<Vj)
    tIk.mat<-t(ifelse(obs.mat>=pts, ifelse(obs.mat<=tt,1,0), 0))  #(kxi) n x n  I(tt>=Vk>=Vi) 
    
    tIkk<-matrix(rep(t(tIk.mat),n), , npts, byrow=T) 
    dim(tIkk)<-c(n,n,npts)  # y: rep time [x x z]=tIk.mat
    Ijj<-array(rep(Ij.mat,each=n),c(n,n,npts)) # x: rep time [y x z]=Ij.mat
    tIkj.mat<-tIkk*Ijj
    
    rm(Ij.mat,tIk.mat,tIkk,Ijj)
    
    v.mat<-tIkj.mat*array(t(wNY.mat$wt*dM)/s0,c(n,n,npts)) # (u, vi, s) - repxsubjxrep
    v1<-(-1)*apply(v.mat[unik,,],3,sum)*s0.p/n 
    J3.1<-cbind(J3.1, colSums((t(dMc)*v1/wpi)[unik,]))
    rm(tIkj.mat, v.mat, v1)
  }  
  
  B<-400
  J.mat<-c()
  for ( b in 1:B)
  {
    Ai<-rnorm(n,0,1)
    J1<-colSums(J1.1*Ai)/sqrt(n)    
    J2<-J2.1%*%colSums(Sgm$wUn*Ai)/sqrt(n)
    J3<-colSums(J3.1*Ai)/sqrt(n)
    J.mat<-cbind(J.mat,(J1+J2+J3)*exp(-Ht[tpp]))
    rm(J1,J2,J3)
  }
  (sd.CIF<-sqrt(rowMeans(J.mat^2)/n))
  rm(J1.1, J2.1, J3.1, J.mat)
  return(list(Ft=Ft, Ht=Ht, sd.Ft=sd.CIF))
}


CIF.fg.tw<-function(time, event, X, x, beta.tw, NY.mat.all, Sgm, omg, predt)
{
  x1.sub<-X[,1]
  x2.sub<-X[,2]
  pts<-time
  n<-length(time)
  npts<-length(pts)
  unik<-!duplicated(time) 
  obs.mat<-matrix(time,npts,n,byrow=T)
  Ntt2.mat<-t(ifelse(obs.mat<=pts, 1, 0))  #(ixj) n x n  I(Vi<=Vj)
  
  # Time-varying weight #
  s0.p<-exp(x%*%beta.tw)
  s0<-colMeans(NY.mat.all$twt*NY.mat.all$twYt*c(exp(beta.tw%*%t(X))))
  h0t<-colMeans(NY.mat.all$twt*NY.mat.all$dtwNt/s0)
  Ht<-colSums((Ntt2.mat*h0t)[unik,])*s0.p
  Ft<-1-exp(-Ht)
  
  dM<-Sgm$dM 
  wpi<-Sgm$wpi 
  dMc<-Sgm$dMc 
  
  tpp<-c()
  for (k in 1:length(predt) )
  {
    tpp<-c(tpp, order(ifelse((predt[k]-time)<0,10,predt[k]-time))[1])
  }
  
  j11<-t(t(NY.mat.all$twt*dM)/s0)*c(s0.p)
  J11<-array(j11,c(n, n, n))
  J12<-matrix(rep(Ntt2.mat, each=n),,npts,byrow=F)
  dim(J12)<-c(n,n,npts) 
  J1.1<-apply((J11*J12)[,unik,tpp],c(1,3),sum)    
  rm(j11,J11,J12)
  
  hx0.1<-(matrix(x,n,2,byrow=T)-Sgm$xbar)*h0t
  hx0<- cbind(colSums((hx0.1[,1]*Ntt2.mat)[unik,tpp]), colSums((hx0.1[,2]*Ntt2.mat)[unik,tpp]))*c(s0.p)
  J2.1<-hx0%*%solve(omg)
  
  
  J3.1<-c()
  for (tt in predt)
  {
    Ij.mat<-t(ifelse(obs.mat<pts,  1, 0))  #(jxi) n x n  I(Vj<Vi)
    tIk.mat<-t(ifelse(obs.mat>=pts, ifelse(obs.mat<=tt,1,0), 0))  #(kxi) n x n  I(tt>=Vk>=Vi) 
    
    tIkk<-matrix(rep(t(tIk.mat),n), , npts, byrow=T) 
    dim(tIkk)<-c(n,n,npts) 
    Ijj<-array(rep(Ij.mat,each=n),c(n,n,npts))
    tIkj.mat<-tIkk*Ijj
    
    rm(Ij.mat,tIk.mat,tIkk,Ijj)
    
    v.mat<-tIkj.mat*array(t(NY.mat.all$twt*dM)/s0,c(n,n,npts))
    v1<-(-1)*apply(v.mat[unik,,],3,sum)*s0.p/n 
    J3.1<-cbind(J3.1, colSums((ifelse(t(dMc)*v1==0, 0, t(dMc)*v1/wpi))[unik,]))
    rm(tIkj.mat, v.mat, v1)
  }  
  
  B<-400
  J.mat<-c()
  for ( b in 1:B)
  {
    Ai<-rnorm(n,0,1)
    J1<-colSums(J1.1*Ai)/sqrt(n)    
    J2<-J2.1%*%colSums(Sgm$wUn*Ai)/sqrt(n)
    J3<-colSums(J3.1*Ai)/sqrt(n)
    J.mat<-cbind(J.mat,(J1+J2+J3)*exp(-Ht[tpp]))
    rm(J1,J2,J3)
  }
  sd.CIF<-sqrt(rowMeans(J.mat^2)/n)
  rm(J1.1,J2.1,J3.1,J.mat)
  return(list(Ft=Ft, Ht=Ht, sd.Ft=sd.CIF))
  
}


est.wFG<-function(dataset, id, Q, predx=c(xs1,xs2), predt=c(t1,t2))
{
  subdata<-dataset$Data.all[id, ] 
  nsub<-length(id)
  Q.sub<-Q[id]
  time<-subdata$time
  event<-subdata$event
  X<-as.matrix(subdata$X,nsub,)[ ,-1] 
  np<-dim(X)[2]
  
  Gt.Lambda<-wSurvKM.CHazNA(time, event, case=0, Q.sub)
  wNY.mat<-wNt.wYt(time, event, case=1, Q.sub, Gt.Lambda$Sx)
  estb.fg.w<-multiroot(Ub.fg, start=rep(0.3, np), X=X, wt=wNY.mat$wt, Yt=wNY.mat$wYt, dNt=wNY.mat$dwNt, maxiter = 100)$root
  omg<-Omega.fg(estb.fg.w, X, wt=wNY.mat$wt, Yt=wNY.mat$wYt, dNt=wNY.mat$dwNt)
  Sgm<-Sigma.fg(time, event, case=1, beta=estb.fg.w, X, wNY.mat, weight=Q.sub, Gt.Lambda)
  estsd<-sqrt(diag(solve(omg)%*%Sgm$Sigma%*%solve(omg)/nsub))
  
  Ft.fg.all<-CIF.fg.w(time, event, X, x=predx, estb.fg.w, wNY.mat, Sgm, omg, predt)
  Ft.fg.w<-Ft.fg.all$Ft[which(event==1)][order(time[event==1])]
  sd.Ft<-Ft.fg.all$sd.Ft
  time.fg<-sort(time[event==1])
  
  return(list(estb=estb.fg.w, p1=Ft.fg.w, time=time.fg, omg=omg, estb.se=estsd, p1.se=sd.Ft, sgm=Sgm))
}

est.twFG<-function(dataset, Q, predx=c(xs1,xs2), predt=c(t1,t2))
{
  subdata<-dataset$Data.all # for time-varying 
  nsub<-length(Q)
  time<-subdata$time
  event<-subdata$event
  Rtime<-subdata$Rtime
  #x1.sub<-subdata$x1
  #x2.sub<-subdata$x2
  #X<-cbind(x1.sub,x2.sub)
  X<-as.matrix(subdata$X,nsub,)[ ,-1] 
  
  Gt.Lambda.all<-twSurvKM.CHazNA(time, event, case=0, Q, Rtime) 
  NY.mat.all<-twNt.twYt(time, event, case=1, Q, Gt.Lambda.all$Sx, Rtime)
  #twY.mat<-NY.mat.all$twt*NY.mat.all$twYt
  
  estb.fg.tw<-multiroot(Ub.fg, start=c(0.3,0.3), X=X, wt=NY.mat.all$twt, Yt=NY.mat.all$twYt, dNt=NY.mat.all$dtwNt, maxiter = 100)$root
  omg<-Omega.fg(estb.fg.tw, X, wt=NY.mat.all$twt, Yt=NY.mat.all$twYt, dNt=NY.mat.all$dtwNt)
  Sgm<-Sigma.fg.tw(time, event, case=1, estb.fg.tw, X, NY.mat.all, Q, Gt.Lambda.all, Rtime)
  estsd<-sqrt(diag(solve(omg)%*%Sgm$Sigma%*%solve(omg)/nsub))
  
  Ft.fg.all<-CIF.fg.tw(time, event, X, x=predx, estb.fg.tw, NY.mat.all, Sgm, omg, predt)
  Ft.fg.tw<-Ft.fg.all$Ft[which(event==1)][order(time[event==1])]
  sd.Ft.tw<-Ft.fg.all$sd.Ft
  time.fg.tw<-sort(time[event==1])
  
  return(list(estb=estb.fg.tw, p1=Ft.fg.tw, time=time.fg.tw, omg=omg, estb.se=estsd, p1.se=sd.Ft.tw))
}



CIF.t0<-function(fg, time.seq)
{
  npt<-length(fg$time)
  f1<-c()
  for (tt in time.seq)
  {
    tpt<-max(ifelse(fg$time<=tt,1,0)*c(1:npt))
    f1<-c(f1, ifelse(tpt>0, fg$p1[tpt],0))
  }
  return(f1)
}
