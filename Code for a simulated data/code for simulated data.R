library(rootSolve)
library(cmprsk)
library(timereg)

source("wfg-twfg-functions.R")
source("wsc-twsc-functions.R")

data<-read.csv("exampledata.csv",header=T,sep=",")

##########################################################################
## The data in this example has the same format as that in the real neuroblastoma data reported in the data analysis section of our paper. More specifically,
## There are two initial treatments A1 and A2 with 5 covariates, such as x1 (continuous), and x2-x5(binary) at baseline.
## Participants who had responded to the initial treatment Aj, j = 1 or 2, were further randomized to treatment B1 or B2.
## Participants who had not responded to the initial treatment received no further treatment, denoted as Bp.
## Therefore, this study design contains 4 dynamic treatment regimens : A1B1Bp, A1B2Bp, A2B1Bp and A2B2Bp.
## Similar to the real data analysis section, we first used the fixed-weight Scheike model.
## After the multiple test adjustment, we then used the fixed-weight Fine and Gray model as our final model.
## The following code will produce the graph for the estimated CIF using the WFG model under different covariates values.
##########################################################################

A1<-data[data$initrt_A==1,]
A2<-data[data$initrt_A==2,]

R.A1<-A1$R
Z1.A1<-A1$Z1
Z2.A1<-A1$Z2

R.A2<-A2$R
Z1.A2<-A2$Z1
Z2.A2<-A2$Z2

#  for all subj
regimeind1.A1<-R.A1*Z1.A1+Z2.A1      # A1B1Bp
regimeind2.A1<-R.A1*(1-Z1.A1)+Z2.A1  # A1B2Bp   

regimeind1.A2<-R.A2*Z1.A2+Z2.A2       # A2B1Bp        
regimeind2.A2<-R.A2*(1-Z1.A2)+Z2.A2   # A2B2Bp      

id1.A1<-which(regimeind1.A1==1)
id2.A1<-which(regimeind2.A1==1)
id1.A2<-which(regimeind1.A2==1)
id2.A2<-which(regimeind2.A2==1)

# calculate hat_Q for two regimes 
prhat.A1<-mean(R.A1)
p1hat.A1<-sum(Z1.A1)/sum(R.A1) 
p1p.hat.A1<-1-p1hat.A1 

prhat.A2<-mean(R.A2)
p1hat.A2<-sum(Z1.A2)/sum(R.A2) 
p1p.hat.A2<-1-p1hat.A2 

# weights for all subj
Q1.A1<-R.A1*Z1.A1/p1hat.A1+Z2.A1
Q2.A1<-R.A1*(1-Z1.A1)/p1p.hat.A1+Z2.A1

Q1.A2<-R.A2*Z1.A2/p1hat.A2+Z2.A2
Q2.A2<-R.A2*(1-Z1.A2)/p1p.hat.A2+Z2.A2

# subj for specific path
a1b1.A1<-R.A1*Z1.A1           #  a1b1.A1 
a1b2.A1<-R.A1*(1-Z1.A1)       #  a1b2.A1

id.a1b1.A1<-which(a1b1.A1==1)
id.a1b2.A1<-which(a1b2.A1==1)  
id.a1bp.A1<-which(Z2.A1==1)  

a1b1.A2<-R.A2*Z1.A2           #  a1b1.A2 
a1b2.A2<-R.A2*(1-Z1.A2)       #  a1b2.A2 

id.a1b1.A2<-which(a1b1.A2==1)
id.a1b2.A2<-which(a1b2.A2==1)  
id.a1bp.A2<-which(Z2.A2==1)  


## Outcome and Event ##
time.A1<-A1$time
event.A1<-A1$event 
table(event.A1) # 59 had event of interest, 9 died, 32 censored

time.A2<-A2$time
event.A2<-A2$event 
table(event.A2) # 91 had event of interest, 5 died, 24 censored

## covariates ##
A1.cov<-A1[ ,6:10]
A2.cov<-A2[ ,6:10]

######################################
### multivariate analysis ############
Data.all.A1<-data.frame(time.A1, event.A1, 1, 10000)
colnames(Data.all.A1)<-c("time", "event", "S", "Rtime")
Data.all.A1$X<-cbind(1,A1.cov)
data.A1<-list( Data.all=Data.all.A1 )

Data.all.A2<-data.frame(time.A2, event.A2, 1, 10000)
colnames(Data.all.A2)<-c("time", "event", "S", "Rtime")
Data.all.A2$X<-cbind(1,A2.cov)
data.A2<-list( Data.all=Data.all.A2 )

t0<-c(1:3)*365

#### Fixed-weight Scheike model  ########################
##### For A1B1Bp 
multdata.A1.1<-data.A1$Data.all[id1.A1,]
multdata.A1.1$Qseq<-Q1.A1[id1.A1]

sc.A1.1<-comp.risk(Event(time, event)~X$x1+X$x2+X$x3+X$x4+X$x5, data=multdata.A1.1, cause=1, n.sim =100, n.times = NULL, model="prop")
intersect.time<-sc.A1.1$cum[!is.na(sc.A1.1$cum[,2]),1] # unique time to use

predx<-round(apply(data.matrix(multdata.A1.1$X),2,median))
wsc.A1.1<-est.wSC(multdata.A1.1, predx, predt=t0, int.time=intersect.time)
test11<-rbind(wsc.A1.1$pvalue.KS[-1],wsc.A1.1$pvalue.CvM[-1])
row.names(test11)<-c("KS","CvM")
test11

################################
######### WFG models ###########
################################

x0<-round(apply(data.matrix(data.A1$Data.all$X),2,median))[-1]

# Fixed weight FG function for est b and CIF for event 1 
wfg.A1.1<-est.wFG(dataset=data.A1, id=id1.A1, Q=Q1.A1, predx=x0, predt=t0)  
wfg.A1.2<-est.wFG(dataset=data.A1, id=id2.A1, Q=Q2.A1, predx=x0, predt=t0)  
wfg.A2.1<-est.wFG(dataset=data.A2, id=id1.A2, Q=Q1.A2, predx=x0, predt=t0)   
wfg.A2.2<-est.wFG(dataset=data.A2, id=id2.A2, Q=Q2.A2, predx=x0, predt=t0) 
pvalue<-round(rbind(2*(1-pnorm(abs(wfg.A1.1$estb)/wfg.A1.1$estb.se)),
                    2*(1-pnorm(abs(wfg.A1.2$estb)/wfg.A1.2$estb.se)),
                    2*(1-pnorm(abs(wfg.A2.1$estb)/wfg.A2.1$estb.se)),
                    2*(1-pnorm(abs(wfg.A2.2$estb)/wfg.A2.2$estb.se))),3)
Est.b.sd<-rbind(paste0(round(wfg.A1.1$estb,2),"(",round(wfg.A1.1$estb.se,2),")")
                ,paste0(round(wfg.A1.2$estb,2),"(",round(wfg.A1.2$estb.se,2),")")
                ,paste0(round(wfg.A2.1$estb,2),"(",round(wfg.A2.1$estb.se,2),")")
                ,paste0(round(wfg.A2.2$estb,2),"(",round(wfg.A2.2$estb.se,2),")"))

rownames(pvalue)<-rownames(Est.b.sd)<-c("A1B1Bp","A1B2Bp","A2B1Bp","A2B2Bp")
colnames(pvalue)<-colnames(Est.b.sd)<-c("x1","x2","x3","x4","x5")
Est.b.sd
pvalue
