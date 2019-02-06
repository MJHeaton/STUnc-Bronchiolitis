####################################################
## Code to Fit ST Unc Model to Bronchiolitis Data ##
####################################################

###################
## Preliminaries ##
###################

##Load Libraries
library(LatticeKrig)
library(data.table)
library(parallel)

## Set jittering length and grid size
t.jit.len <- 7 #Not true value (masked for confidentiality)
grid.box.len <- 0.0018/2 #1/2 side length of spatial discretized grid box
s.jit.len <- 0.02 #Not true value (masked for confidentiality)

## Which model to run?
include.X.unc <- TRUE
include.ST.unc <- TRUE

## Set Number of cores to use for parallel processing
ncores.unc.update <- 1
ncores.param.update <- 1
one.year <- 24  #Number of time periods per year

## Load Census information
load("CensusData.RData") #Water locations are NA

############################################################
## Load City Grid  and keep only grid cells that are land ##
############################################################
load('Norfolk200mgrid.Rdata')
water <- which(is.na(rowSums(NorfolkCensus)))
lats <- c(36.82, 36.98)
lon <- c(-76.346, -76.133)
kp.grid <- which(apply(grid.cent,1,function(x){
  length(which(x[1]==grid.cent.land[,1] & x[2]==grid.cent.land[,2]))>0
}))
lon.locs <- unique(grid.cent[,1])
lat.locs <- unique(grid.cent[,2])

####################################
## Format PM25/SO2 estimates ##
####################################
load('NorfolkFullPredictionsSO2.Rdata')
load('NorfolkFullPredictionsPM25.Rdata')

D <- rdist(grid.cent.land,coords) #assign to closest of grid coordinates
grp <- apply(D,1,which.min)

tm <- c(which(as.Date(dates)==as.Date("2003-09-01")):248)

so2 <- so2[grp,tm]
so2.var <- so2.var[grp,tm]
pm25 <- pm25[grp,tm]
pm25.var <- pm25.var[grp,tm]
so2 <- so2[,-(1:one.year)]
so2.var <- so2.var[,-(1:one.year)]
pm25 <- pm25[,-(1:one.year)]
pm25.var <- pm25.var[,-(1:one.year)]

#####################################
## Remove the first year of cases  ##
## Because we won't have counts of ##
## number of susceptibles          ## 
#####################################

## Temporal Windows
strt.dates <- c("01/01/","01/15/","02/01/","02/15/","03/01/",
                "03/15/","04/01/","04/15/","05/01/","05/15/","06/01/","06/15/","07/01/",
                "07/15/","08/01/","08/15/","09/01/","09/15/","10/01/",
                "10/15/","11/01/","11/15/","12/01/","12/15/")
fin.dates <- c("01/14/","01/31/","02/14/","02/28/","03/14/","03/31/",
               "04/14/","04/30/","05/14/","05/31/","06/14/","06/30/","07/14/","07/31/",
               "08/14/","08/31/","09/14/","09/30/","10/14/","10/31/",
               "11/14/","11/30/","12/14/","12/31/")
win.strt <- c()
win.end <- c()
for(y in 2003:2013){
  win.strt <- c(win.strt,paste(strt.dates,y,sep=""))
  win.end <- c(win.end,paste(fin.dates,y,sep=""))
}
date.cutoff <- which(win.strt=="05/01/2013")
date.strt <- which(win.strt=="08/15/2003")
win.strt <- win.strt[-(date.cutoff:length(win.strt))]
win.end <- win.end[-(date.cutoff:length(win.end))]
win.strt <- win.strt[-(1:date.strt)]
win.end <- win.end[-(1:date.strt)]
win.strt.cases <- win.strt[-(1:one.year)]
win.end.cases <- win.end[-(1:one.year)]

########################################
## Load in Data on Cases and Controls ##
########################################
cases <- read.table("./SimulatedCasesInfo.txt", header=TRUE)
cases$JitteredBDay <- as.Date(cases$JitteredBDay)
cases$JitteredInfection <- as.Date(cases$JitteredInfection)
controls <- read.table("./SimulatedControlInfo.txt", header=TRUE)
controls$JitteredBDay <- as.Date(controls$JitteredBDay)
N <- nrow(cases) + nrow(controls)

###############################################
## Calculate Overlap Probabilities for cases ##
###############################################
get.case.overlap <- function(obj){
  
  ## Spatial Overlap
  glocs <- which(rdist(obj[,c("JitteredLon","JitteredLat")],
                       grid.cent.land) < (sqrt(2)*(grid.box.len+s.jit.len)))
  glocs <- glocs[!glocs%in%water]
  gwgts <- rep(0,length(glocs))
  for(g in glocs){
    unif.pts <- cbind(runif(1000,grid.cent.land[g,1]-grid.box.len,grid.cent.land[g,1]+grid.box.len),
                      runif(1000,grid.cent.land[g,2]-grid.box.len,grid.cent.land[g,2]+grid.box.len))
    gwgts[g==glocs] <- mean(rdist(obj[,c("JitteredLon","JitteredLat")],unif.pts)<=s.jit.len)
  }
  gcell <- glocs[which.max(gwgts)]
  # plot(grid.cent.land[,1], grid.cent.land[,2])
  # points(grid.cent.land[glocs,1], grid.cent.land[glocs,2], pch=19, col="blue")
  # points(obj[,1], obj[,2], pch=19, col="red")
  
  ## Temporal Overlap
  bday.jit.set <- seq.Date(obj[,"JitteredBDay"]-t.jit.len, obj[,"JitteredBDay"]+t.jit.len, by="day")
  inf.jit.set <- seq.Date(obj[,"JitteredInfection"]-t.jit.len, obj[,"JitteredInfection"]+t.jit.len, by="day")
  wlocs <- which(abs(obj[,"JitteredBDay"]-as.Date(win.strt, format="%m/%d/%Y")) <= 8 |
                   abs(obj[,"JitteredBDay"]-as.Date(win.end, format="%m/%d/%Y")) <= 8)
  bt.wgt <- matrix(0, nrow=length(wlocs), ncol=one.year)
  for(w in wlocs){
    dseq <- seq.Date(as.Date(win.strt[w], format="%m/%d/%Y"),
                     as.Date(win.end[w], format="%m/%d/%Y"), by="day")
    doverlap <- length(base::intersect(dseq,bday.jit.set))/length(dseq)
    for(w2 in 0:(one.year-1)){
      if((w+w2)<=length(win.strt)){
        dseq <- seq.Date(as.Date(win.strt[w+w2], format="%m/%d/%Y"),
                         as.Date(win.end[w+w2], format="%m/%d/%Y"), by="day")
        bt.wgt[w==wlocs,w2+1] <- doverlap*length(base::intersect(dseq,inf.jit.set))/length(dseq)
      }
    }
  }
  
  ## Choose birthday and infection window
  bwin <- arrayInd(which.max(bt.wgt), .dim=dim(bt.wgt))
  twin <- wlocs[bwin[1]] + bwin[2] - 1
  bwin <- wlocs[bwin[1]]
  
  ## Return all the info
  return(list(glocs=glocs, gwgts=gwgts, bt.wgt=bt.wgt, gcell=gcell, bwin=bwin, twin=twin))
}
cases.sbt <- mclapply(split(cases, 1:nrow(cases)), get.case.overlap, mc.cores=ncores.param.update)

###############################################
## Calculate Overlap Probabilities for cases ##
###############################################
get.control.overlap <- function(obj){
  
  ## Spatial Overlap
  glocs <- which(rdist(obj[,c("JitteredLon","JitteredLat")],
                       grid.cent.land) < (sqrt(2)*(grid.box.len+s.jit.len)))
  glocs <- glocs[!glocs%in%water]
  gwgts <- rep(0,length(glocs))
  for(g in glocs){
    unif.pts <- cbind(runif(1000,grid.cent.land[g,1]-grid.box.len,grid.cent.land[g,1]+grid.box.len),
                      runif(1000,grid.cent.land[g,2]-grid.box.len,grid.cent.land[g,2]+grid.box.len))
    gwgts[g==glocs] <- mean(rdist(obj[,c("JitteredLon","JitteredLat")],unif.pts)<=s.jit.len)
  }
  gcell <- glocs[which.max(gwgts)]
  
  ## Temporal Overlap
  bday.jit.set <- seq.Date(obj[,"JitteredBDay"]-t.jit.len, obj[,"JitteredBDay"]+t.jit.len, by="day")
  wlocs <- which(abs(obj[,"JitteredBDay"]-as.Date(win.strt, format="%m/%d/%Y")) <= 8 |
                   abs(obj[,"JitteredBDay"]-as.Date(win.end, format="%m/%d/%Y")) <= 8)
  wwgts <- rep(0, length(wlocs))
  for(w in wlocs){
    dseq <- seq.Date(as.Date(win.strt[w], format="%m/%d/%Y"),
                     as.Date(win.end[w], format="%m/%d/%Y"), by="day")
    wwgts[w==wlocs] <- length(base::intersect(dseq,bday.jit.set))/length(dseq)
  }
  bwin <- wlocs[which.max(wwgts)]
  
  ## Return all the info
  return(list(glocs=glocs, gwgts=gwgts, wwgts=wwgts, gcell=gcell, bwin=bwin))
  
}
controls.sb <- mclapply(split(controls, 1:nrow(controls)), get.control.overlap, mc.cores=ncores.param.update)


##################################
## Useful Numbers for Reference ##
##################################
G <- nrow(grid.cent.land)
NT <- length(win.strt)
NT.case <- NT-one.year

####################
## Initiate Delta ##
####################
delta <- length(controls.sb)+length(cases.sbt)

#####################################
## Function for Quicker Tabulating ##
#####################################
quick.tab <- function(bdays,cells){
  tab <- matrix(0,nrow=NT,ncol=G)
  DT <- data.table(s=cells,b=bdays)
  agg.cnts <- DT[,.N,by=names(DT)]
  tab[(agg.cnts$s-1)*NT+agg.cnts$b] <- agg.cnts$N
  return(tab)
}
b.case <- sapply(cases.sbt, function(x){x$bwin})
b.cntrl <- sapply(controls.sb, function(x){x$bwin})
s.case <- sapply(cases.sbt, function(x){x$gcell})
s.cntrl <- sapply(controls.sb, function(x){x$gcell})
cntrl.cnts <- quick.tab(c(b.case,b.cntrl),c(s.case,s.cntrl))
N0 <- matrix(0,nrow=NT.case,ncol=G)
for(i in 1:nrow(N0)){
  N0[i,] <- colSums(cntrl.cnts[(i+1):(i+one.year),])
}
t.case <- sapply(cases.sbt, function(x){x$twin})
case.cnts <- quick.tab(t.case,s.case)
case.cnts <- case.cnts[-(1:one.year),]
kp.spat <- which(colSums(N0)>0)

###########################
## Define Basis matrices ##
###########################

## M matrix for Spatial Piece
source("GetAdjMatrix.R")

## Moran Basis Functions
A <- AdjMat(length(lat.locs),length(lon.locs))
unitnum <- A$unitlabel
A <- A$A[kp.grid,kp.grid]
ones <- matrix(1,nrow=nrow(A),ncol=1)
P.orth <- diag(nrow(A))-(ones%*%t(ones))/sum(ones)
M <- eigen(P.orth%*%A%*%P.orth)
M$vectors <- M$vectors[,order(M$values,decreasing=TRUE)]
M$values <- sort(M$values,decreasing=TRUE)
which.pos <- which(M$values>0)
Xg <-M$vectors[,1:100]
priprec.g <- t(Xg)%*%(diag(rowSums(A))-A)%*%Xg
tau.g <- 0.01

## M matrix for Temporal Effects
Aw <- matrix(0,nrow=NT-one.year,ncol=NT-one.year)
for(r in 1:(nrow(Aw)-1)){
  Aw[r,r+1] <- 1
}
Aw <- Aw+t(Aw)
ones <- matrix(1,nrow=nrow(Aw),ncol=1)
P.orth <- diag(nrow(Aw))-(ones%*%t(ones))/sum(ones)
M <- eigen(P.orth%*%Aw%*%P.orth)
M$vectors <- M$vectors[,order(M$values,decreasing=TRUE)]
M$values <- sort(M$values,decreasing=TRUE)
Xw <-M$vectors[,1:75]
priprec.w <- t(Xw)%*%(diag(rowSums(Aw))-Aw)%*%Xw
tau.w <- 0.001

#####################
## Starting Values ##
#####################
beta.vec <- matrix(c(-5,rep(0,2+ncol(NorfolkCensus))),ncol=1)
tm.log.odds <- rowSums(case.cnts)/rowSums(N0)+.001 
tm.log.odds <- log(tm.log.odds/(1-tm.log.odds))-beta.vec[1]
eta.vec <- matrix(solve(t(Xw)%*%Xw+tau.w*priprec.w)%*%t(Xw)%*%tm.log.odds,nrow=ncol(Xw),ncol=1)
spat.log.odds <- colSums(case.cnts)/colSums(N0)+0.001
spat.log.odds <- log(spat.log.odds/(1-spat.log.odds))-beta.vec[1]
psi.vec <- matrix(solve(t(Xg[kp.spat,])%*%Xg[kp.spat,]+tau.g*priprec.g)%*%t(Xg[kp.spat,])%*%spat.log.odds[kp.spat],nrow=ncol(Xg),ncol=1)

############################################
## Functions to Evaluate Like within each ##
############################################

control.like.diff <- function(x,bvec=beta.vec,eta=eta.vec,psi=psi.vec){
  ## Find valid times
  times <- which(win.strt.cases%in%win.strt[x$bwin+(0:(one.year-1))])
  
  ## Calculate logit(prob) from proposal
  logit.prob.prop <- bvec$prop[1]+bvec$prop[2]*cur.so2[x$gcell,times]+bvec$prop[3]*cur.pm25[x$gcell,times] + 
    c(bvec$prop[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,])) +
    Xw[times,]%*%eta$prop+sum(Xg[x$gcell,]*psi$prop)
  
  ## Calculate logit(prob) from current
  logit.prob.cur <- bvec$cur[1]+bvec$cur[2]*cur.so2[x$gcell,times]+bvec$cur[3]*cur.pm25[x$gcell,times] + 
    c(bvec$cur[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,])) +
    Xw[times,]%*%eta$cur+sum(Xg[x$gcell,]*psi$cur)
  
  ## Back transform to probability
  prob.prop <- exp(logit.prob.prop)/(1+exp(logit.prob.prop))
  prob.cur <- exp(logit.prob.cur)/(1+exp(logit.prob.cur))
  
  ## Sum the 1-prob because we didn't see a case here
  return(sum(log(1-prob.prop))-sum(log(1-prob.cur)))
}

case.like.diff <- function(x,bvec=beta.vec,eta=eta.vec,psi=psi.vec){
  ## Find valid times
  times <- which(win.strt.cases%in%win.strt[x$bwin+(0:(one.year-1))])
  if(length(times)==0){
    return(0) 
  } else {
    case.time <- which(win.strt.cases%in%win.strt[x$twin]) #-24
    if(length(case.time)>0){
      times <- times[times!=case.time]
    }
    
    ## Calculate logit(prob) under proposed values
    logit.prob.prop <- bvec$prop[1]+bvec$prop[2]*cur.so2[x$gcell,times]+bvec$prop[3]*cur.pm25[x$gcell,times]+
      c(bvec$prop[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,]))+
      Xw[times,]%*%eta$prop+sum(Xg[x$gcell,]*psi$prop)
    if(length(case.time)!=0){
      c.logit.prob.prop <- bvec$prop[1]+bvec$prop[2]*cur.so2[x$gcell,case.time]+bvec$prop[3]*cur.pm25[x$gcell,case.time]+
        c(bvec$prop[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,])) +
        sum(Xw[case.time,]*eta$prop)+sum(Xg[x$gcell,]*psi$prop)
      c.prob.prop <- exp(c.logit.prob.prop)/(1+exp(c.logit.prob.prop))
    } else {
      c.prob.prop <- 1
    }
    
    ## Calculate logit(prob) under current values
    logit.prob.cur <- bvec$cur[1]+bvec$cur[2]*cur.so2[x$gcell,times]+bvec$cur[3]*cur.pm25[x$gcell,times]+
      c(bvec$cur[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,]))+
      Xw[times,]%*%eta$cur+sum(Xg[x$gcell,]*psi$cur)
    if(length(case.time)!=0){
      c.logit.prob.cur <- bvec$cur[1]+bvec$cur[2]*cur.so2[x$gcell,case.time]+bvec$cur[3]*cur.pm25[x$gcell,case.time]+
        c(bvec$cur[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,])) +
        sum(Xw[case.time,]*eta$cur)+sum(Xg[x$gcell,]*psi$cur)
      c.prob.cur <- exp(c.logit.prob.cur)/(1+exp(c.logit.prob.cur))
    } else {
      c.prob.cur <- 1
    }
    
    ## Back transform to probability
    prob.prop <- exp(logit.prob.prop)/(1+exp(logit.prob.prop))
    prob.cur <- exp(logit.prob.cur)/(1+exp(logit.prob.cur))
    
    ## Sum the log(prob) @ time + log(1-prob) because we didn't see a case on other time periods
    return(log(c.prob.prop)+sum(log(1-prob.prop))-(log(c.prob.cur)+sum(log(1-prob.cur))))
  }
}

###################
## Update Lambda ##
###################
update.lam <- function(cnts){
  gamvars <- matrix(rgamma(length(cnts),c(cnts)+1,1),nrow=nrow(cnts),ncol=ncol(cnts))
  return(gamvars/sum(gamvars))
}
system.time(Lambda <- update.lam(cntrl.cnts))

#########################################
## Draw time windows/regions for cases ##
#########################################
update.both.case <- function(x,bvec=beta.vec,eta=eta.vec,psi=psi.vec) { 
  ## Update Location
  case.time <- x$twin-24
  if (case.time <= 0) return(x) #if the case before 9/1/04
  logit.probs <- bvec[1]+bvec[2]*cur.so2[x$glocs,case.time]+bvec[3]*cur.pm25[x$glocs,case.time]+
    t(bvec[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$glocs,]))+
    (Xw[case.time,]%*%eta)[1]+(Xg[x$glocs,]%*%psi)
  if(length(Xw[case.time,]%*%eta)>1) cat("Error: too long")
  probs <- (exp(logit.probs)/(1+exp(logit.probs)))*Lambda[x$bwin,x$glocs]*x$gwgts
  x$gcell <- sample(x$glocs,size=1,prob=probs)
  
  ##Update b/t
  if(length(x$wlocs)>1) {
    z <- c(x$wlocs[1]+0:23,x$wlocs[2]+0:23)
  } else { 
    z <- x$wlocs+0:23
  }
  
  ## get valid case/control times
  z.case <- z-24
  stf <- which(z.case>0 & z.case < 209) #only sample cases after 9/1/04
  z.case <- z.case[stf]
  if (length(z.case) < 2) return(x) #if only 1 time or no possible times
  
  stf2 <- which(z < 233) #only births before 5/1/13
  z <- z[stf2]
  
  ## get thetas/lambdas
  logit.theta <- bvec[1]+bvec[2]*cur.so2[x$gcell,z.case]+bvec[3]*cur.pm25[x$gcell,z.case]+
    bvec[4:(3+ncol(NorfolkCensus))]%*%t(NorfolkCensus[x$gcell,])+
    (Xw[z.case,]%*%eta)[1]+Xg[x$gcell,]%*%psi
  theta <- exp(logit.theta)/(1+exp(logit.theta))
  lambda <- Lambda[z,x$gcell]
  
  thta <- rep(0,48) #the impossible cases (to early/late) get probability of 0
  thta[stf] <- theta
  
  lmbda <- rep(0,48)
  lmbda[stf2] <- lambda 
  
  probs <- lmbda*thta*as.vector(t(x$bt.wgt))
  draw <- sample(1:length(probs),size=1,prob=probs) 
  ind <- arrayInd(draw,c(24,length(x$wlocs)))
  x$bwin <- x$wlocs[ind[2]]
  x$twin <- x$bwin+ind[1]-1
  
  return(x)
}

############################################
## Draw time windows/regions for controls ##
############################################

update.both.cntrl <- function(x) {
  probs <- x$wwgts*Lambda[x$wlocs,x$gcell]
  if (length(x$wlocs) > 1) x$bwin <- sample(x$wlocs,size=1,prob=probs) #times lambda?
  
  probs <- x$gwgts*Lambda[x$bwin,x$glocs]
  if (length(x$glocs) > 1) x$gcell <- sample(x$glocs,size=1,prob=probs)
  return(x)
}

# system.time(tst <- lapply(controls.sb,update.both.cntrl))

###################
## MCMC Settings ##
###################
source("AMCMCUpdate.R")
burn <- 1000000
num.it <- 1000000
thin <- 1
beta.ind <- 1:length(beta.vec)
eta.ind <- length(beta.vec)+(1:length(eta.vec))
psi.ind <- max(eta.ind)+(1:length(psi.vec))
amcmc <- list(mn=matrix(0,nrow=max(psi.ind),ncol=1),
              var=matrix(0,nrow=max(psi.ind),ncol=max(psi.ind)))
amcmc.it <- 250
kp <- 0
kpseq <- seq(burn+thin,burn+thin*num.it,by=thin)
save.seq <- round(seq(.02,.98,by=.02)*(burn))
printseq <- kpseq[round(seq(1,length(kpseq),length=20))]

############################
## Matrices to Hold Draws ##
############################
beta.draws <- matrix(0,nrow=num.it,ncol=length(beta.vec))
eta.draws <- matrix(0,nrow=num.it,ncol=nrow(eta.vec))
psi.draws <- matrix(0,nrow=num.it,ncol=nrow(psi.vec))

#####################
## Begin MCMC Loop ##
#####################
system.time({
  for(it in 1:(burn+thin*num.it)){
    
    if(include.X.unc){
      ## Draw new pollution - X matrix
      cur.so2 <- matrix(so2 + sqrt(so2.var)*rnorm(length(so2)),
                        nrow=nrow(so2),ncol=ncol(so2))
      cur.pm25 <- matrix(pm25 + sqrt(pm25.var)*rnorm(length(pm25)),
                         nrow=nrow(so2),ncol=ncol(so2))
    }
    
    ## Update lambda
    s.case <- sapply(cases.sbt,function(x) x$gcell)
    b.case <- sapply(cases.sbt,function(x) x$bwin)
    s.cntrl <- sapply(controls.sb,function(x) x$gcell)
    b.cntrl <- sapply(controls.sb,function(x) x$bwin) 
    cntrl.cnts <- quick.tab(c(b.case,b.cntrl),c(s.case,s.cntrl))
    Lambda <- update.lam(cntrl.cnts)
    
    if(include.ST.unc){
      ## Update Control Spatial Locations and Bdays
      controls.sb <- mclapply(controls.sb,update.both.cntrl,mc.cores=ncores.unc.update)
    
      ## Update Cases Spatial Locations, BDays, and RSV Dates jointly
      cases.sbt <- mclapply(cases.sbt,update.both.case,mc.cores=ncores.unc.update)
    }
    
    ## Update beta, eta and psi jointly
    prop.var <- (0.00001^2)*diag(nrow(amcmc$mn))
    if(it>amcmc.it){
      prop.var <- (2.4^2/nrow(amcmc$mn))*(amcmc$var+prop.var)
    }
    all.prop <- c(beta.vec,eta.vec,psi.vec)+t(chol(prop.var))%*%rnorm(nrow(amcmc$mn))
    prop.beta <- matrix(all.prop[beta.ind],ncol=1)
    prop.eta <- matrix(all.prop[eta.ind],ncol=1)
    prop.psi <- matrix(all.prop[psi.ind],ncol=1)
    bvecs <- list(prop=prop.beta,cur=beta.vec)
    etavecs <- list(prop=prop.eta,cur=eta.vec)
    psivecs <- list(prop=prop.psi,cur=psi.vec)
    cntrl.llike.diff <- sum(as.numeric(mclapply(controls.sb,control.like.diff,bvec=bvecs,eta=etavecs,psi=psivecs,mc.cores=ncores.param.update)))
    case.llike.diff <- sum(as.numeric(mclapply(cases.sbt,case.like.diff,bvec=bvecs,eta=etavecs,psi=psivecs,mc.cores=ncores.param.update)))
    prior.prop <- sum(dnorm(prop.beta,0,5,log=TRUE))-(tau.w/2)*t(prop.eta)%*%priprec.w%*%prop.eta-
      (tau.g/2)*t(prop.psi)%*%priprec.g%*%prop.psi
    prior.cur <- sum(dnorm(beta.vec,0,5,log=TRUE))-(tau.w/2)*t(eta.vec)%*%priprec.w%*%eta.vec-
      (tau.g/2)*t(psi.vec)%*%priprec.g%*%psi.vec
    MH.ratio <- cntrl.llike.diff+case.llike.diff+prior.prop-prior.cur
    if(log(runif(1,0,1))<MH.ratio){
      beta.vec <- prop.beta
      eta.vec <- prop.eta
      psi.vec <- prop.psi
    }
    amcmc <- AMCMC.update(rbind(beta.vec,eta.vec,psi.vec),amcmc$mn,amcmc$var,it)
    
    
    ## Update Delta
    delta <- rgamma(1,N+1/2,1)
    
    ## Retain Draw if Necessary
    if(it%in%kpseq){
      kp <- kp+1
      if(it%in%printseq){
        cat(paste(round(100*kp/num.it,2),"of Draws Obtained (Post Burn)\n"))
      }
      beta.draws[kp,] <- beta.vec
      eta.draws[kp,] <- eta.vec
      psi.draws[kp,] <- psi.vec
      save(file="./RSVResults.RData",list=c("beta.draws","eta.draws","psi.draws","Xg","Xw"))
    }
    
    if(it%in%save.seq){
      ## Save current draw so can start where left off if necessary
      save(file="./CurrentDraw.RData",list=c("beta.vec","psi.vec","eta.vec"))
    }
    
    
  }
})









