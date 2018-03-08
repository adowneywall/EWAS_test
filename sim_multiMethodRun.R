
### Libraries ####
library(cate)
library(dplyr)
library(reshape2)
source("reFactor.R")

### Read in Sims ####

sim.list<-read.csv("DATA/EWAS_Sims/FirstSimTestRun_2018-02-26/SimulationParameterList.csv")

SimRun<-"FirstSimTestRun_2018-02-26/"

# Options for selecting a subset of your simulation list for testing
sub.list<-subset(sim.list,sim.list$mean.B == 4)
#simRange<-seq(from=1,to=max(sim.list$X))
simRange<-sub.list$X


simRange<-c(1,2)
t1<-simTest(SimRun,simRange,sim.list,K.est = T,lfmm.test =T,cate.test = T,RFEM.test=T,refac.test = F)
sim.data<-simRead(SimRun = SimRun,sim.list[simRange[2],])
t1$refac.performance
  # est.confounder.num(~V1, sim.data2$X, as.matrix(sim.data2$Y), method = "bcv",bcv.plot = T,nRepeat = 30)$r
  # sva::num.sv(t(as.matrix(sim.data2$Y)),as.matrix(sim.data2$X), method = "be")
  # sva::num.sv(t(as.matrix(sim.data2$Y)),as.matrix(sim.data2$X), method = "leek")

#### Main Function: simTest ####
simTest <- function(SimRun,simRange,simList,K.est = FALSE,
                    lfmm.test=FALSE,cate.test=FALSE,refactor.test=FALSE,RFEM.test=FALSE,refac.test = FALSE,
                    which.K="K.bcv"){
  j<-0
  p.max<-max(simList[simRange,]$p)+2
  
  if(lfmm.test == T){
    lfmm.vals<-array(dim=c(p.max,4,length(simRange)))
    lfmm.performance<-matrix(ncol=3,nrow=length(simRange))
  }
  if(cate.test == T){
    cate.vals<-array(dim=c(p.max,4,length(simRange)))
    cate.performance<-matrix(ncol=3,nrow=length(simRange))
  }
  if(RFEM.test == T){
    RFEM.vals<-array(dim=c(p.max,4,length(simRange)))
    RFEM.performance<-matrix(ncol=3,nrow=length(simRange))
  }
  if(refac.test == T){
    refac.vals<-array(dim=c(p.max,4,length(simRange)))
    refac.performance<-matrix(ncol=3,nrow=length(simRange))
  }
  
  for(i in 1:length(simRange)){
    sim<-simList[simRange[i],]
    sim.data<-simRead(SimRun,sim)
    
    ### K estimates ####
    if(K.est == TRUE){
      ### Estimate the number of latent factors
      ### estimating latent factors
      K.ed.error<-is.numeric(try(K.ed<-est.confounder.num(~V1,sim.data$X, sim.data$Y, method = "ed"),silent = T))
      if(K.ed.error == FALSE){K.ed <- NA}
      K.bcv.error<-is.numeric(try(K.bcv<-est.confounder.num(~V1,sim.data$X,sim.data$Y, method = "bcv",bcv.plot = F)$r,silent = T))
      if(K.bcv.error == FALSE){K.bcv <- NA}
      
      ### Estimating confounders based on sva package
      num.sva.be     <- sva::num.sv(t(sim.data$Y),as.matrix(sim.data$X), method = "be")
      num.sva.leek   <- sva::num.sv(t(sim.data$Y),as.matrix(sim.data$X), method = "leek")
      
      if(i == 1){
        K.vals<-c(K.ed,K.bcv,num.sva.be,num.sva.leek)
      }else{K.vals<-rbind(K.vals,c(K.ed,K.bcv,num.sva.be,num.sva.leek))}
    }
  
    ### lfmm ####
    if(lfmm.test == TRUE){
      mod.lfmm <- lfmm::lfmm_ridge(sim.data$Y, sim.data$X, K = K.bcv)
      ## tests based on GLMs
      scores<-binomal.glm.test(sim.data,mod.lfmm$U)
      # power
      power<-mean(which(-log10(scores[,3]) > 4) %in% sim.data$CL)
      # Type one  error
      error<-mean(!(sim.data$CL %in% which(-log10(scores[,3]) > 4)))
      
      fill<-matrix(ncol=4,nrow=p.max-ncol(sim.data$Y))
      dat<-cbind(rep(x = simRange[i],times=ncol(sim.data$Y)),scores)
      
      lfmm.vals[,,i] <-rbind(dat,fill)
      lfmm.performance[i,]<-c(simRange[i],power,error)
    }
    
    ### cate ####
    if(cate.test == T){
      mod.cate <- cate(~ V1,sim.data$X, sim.data$Y, r = K.bcv,adj.method = "rr")
      ## tests based on GLMs
      scores<-binomal.glm.test(sim.data,mod.cate$Z)
      # power
      power<-mean(which(-log10(scores[,3]) > 4) %in% sim.data$CL)
      # Type one  error
      error<-mean(!(sim.data$CL %in% which(-log10(scores[,3]) > 4)))
      
      fill<-matrix(ncol=4,nrow=p.max-ncol(sim.data$Y))
      dat<-cbind(rep(x = simRange[i],times=ncol(sim.data$Y)),scores)
      
      cate.vals[,,i] <-rbind(dat,fill)
      cate.performance[i,]<-c(simRange[i],power,error)
      
    }
    
    ### reFACTor ####
    if(refac.test == T){
      mod.refac<-re(t(sim.data$Y),as.matrix(sim.data$X),K=5)
      h<-refactor(sim.data$Y,K.bcv,data_raw = T,covarfile = NULL, out = "demo_refactor")
      dim(sim.data$X)
      scores<-binomal.glm.test(sim.data,mod.refac$refactor_components)
      # power
      power<-mean(which(-log10(scores[,3]) > 4) %in% sim.data$CL)
      # Type one  error
      error<-mean(!(sim.data$CL %in% which(-log10(scores[,3]) > 4)))
      
      fill<-matrix(ncol=4,nrow=p.max-ncol(sim.data$Y))
      dat<-cbind(rep(x = simRange[i],times=ncol(sim.data$Y)),scores)
      
      refac.vals[,,i] <-rbind(dat,fill)
      refac.performance[i,]<-c(simRange[i],power,error)
    }
    ### RefFreeEWAS ####
    if(RFEM.test == T){
      mod.RFEM<-RefFreeEwasModel(t(sim.data$Y),as.matrix(sim.data$X),K=K.bcv)
      scores<-binomal.glm.test(sim.data,mod.RFEM$U[-1,])
      # power
      power<-mean(which(-log10(scores[,3]) > 4) %in% sim.data$CL)
      # Type one  error
      error<-mean(!(sim.data$CL %in% which(-log10(scores[,3]) > 4)))
      
      fill<-matrix(ncol=4,nrow=p.max-ncol(sim.data$Y))
      dat<-cbind(rep(x = simRange[i],times=ncol(sim.data$Y)),scores)
      
      RFEM.vals[,,i] <-rbind(dat,fill)
      RFEM.performance[i,]<-c(simRange[i],power,error)
      
    } 
    
    j<-j+1
    print(paste("Sim ",j," of ", length(simRange)," done!",sep=""))
  }
  return(list(
    K.vals=K.vals,
    lfmm.vals=lfmm.vals,
    lfmm.performance = lfmm.performance,
    cate.vals = cate.vals,
    cate.performance = cate.performance,
    RFEM.vals = RFEM.vals,
    RFEM.performance = RFEM.performance #,
    # refac.vals = refac.vals,
    # refac.performance = refac.performance
  ))
}

#### Short Functions ####
simRead <- function(SimRun,sim){
  local.folder<-paste("DATA/EWAS_Sims/",SimRun,"Sim",sim$X,sep="")
  U<-as.matrix(data.table::fread(Sys.glob(paste(local.folder,"/U*.csv",sep="")),drop = T))
  B<-as.matrix(data.table::fread(Sys.glob(paste(local.folder,"/B*.csv",sep="")),drop = T))
  V<-as.matrix(data.table::fread(Sys.glob(paste(local.folder,"/V*.csv",sep="")),drop = T))
  X<-data.table::fread(Sys.glob(paste(local.folder,"/X*.csv",sep="")),drop = T)
  Y<-as.matrix(data.table::fread(Sys.glob(paste(local.folder,"/Y*.csv",sep="")),drop = T))
  causalLoci<-data.table::fread(Sys.glob(paste(local.folder,"/causalLoci*.csv",sep="")),drop = T)$x
  covar<-as.matrix(data.table::fread(Sys.glob(paste(local.folder,"/covar*.csv",sep="")),drop = T))
  return(list(U=U,
              B=B,
              V=V,
              X=X,
              Y=Y,
              CL=causalLoci,
              CV=covar))
}
binomal.glm.test <- function(sim.data,lf){
  p <- NULL
  z <- NULL
  for (k in 1:ncol(sim.data$Y)){
    dim(sim.data$Y)
    mod.glm <- glm(sim.data$Y[,k] ~ ., 
                   data = data.frame(sim.data$X,lf),
                   binomial(link = "probit"))
    p[k] <- summary(mod.glm)$coeff[2,4]
    z[k] <- summary(mod.glm)$coeff[2,3] 
  }
  gif <- median(z^2)/0.456
  pcal <- pchisq(z^2/gif , df = 1, low = F)
  return(cbind(p,z,pcal))
}

#### Basic Analysis of SimTest01 on lfmm ####

#This will take a long time
#SimTest01<-simTest(SimRun,simRange,sim.list,K.est = T,lfmm.test =T)

SimTest01.list<-sim.list[simRange,]
SimTestPerfOut<-merge(SimTest01.list,data.frame(SimTest01$lfmm.performance),by.x=c("X"),by.y=c("X1"))

SimTest01$K.vals

library(ggplot2)
plotlfmm <- ggplot(SimTestPerfOut,aes(x=prop.variance,y=X3,colour=as.factor(p),shape=as.factor(mean.B),
                                      group = interaction(as.factor(p),as.factor(mean.B)))) 
plotlfmm + geom_line() + geom_point() + facet_grid(.~as.factor(n))

model.lfmm.full<-glm(X3~prop.variance+n+p+mean.B,data=SimTestPerfOut)
summary(model.lfmm.full)
model.lfmm.1<-glm(X3~prop.variance+n+mean.B,data=SimTestPerfOut)
summary(model.lfmm.1)


#### SimTest02 - multimethod

sub.list<-subset(sim.list,sim.list$mean.B == 4)
simRange<-sub.list$X

SimTest02<-simTest(SimRun,simRange,sim.list,K.est = T,lfmm.test =T,cate.test = T,RFEM.test=T)

SimTest02.list<-sim.list[simRange,]
SimTest02PerfOut.lfmm<-merge(SimTest02.list,data.frame(SimTest02$lfmm.performance),by.x=c("X"),by.y=c("X1"))

# error
plotlfmm <- ggplot(SimTest02PerfOut.lfmm,aes(x=prop.variance,y=X2,colour=as.factor(p))) 
plotlfmm + geom_line() + geom_point() + facet_grid(.~as.factor(n))
# power
plotlfmm <- ggplot(SimTest02PerfOut.lfmm,aes(x=prop.variance,y=X3,colour=as.factor(p))) 
plotlfmm + geom_line() + geom_point() + facet_grid(.~as.factor(n))

## All models
lfmm.o<-cbind("lfmm",data.frame(SimTest02$lfmm.performance)) 
cate.o<-cbind("cate",data.frame(SimTest02$cate.performance))
RFEM.o<-cbind("RFEM",data.frame(SimTest02$RFEM.performance))
coln<-c("test","sim","error","power")
colnames(lfmm.o)<-coln
colnames(cate.o)<-coln
colnames(RFEM.o)<-coln

SimT02.perf<-rbind(lfmm.o,cate.o,RFEM.o)
SimTest02PerfOut.mutl<-merge(SimTest02.list,SimT02.perf,by.x=c("X"),by.y=c("sim"))
SimT.multi<-subset(SimTest02PerfOut.mutl,SimTest02PerfOut.mutl$p == 5000)

plotlfmm <- ggplot(SimT.multi,aes(x=prop.variance,y=error,colour=as.factor(test))) 
plotlfmm + geom_line() + geom_point() + facet_grid(.~as.factor(n))



plotlfmm <- ggplot(SimT.multi,aes(x=prop.variance,y=X2,colour=as.factor(p),shape=as.factor(mean.B),
                                             group = interaction(as.factor(p),as.factor(mean.B)))) 
plotlfmm + geom_line() + geom_point() + facet_grid(.~as.factor(n))







