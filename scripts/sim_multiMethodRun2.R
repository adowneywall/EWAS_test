
### Libraries ####
library(cate)
library(dplyr)
library(reshape2)
library(dSVA)
library(fdrtool)
source("reFactor.R")

### Read in Sims ####
sim.folder<-"DATA/EWAS_Sims/Sim_Formal_2018-04-05/"

sMeta<-simMeta(sim.folder) # enter the name of folder w/ directory, returns parameters list [param] and sim folder name

# Options for selecting a subset of your simulation list for testing
sub.list<-subset(sMeta$param,sMeta$param$Sim == c(1:2))
#simRange<-seq(from=1,to=max(sim.list$X))
simRange<-sub.list$Sim

simList=sMeta$param
SimRun=sMeta$sim
t1<-simTest(sMeta$sim,simRange,sMeta$param,K.est = T,lfmm.test = T,cate.test = T,RFEM.test=F,
            refac.test = F,rep=F,rep.times = 2)

ss<-t1$full.df
me<-t1$m.e
ms<-t1$eval.tab

mean(t1$full.df$K.bcv)
mean(t1$full.df$K.be)
#### Main Function: simTest ####
simList=sMeta$param
simTest <- function(SimRun,simRange,simList,K.est = FALSE,
                    lfmm.test=FALSE,cate.test=FALSE,refactor.test=FALSE,RFEM.test=FALSE,refac.test = FALSE,
                    which.K="K.bcv",rep=T,rep.times=1){
  if(rep == F){
    if(rep.times>simList$Rep[1]){
      print("User input for # of replicates to be analyzed larger than # of replicates within simulation.")
      print(paste("Please select 'rep.times' value <= ",simList$Rep[1],sep=""))
      exit()
    } else{reps<-rep.times}
  } else{reps<-simList$Rep[1]}

  val.dir<-expand.grid(simList[simRange,]$Sim,1:reps) 
  val.dir<-val.dir[order(val.dir[,1]),]
  colnames(val.dir)<-c("Sim","Rep")
  
  for(i in 1:length(simRange)){
    print(paste("Sim ",i," of ",length(simRange),sep=""))
    
    #sim<-simList[simRange[1],]
    #reps=5
    
    if(exists("eval.tab")){rm(eval.tab)} #summary statistics
    sim<-simList[simRange[i],]
    
    print(".reading in simulation...")
    sim.data<-simRead(SimRun,sim,reps)
    
    print(".evaluating simulation...")
    temp.lfmm<-list(NULL)
    temp.cate<-list(NULL)
    
    for(k in 1:length(sim.data$X.list)){
      temp.lfmm<-list(NULL)
      print(paste("..rep ",k," of ",length(sim.data$X.list),sep=""))
      if(exists("m.e")){rm(m.e)}
      if(exists("k.vals")){rm(k.vals)}
      ### K estimates ####
      if(K.est == TRUE){
        print("...estimating K values...")
        ### Estimate the number of latent factors
        ### estimating latent factors
        K.ed.error<-is.numeric(try(K.ed<-est.confounder.num(~V1,sim.data$X.list[[k]], sim.data$Y.list[[k]], method = "ed"),silent = T))
        if(K.ed.error == FALSE){K.ed <- NA}
        K.bcv.error<-is.numeric(try(K.bcv<-est.confounder.num(~V1,sim.data$X.list[[k]],sim.data$Y.list[[k]], method = "bcv",bcv.plot = F)$r,silent = T))
        if(K.bcv.error == FALSE){K.bcv <- NA}
        ### Estimating confounders based on sva package
        K.be.error<-is.numeric(try(num.sva.be  <- sva::num.sv(t(sim.data$Y.list[[k]]),as.matrix(sim.data$X.list[[k]]), method = "be")))
        if(K.be.error == FALSE){num.sva.be <- NA}
        K.le.error<-is.numeric(try(num.sva.leek<- sva::num.sv(t(sim.data$Y.list[[k]]),as.matrix(sim.data$X.list[[k]]), method = "leek")))
        if(K.le.error == FALSE){num.sva.leek <- NA}
        
        k.vals<-cbind.data.frame(K.ed=K.ed,K.bcv=K.bcv,K.be=num.sva.be,K.leek=num.sva.leek)
      }
      
      #### lfmm ####
      if(lfmm.test == TRUE){
        print("...fitting with lfmm...")
        #mod.lfmm <- lfmm::lfmm_ridge(sim.data$Y.list[[k]], sim.data$X.list[[k]], K = 4)
        ## tests based on GLMs
        mod.lfmm <- lfmm::lfmm_ridge(sim.data$Y.list[[k]], sim.data$X.list[[k]], K = 5) ### FIX this process to allow for dynamic Ks
        scores<-binomal.glm.test(sim.data$Y.list[[k]],sim.data$X.list[[k]],mod.lfmm$U)
        fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = F)
        # power
        error<-1-mean(which(-log10(fdr.adj$qval) > 3.2) %in% sim.data$CL.list[[k]])
        # Type one  error
        power<-1-mean(!(sim.data$CL.list[[k]] %in% which(-log10(fdr.adj$qval) > 3.2)))
        
        temp.lfmm[[1]][[k]]<-list(cbind(scores,fdr.adj$pval,fdr.adj$qval))
      
        if(k == length(sim.data$X.list)){
          if(i==1){lfmm.vals<-temp.lfmm}
          else{lfmm.vals<-append(lfmm.vals,temp.lfmm)}
        }
        
        if(exists("m.e")){
        m.e<-rbind(m.e,cbind(method="lfmm",rep=k,et0=fdr.adj$param[3],power=power,error=error))}
        else{
          m.e<-cbind.data.frame(method="lfmm",rep=k,et0=fdr.adj$param[3],power=power,error=error)}
      }

      ### cate ####
      if(cate.test == T){
        print("...fitting with cate...")
        mod.cate <- cate(~ V1,sim.data$X.list[[k]], sim.data$Y.list[[k]], r = K.bcv,adj.method = "rr")
        ## tests based on GLMs
        scores<-binomal.glm.test(sim.data$Y.list[[k]],sim.data$X.list[[k]],mod.cate$Z)
        fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = F)
        # power
        error<-1-mean(which(-log10(fdr.adj$qval) > 3.2) %in% sim.data$CL.list[[k]])
        # Type one  error
        power<-1-mean(!(sim.data$CL.list[[k]] %in% which(-log10(fdr.adj$qval) > 3.2)))

        temp.cate[[1]][[k]]<-list(cbind(scores,fdr.adj$pval,fdr.adj$qval))
        
        if(k == length(sim.data$X.list)){
          if(i==1){cate.vals<-temp.cate}
          else{cate.vals<-append(cate.vals,temp.cate)}
        }
        if(exists("m.e")){
          m.e<-rbind(m.e,cbind(method="cate",rep=k,et0=fdr.adj$param[3],power=power,error=error))}else{
            m.e<-cbind.data.frame(method="cate",rep=k,et0=fdr.adj$param[3],power=power,error=error)}
      }
      # 
      # ### reFACTor ####
      # if(refac.test == T){
      #   mod.refac<-re(t(sim.data$Y),as.matrix(sim.data$X),K=5)
      #   h<-refactor(sim.data$Y,K.bcv,data_raw = T,covarfile = NULL, out = "demo_refactor")
      #   dim(sim.data$X)
      #   scores<-binomal.glm.test(sim.data,mod.refac$refactor_components)
      #   # power
      #   power<-mean(which(-log10(scores[,3]) > 4) %in% sim.data$CL)
      #   # Type one  error
      #   error<-mean(!(sim.data$CL %in% which(-log10(scores[,3]) > 4)))
      #   
      #   fill<-matrix(ncol=4,nrow=p.max-ncol(sim.data$Y))
      #   dat<-cbind(rep(x = simRange[i],times=ncol(sim.data$Y)),scores)
      #   
      #   refac.vals[,,i] <-rbind(dat,fill)
      #   refac.performance[i,]<-c(simRange[i],power,error)
      # }
      # ### RefFreeEWAS ####
      # if(RFEM.test == T){
      #   mod.RFEM<-RefFreeEwasModel(t(sim.data$Y),as.matrix(sim.data$X),K=K.bcv)
      #   scores<-binomal.glm.test(sim.data,mod.RFEM$U[-1,])
      #   # power
      #   power<-mean(which(-log10(scores[,3]) > 4) %in% sim.data$CL)
      #   # Type one  error
      #   error<-mean(!(sim.data$CL %in% which(-log10(scores[,3]) > 4)))
      #   
      #   fill<-matrix(ncol=4,nrow=p.max-ncol(sim.data$Y))
      #   dat<-cbind(rep(x = simRange[i],times=ncol(sim.data$Y)),scores)
      #   
      #   RFEM.vals[,,i] <-rbind(dat,fill)
      #   RFEM.performance[i,]<-c(simRange[i],power,error)
      #   
      # } 
      if(exists("eval.tab")){
        eval.tab<-rbind(eval.tab,cbind.data.frame(k.vals,m.e))}
      else{eval.tab<-cbind.data.frame(k.vals,m.e)}
    }
    if(exists("full.df")){
      full.df<-rbind(full.df,cbind.data.frame(sim,eval.tab))
    }
    else{full.df<-cbind.data.frame(sim,eval.tab)}
    
    print(paste("Sim ",i," of ",length(simRange)," complete",sep=""))
    print(paste("Project is ",(i/length(simRange))*100,"% complete."))
  }
  return(list(
    full.df=full.df,
    m.e=m.e,
    eval.tab=eval.tab,
    val.dir=val.dir,
    lfmm.vals=lfmm.vals,
    cate.vals = cate.vals
    # RFEM.vals = RFEM.vals,
    # RFEM.performance = RFEM.performance #,
    # refac.vals = refac.vals,
    # refac.performance = refac.performance
  ))
}

#### Short Functions ####
simMeta<-function(x){
  len=length(strsplit(sim.folder,split = "/")[[1]])
  sim=paste(strsplit(sim.folder,split = "/")[[1]][len],"/",sep="")
  return(
    list(
      param=read.csv(paste(sim.folder,"SimulationParameterList.csv",sep="")),
      sim=sim
    )
  )
}

simRead <- function(SimRun,sim,reps){
  U.list<-list()
  B.list<-list()
  V.list<-list()
  X.list<-list()
  Y.list<-list()
  causalLoci.list<-list()
  covar.list<-list()
  
    for(i in 1:reps){
      local.folder<-paste("DATA/EWAS_Sims/",SimRun,"Sim",sim$Sim,"/Rep",i,sep="")
      U.list[[i]]<-as.matrix(read.csv(Sys.glob(paste(local.folder,"/U*.csv",sep=""))))
      #read.csv(Sys.glob(paste("DATA/EWAS_Sims/Sim_TestRun_2018-04-04/Sim1/Rep1","/U*.csv",sep="")))
      U.list[[i]]<-as.matrix(read.csv(Sys.glob(paste(local.folder,"/U*.csv",sep=""))))
      B.list[[i]]<-as.matrix(read.csv(Sys.glob(paste(local.folder,"/B*.csv",sep=""))))
      V.list[[i]]<-as.matrix(read.csv(Sys.glob(paste(local.folder,"/V*.csv",sep=""))))
      X.list[[i]]<-read.csv(Sys.glob(paste(local.folder,"/X*.csv",sep="")))
      Y.list[[i]]<-as.matrix(read.csv(Sys.glob(paste(local.folder,"/Y*.csv",sep=""))))
      causalLoci.list[[i]]<-read.csv(Sys.glob(paste(local.folder,"/causalLoci*.csv",sep="")))$x
      covar.list[[i]]<-as.matrix(read.csv(Sys.glob(paste(local.folder,"/covar*.csv",sep=""))))
      #B.list[[i]]<-as.matrix(data.table::fread(Sys.glob(paste(local.folder,"/B*.csv",sep="")),drop = T))
      # V.list[[i]]<-as.matrix(data.table::fread(Sys.glob(paste(local.folder,"/V*.csv",sep="")),drop = T))
      # X.list[[i]]<-data.table::fread(Sys.glob(paste(local.folder,"/X*.csv",sep="")),drop = T)
      # Y.list[[i]]<-as.matrix(data.table::fread(Sys.glob(paste(local.folder,"/Y*.csv",sep="")),drop = T))
      # causalLoci.list[[i]]<-data.table::fread(Sys.glob(paste(local.folder,"/causalLoci*.csv",sep="")),drop = T)$x
      # covar.list[[i]]<-as.matrix(data.table::fread(Sys.glob(paste(local.folder,"/covar*.csv",sep="")),drop = T))
    }
    return(list(U.list=U.list,
                B.list=B.list,
                V.list=V.list,
                X.list=X.list,
                Y.list=Y.list,
                CL.list=causalLoci.list,
                CV.list=covar.list))
}

binomal.glm.test <- function(data.y,data.x,lf){
  p <- NULL
  z <- NULL
  for (k in 1:ncol(data.y)){
    dim(sim.data$Y)
    mod.glm <- glm(data.y[,k] ~ ., 
                   data = data.frame(data.x,lf),
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







