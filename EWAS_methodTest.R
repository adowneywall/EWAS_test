
### Libraries ####
# library(cate)
# library(dplyr)
# library(reshape2)
# library(dSVA)
# library(sva)
# library(fdrtool)
# library(RefFreeEWAS)
# source("reFactor.R")
# install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock'))

#### Main Function: simTest ####

simTest<- function(SimRun,simRange,simList,K.est = FALSE,
                    lfmm.test=FALSE,cate.test=FALSE,refactor.test=FALSE,RFEM.test=FALSE,refac.test = FALSE,
                    SVA.test = FALSE, dSVA.test = FALSE,oracle.test=FALSE,glm.test=FALSE,
                    which.K="K.bcv",rep=T,rep.times=1,gif=3.5,fdr=0.05){
  if(rep == T){
    if(rep.times>simList$Rep[1]){
      print("User input for # of replicates to be analyzed larger than # of replicates within simulation.")
      print(paste("Please select 'rep.times' value <= ",simList$Rep[1],sep=""))
      exit()
    } else{reps<-rep.times}
  } else{reps<-simList$Rep[1]}
  
  val.dir<-expand.grid(simList[simRange,]$Sim,1:reps) 
  val.dir<-val.dir[order(val.dir[,1]),]
  colnames(val.dir)<-c("Sim","Rep")
  if(exists("full.df")){rm(full.df)}
  
  for(i in 1:length(simRange)){
    print(paste("Sim ",i," of ",length(simRange),sep=""))
    
    #sim<-simList[simRange[1],]
    #reps=2
    
    sim<-simList[simRange[i],]
    
    print(".reading in simulation...")
    sim.data<-NULL
    sim.data<-simRead(SimRun,sim,reps)
    print(".evaluating simulation...")
    temp.lfmm<-list(NULL)
    temp.cate<-list(NULL)
    temp.sva<-list(NULL)
    temp.dsva<-list(NULL)
    temp.RFEM<-list(NULL)
    temp.oracle<-list(NULL)
    temp.glm<-list(NULL)
    
    if(exists("eval.tab")){rm(eval.tab)} #summary statistics
    for(k in 1:length(sim.data$X.list)){
      print(paste("..rep ",k," of ",length(sim.data$X.list),sep=""))
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
      
      if(exists("m.e")){rm(m.e)}
      #### lfmm ####
      if(lfmm.test == TRUE){
        print("...fitting with lfmm...")
        ## tests based on GLMs
        
        mod.lfmm <- lfmm::lfmm_ridge(sim.data$Y.list[[k]], sim.data$X.list[[k]], K = sim$K) ### FIX this process to allow for dynamic Ks
        #mod.lfmm <- lfmm::lfmm_ridge(sim.data$Y.list[[k]], sim.data$X.list[[k]], K = num.sva.be) ### FIX this process to allow for dynamic Ks
        scores<-binomal.glm.test(sim.data$Y.list[[k]],sim.data$X.list[[k]],mod.lfmm$U)
        fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = F)
        
        pow_err<-powErr(fdr.adj,sim.data$CL.list[[k]],gif,fdr)
        
        temp.lfmm[[1]][[k]]<-list(cbind(scores,fdr.adj$pval,fdr.adj$qval))
      
        if(k == length(sim.data$X.list)){
          if(i==1){lfmm.vals<-temp.lfmm}
          else{lfmm.vals<-append(lfmm.vals,temp.lfmm)}
        }
        
        if(exists("m.e")){
        m.e<-rbind(m.e,cbind(method="lfmm",rep=k,et0=fdr.adj$param[3],
                             power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4]))}
        else{
          m.e<-cbind.data.frame(method="lfmm",rep=k,et0=fdr.adj$param[3],
                                power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4])}
      }
      ### cate ####
      if(cate.test == T){
        print("...fitting with cate...")
        #mod.cate <- cate(~ V1,sim.data$X.list[[k]], sim.data$Y.list[[k]], r = num.sva.be,adj.method = "rr")
        mod.cate <- cate(~ V1,sim.data$X.list[[k]], sim.data$Y.list[[k]], r = sim$K,adj.method = "rr")
        ## tests based on GLMs
        scores<-binomal.glm.test(sim.data$Y.list[[k]],sim.data$X.list[[k]],mod.cate$Z)
        fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = F)
        pow_err<-powErr(fdr.adj,sim.data$CL.list[[k]],gif,fdr)

        temp.cate[[1]][[k]]<-list(cbind(scores,fdr.adj$pval,fdr.adj$qval))
        
        if(k == length(sim.data$X.list)){
          if(i==1){cate.vals<-temp.cate}
          else{cate.vals<-append(cate.vals,temp.cate)}
        }
        if(exists("m.e")){
          m.e<-rbind(m.e,cbind(method="cate",rep=k,et0=fdr.adj$param[3],
                               power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4]))}
        else{m.e<-cbind.data.frame(method="cate",rep=k,et0=fdr.adj$param[3],
                                   power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4])}
      }
      
      #### SVA ####
      if(SVA.test == T){
        print("...fitting with SVA...")
        ## Full model
        # variables of interest + adjustment variables
        mod <- model.matrix(~V1, data=sim.data$X.list[[k]])#in this case we only have cancer (variable of interest)
        # NOTE: cancer is included as.factor() here since it has multiple levels

        ## Null model
        # adjustment variables
        mod0 <- model.matrix(~1,data=sim.data$X.list[[k]])

        #Estimate the surrogate variables using sva
        mod.sva <-sva(t(sim.data$Y.list[[k]]),mod,mod0,n.sv=sim$K)
        scores<-binomal.glm.test(sim.data$Y.list[[k]],sim.data$X.list[[k]],mod.sva$sv)

        fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = F)
        pow_err<-powErr(fdr.adj,sim.data$CL.list[[k]],gif,fdr)
        
        temp.sva[[1]][[k]]<-list(cbind(scores,fdr.adj$pval,fdr.adj$qval))

        if(k == length(sim.data$X.list)){
          if(i==1){sva.vals<-temp.sva}
          else{sva.vals<-append(sva.vals,temp.sva)}
        }
        if(exists("m.e")){
          m.e<-rbind(m.e,cbind(method="sva",rep=k,et0=fdr.adj$param[3],
                               power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4]))}
        else{m.e<-cbind.data.frame(method="sva",rep=k,et0=fdr.adj$param[3],
                                   power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4])}
      }
      #### dSVA ####
      if(dSVA.test == T){
        print("...fitting with dSVA...")
        mod.dSVA<-dSVA(Y = sim.data$Y.list[[k]],X = as.matrix(sim.data$X.list[[k]]),ncomp = sim$K)
        
        scores<-binomal.glm.test(sim.data$Y.list[[k]],sim.data$X.list[[k]],mod.dSVA$Z)
        
        fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = F)
        pow_err<-powErr(fdr.adj,sim.data$CL.list[[k]],gif,fdr)
        
        temp.dsva[[1]][[k]]<-list(cbind(scores,fdr.adj$pval,fdr.adj$qval))
        #print(temp.dsva)
        if(k == length(sim.data$X.list)){
          if(i==1){dsva.vals<-temp.dsva}
          else{dsva.vals<-append(dsva.vals,temp.dsva)}
        }
        if(exists("m.e")){
          m.e<-rbind(m.e,cbind(method="dsva",rep=k,et0=fdr.adj$param[3],
                               power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4]))}
        else{m.e<-cbind.data.frame(method="dsva",rep=k,et0=fdr.adj$param[3],
                                   power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4])}
      }
      #### reFACTor ####
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
      #### RefFreeEWAS ####
      if(RFEM.test == T){
        print("...fitting with RFEM...")
        mod.RFEM<-RefFreeEwasModel(t(sim.data$Y.list[[k]]),as.matrix(sim.data$X.list[[k]]),K=sim$K)
        scores<-binomal.glm.test(sim.data$Y.list[[k]],sim.data$X.list[[k]],mod.RFEM$U[-1,])

        fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = F)
        pow_err<-powErr(fdr.adj,sim.data$CL.list[[k]],gif,fdr)
        
        temp.RFEM[[1]][[k]]<-list(cbind(scores,fdr.adj$pval,fdr.adj$qval))
        
        if(k == length(sim.data$X.list)){
          if(i==1){RFEM.vals<-temp.RFEM}
          else{RFEM.vals<-append(RFEM.vals,temp.RFEM)}
        }
        if(exists("m.e")){
          m.e<-rbind(m.e,cbind(method="RFEM",rep=k,et0=fdr.adj$param[3],
                               power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4]))}
        else{m.e<-cbind.data.frame(method="RFEM",rep=k,et0=fdr.adj$param[3],
                                   power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4])}

      }
      #### Oracle ####
      if(oracle.test == T){
        print("...fitting with oracle...")
        mod.oracle<-sim.data$U.list[[k]]
        
        scores<-binomal.glm.test(sim.data$Y.list[[k]],sim.data$X.list[[k]],mod.oracle)
        fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = F)
        pow_err<-powErr(fdr.adj,sim.data$CL.list[[k]],gif,fdr)
        
        temp.oracle[[1]][[k]]<-list(cbind(scores,fdr.adj$pval,fdr.adj$qval))
        
        if(k == length(sim.data$X.list)){
          if(i==1){oracle.vals<-temp.oracle}
          else{oracle.vals<-append(oracle.vals,temp.oracle)}
        }
        
        if(exists("m.e")){
          m.e<-rbind(m.e,cbind(method="oracle",rep=k,et0=fdr.adj$param[3],
                               power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4]))}
        else{m.e<-cbind.data.frame(method="oracle",rep=k,et0=fdr.adj$param[3],
                                   power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4])}
        
      }
      #### glm ####
      if(glm.test == T){
        print("...fitting with glm...")
        
        p <- NULL
        z <- NULL
        for (l in 1:ncol(sim.data$Y.list[[k]])){
          mod.glm <- glm(sim.data$Y.list[[k]][,l] ~ ., 
                         data = data.frame(sim.data$X.list[[k]]),
                         binomial(link = "probit"))
          p[l] <- summary(mod.glm)$coeff[2,4]
          z[l] <- summary(mod.glm)$coeff[2,3] 
        }
        gif <- median(z^2)/0.456
        pcal <- pchisq(z^2/gif , df = 1, low = F)
        scores<-cbind(p,z,pcal)
       
        
        fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = F)
        pow_err<-powErr(fdr.adj,sim.data$CL.list[[k]],gif,fdr)
        
        temp.glm[[1]][[k]]<-list(cbind(scores,fdr.adj$pval,fdr.adj$qval))
        
        if(k == length(sim.data$X.list)){
          if(i==1){glm.vals<-temp.glm}
          else{glm.vals<-append(glm.vals,temp.glm)}
        }
        if(exists("m.e")){
          m.e<-rbind(m.e,cbind(method="glm",rep=k,et0=fdr.adj$param[3],
                               power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4]))}
        else{m.e<-cbind.data.frame(method="glm",rep=k,et0=fdr.adj$param[3],
                                   power_gif=pow_err[1],error_gif=pow_err[2],power_fdr=pow_err[3],error_fdr=pow_err[4])}
      }
 
      if(exists("eval.tab")){
        eval.tab<-rbind(eval.tab,cbind.data.frame(k.vals,m.e))
      }
      else{eval.tab<-cbind.data.frame(k.vals,m.e)}
    }
    
    if(exists("full.df")){
      full.df<-rbind(full.df,cbind.data.frame(sim,eval.tab))
    }
    else{full.df<-cbind.data.frame(sim,eval.tab)}
    
    print(paste("Sim ",i," of ",length(simRange)," complete",sep=""))
    print(paste("Project is ",(i/length(simRange))*100,"% complete."))
    print("===========================================================")
  }
  return(list(
    full.df=full.df,
    if(lfmm.test == T){lfmm.vals=lfmm.vals},
    if(cate.test == T){cate.vals = cate.vals},
    if(SVA.test == T){sva.vals = sva.vals},
    if(dSVA.test == T){dsva.vals = dsva.vals},
    if(RFEM.test == T){RFEM.vals = RFEM.vals},
    if(oracle.test == T){oracle.vals = oracle.vals},
    if(glm.test == T){glm.vals = glm.vals},
    val.dir=val.dir
  ))
}

#### Short Functions ####

simMeta<-function(x){
  len=length(strsplit(x,split = "/")[[1]])
  sim=paste(strsplit(x,split = "/")[[1]][len],"/",sep="")
  return(
    list(
      param=read.csv(paste(x,"SimulationParameterList.csv",sep="")),
      sim=sim
    )
  )
}

powErr <- function(x,y,gif,fdr){
  #Error gif
  ef.gif<-1-mean(which(-log10(x$qval) > gif) %in% y)
  #Power gif
  pf.gif<-1-mean(!(y %in% which(-log10(x$qval) > gif)))
  
  #error fdr
  ef.fdr<-1-mean(which((x$lfdr) < fdr) %in% y)
  #Power fdr
  pf.fdr<-1-mean(!(y %in% which((x$lfdr) < fdr)))
  
  return(c(pf.gif,ef.gif,pf.fdr,ef.fdr))
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
    mod.glm <- glm(data.y[,k] ~ ., 
                   data = data.frame(data.x,lf),
                   binomial(link = "probit"))
    p[k] <- summary(mod.glm)$coeff[2,4]
    z[k] <- summary(mod.glm)$coeff[2,3] 
  }
  gif <- median(z^2)/0.456
  pcal <- pchisq(z^2/gif , df = 1, low = F)
  return(cbind(p,z,pcal,gif))
}

simTestStore<-function(simTest,name=NULL,dir=NULL,verbose=T){
  index<-simTest$full.df %>% select(Sim,rep,method)
  sim.folder.out<-paste(dir,"Outputs/",sep="")
  sim.folder.loc<-paste(sim.folder.out,name,"/",sep="")
  dir.create(sim.folder.loc)
  methods<-length(unique(index$method))
  sims<-length(unique(index$Sim))
  reps<-length(unique(index$rep))
  write.csv(x=simTest$full.df,file = paste(sim.folder.loc,"SimTestList",".csv",sep=""))
  
  for(i in 2:methods){ # output from single method
    for(j in 1:sims){ # loop through sim
      for(l in 1:reps){ # loop through replicates
        r_temp<-data.frame(simTest[[i]][[j]][[l]])
        names(r_temp)<-c("p","z","pcal","fdr.p","fdr.q")
        #paste(unique(index$method)[1],"_",unique(index$Sim)[1],"_",unique(index$rep)[1],".csv",sep="")
        write.csv(x=r_temp,file = paste(sim.folder.loc,unique(index$method)[i-1],"_",unique(index$Sim)[j],"_",unique(index$rep)[l],".csv",sep=""))
      }
      if(verbose == T){print(paste("Storing...",((i-1)*sims+j)/(methods*sims)*100," % complete",sep=""))}
    }
  }
}









