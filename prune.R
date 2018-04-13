#### Pruning Data ####
#Simple script for pruning DNA methylation data sets (provided as N(row) x P(col) matrices)
#Will filter based on 4 primary criteria (arguements in the prune_data function):
#NAr = # of allowed nas for a single probe (default = 10)
#meanU or meanL = cut off threshold for mean probe values (default = 0.01 or 0.99)
#bsd = cut off for standard deviation for each probe (default = 0.01)

#Script can also handle some simple formating issues including metadata columns within the supplied matrix
#First set 'metadata = T', then provide the 'metaCols' arguement with a vector of cols to be removed 
#for the filtering step.
#Removed metadata will be stored in the resulting object as a list labels 'metadata'

#The two final arguement, violin and diagnostics, allows users to seem some basic diagnostic plots to see how
#the data is reduced during each filtering step. (default = T)

## Outputs:
# The function results in a object list with 6 elements (7 if metadata is provided).
# prune data: vector that provides pruning criteria (NAr,bmeanL,bmeanU,bsd) and the unfiltered and filtered loci counts.
# B: the nxp matrix of remaining after filtering 
# mean.vec: vector of mean b-values for each locus
# sd.vec: vector of standard deviations for each locus
# num.vec: vector of original positions for each locus

prune_data<-function(m,NAr=10,bmeanL=0.01,bmeanU=0.99,bsd=0.01,violin=T,metadata=F,metaCols,diagnostics=T){
  if(metadata){
    b<-as.matrix(m[,-metaCols])
    b.meta<-m[,metaCols]
  }
  #Store pruning criteria
  criteria<-c(NAr=NAr,bmeanL=bmeanL,bmeanU=bmeanU,bsd=bsd)
  #Basic Diagnostics
  # Taking probe means
  b.total.means<-rowMeans(b,na.rm=TRUE) #vector of means
  
  ##Remove fixed (0 or 1 beta values)
  rm.fixed<-!(b.total.means >= 1 | b.total.means <= 0)
  b<-b[rm.fixed,] # removes rows with no variation (1 or 0)
  b.meta<-b.meta[rm.fixed,] # removes rows with no variation (1 or 0)
  b.total.means<-b.total.means[rm.fixed]
  
  #Number,mean,sd
  b.total.num<-length(b.total.means)
  m.b.total.means<-mean(b.total.means)
  sd.b.total.means<-sd(b.total.means)
  
  ## Removing NAs (based on threshold defined by user)
  b.NAs<-apply(b,1,function(x){sum(is.na(x))})
  na.rm.index<-which(b.NAs <= NAr)
  b.NAr<-b[na.rm.index,]
  b.meta.NAr<-b.meta[na.rm.index,]
  b.means.NAr<-b.total.means[na.rm.index]#rowMeans(b.NAr,na.rm=TRUE)
  
  #Number,mean,sd
  b.NAr.num<-length(b.means.NAr)
  m.b.means.NAr<-mean(b.means.NAr)
  sd.b.means.NAr<-mean(b.means.NAr)
  
  ## Removed loci with low variation (set by user)
  l.mv.index<-which(b.means.NAr >= bmeanL & b.means.NAr <= bmeanU)
  b.lowVarRm<-b.NAr[l.mv.index,]
  b.meta.lowVarRm<-b.meta.NAr[l.mv.index,]
  bmeans.lowVarRm<-b.means.NAr[l.mv.index]
  
  #Number,mean,sd
  b.lowVarRm.num<-length(bmeans.lowVarRm)
  m.b.means.lowVarRm<-mean(bmeans.lowVarRm)
  sd.b.means.lowVarRm<-mean(bmeans.lowVarRm)
  
  # removes sd with greater than var of 0.01 (not currently doing this)
  b.sd<-apply(b.lowVarRm,1,function(x){sd(x,na.rm=T)})
  l.sdv.index<-which(b.sd >= bsd)
  b.lowsdRm<-b.lowVarRm[l.sdv.index,]
  b.meta.lowsdRm<-b.meta.lowVarRm[l.sdv.index,]
  bmeans.lowsdRm<-bmeans.lowVarRm[l.sdv.index]
  
  #Number,mean,sd
  b.lowsdRm.num<-length(bmeans.lowsdRm)
  m.b.means.lowsdRm<-mean(bmeans.lowsdRm)
  sd.b.means.lowsdRm<-mean(bmeans.lowsdRm)
  
  if(diagnostics){
    print("Total loci (minus fixed probes)....")
    print(paste("Total number of probes:",b.total.num))
    print(paste("Mean:",m.b.total.means))
    print(paste("Standard Deviation:",sd.b.total.means))
    hist(b.total.means)
    
    print("NAs removed....")
    print(paste("Total number of probes:",b.NAr.num))
    print(paste("Mean:",m.b.means.NAr))
    print(paste("Standard Deviation:",sd.b.means.NAr))
    hist(b.means.NAr)
    
    print("Mean low variation removed....")
    print(paste("Total number of probes:",b.lowVarRm.num))
    print(paste("Mean:",m.b.means.lowVarRm))
    print(paste("Standard Deviation:",sd.b.means.lowVarRm))
    hist(bmeans.lowVarRm)
    
    print(paste("Mean and SD low variation removed...."))
    print(paste("Total number of probes:",b.lowsdRm.num))
    print(paste("Mean:",m.b.means.lowsdRm))
    print(paste("Standard Deviation:",sd.b.means.lowsdRm))
    hist(bmeans.lowsdRm)
  }
  if(violin){
    violin<-data.frame(x="Total",y=b.total.means)
    violin<-rbind(violin,data.frame(x="NA.removed",y= b.means.NAr))
    violin<-rbind(violin,data.frame(x="Low Mean.r",y=bmeans.lowVarRm))
    violin<-rbind(violin,data.frame(x="Low SD.r",y=bmeans.lowsdRm))
    g<-ggplot(violin,aes(x=as.factor(x),y=y),environment = environment())
    print(g + geom_violin() + labs(x="Pruning Step",y="mean beta-value"))
  }
  pruneData<-c(criteria,Total_Loci=b.total.num,Remaining_Loci=b.lowsdRm.num)
  return(
    if(metadata){
      list(pruneData = pruneData, # vector of pruning criteria (nas removes, lower mean threshold, upper mean threshold, sd threshold)
           B=b.lowsdRm, # matrix of bvalues x samples
           means=bmeans.lowsdRm, # mean of bvalue across samples
           mean.vec=c(m.b.total.means,m.b.means.NAr,m.b.means.lowVarRm,m.b.means.lowsdRm), # vector of means for each pruning step
           sd.vec=c(sd.b.total.means,sd.b.means.NAr,sd.b.means.lowVarRm,sd.b.means.lowsdRm), # vector of sd for each pruning step
           num.vec=c(b.total.num,b.NAr.num,b.lowVarRm.num,b.lowsdRm.num), #Vector of number of probes after each pruning step
           metadata=b.meta.lowsdRm) #matrix of meta data for each loci after pruning
    }
    else{
      list(pruneData = pruneData, # vector of pruning criteria (nas removes, lower mean threshold, upper mean threshold, sd threshold)
           B=b.lowsdRm,
           means=bmeans.lowsdRm,
           mean.vec=c(m.b.total.means,m.b.means.NAr,m.b.means.lowVarRm,m.b.means.lowsdRm), # vector of means for each pruning step
           sd.vec=c(sd.b.total.means,sd.b.means.NAr,sd.b.means.lowVarRm,sd.b.means.lowsdRm), # vector of sd for each pruning step
           num.vec=c(b.total.num,b.NAr.num,b.lowVarRm.num,b.lowsdRm.num) #Vector of number of probes after each pruning step
      )
    }
  )
}