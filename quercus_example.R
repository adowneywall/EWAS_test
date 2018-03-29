
library(readr)
library(data.table)
library(dplyr)
library(ggplot2)
#q.meth<-read_table2("DATA/Empirical_Data/quercus_example/set3.meth.full.tab.gz",col_names = T)
rm(q.meth)
q_meth.start <- read_table2("DATA/Empirical_Data/quercus_example/set3.meth.full.tab.gz", 
                      col_names = TRUE, col_types = cols(.default = "d",
                                                          chr = "c",
                                                         watson = "c",
                                                         context = "c",
                                                         dinucl = "c"), 
                      n_max =500000 )
cn<-colnames(q_meth.start)
q_meth.add1 <- read_table2("DATA/Empirical_Data/quercus_example/set3.meth.full.tab.gz", 
                       col_names = FALSE, col_types = cols(.default = "d",
                                                          X1 = "c",
                                                          X2 = "c",
                                                          X4 = "c",
                                                          X5 = "c"), 
                       skip = 500001, n_max =500000 )
q_meth.add2 <- read_table2("DATA/Empirical_Data/quercus_example/set3.meth.full.tab.gz", 
                           col_names = FALSE, col_types = cols(.default = "d",
                                                               X1 = "c",
                                                               X2 = "c",
                                                               X4 = "c",
                                                               X5 = "c"), 
                           skip = 1000001, n_max = 500000)
colnames(q_meth.add1)<-cn
colnames(q_meth.add2)<-cn

#CHG
q.meth.s1<-subset(q_meth.start,q_meth.start$context == "CHG")
q.meth.s2<-subset(q_meth.add1,q_meth.add1$context == "CHG")
q.meth.s3<-subset(q_meth.add2,q_meth.add2$context == "CHG")
q.chg<-bind_rows(q.meth.s1,q.meth.s2,q.meth.s3)


#CHH
q.meth.s1<-subset(q_meth.start,q_meth.start$context == "CHH")
q.meth.s2<-subset(q_meth.add1,q_meth.add1$context == "CHH")
q.meth.s3<-subset(q_meth.add2,q_meth.add2$context == "CHH")
q.chh<-bind_rows(q.meth.s1,q.meth.s2,q.meth.s3)

# CpGs
q.meth.s1<-subset(q_meth.start,q_meth.start$context == "CG")
q.meth.s2<-subset(q_meth.add1,q_meth.add1$context == "CG")
q.meth.s3<-subset(q_meth.add2,q_meth.add2$context == "CG")
q.cg<-bind_rows(q.meth.s1,q.meth.s2,q.meth.s3)

#removing large dataset from global environment
rm(q.meth.s1,q.meth.s2,q.meth.s3)
rm(q_meth.start,q_meth.add1,q_meth.add2)



write.csv(file="DATA/Empirical_Data/quercus_example/quercus_meanbvalues.csv",adj.quercus.means,row.names = F)
write.csv(file="DATA/Empirical_Data/quercus_example/processed_CpG_bvalues.csv",q.cg.lowvarOmit,row.names = F)


#output_cpg<-prune_data(q.cg,metadata = T,metaCols = c(1:5),diagnostics = F)
output_cpg<-prune_data(q.cg,metadata=T,metaCols=c(1:5),bmeanL = .01,bmeanU = .99,NAr = 0)
output_chh<-prune_data(q.chh,metadata=T,metaCols=c(1:5),bmeanL = .01,bmeanU = .99,NAr = 0)
output_chg<-prune_data(q.chg,metadata=T,metaCols=c(1:5),bmeanL = .01,bmeanU = .99,NAr = 0)


q.cpg.means<-data.frame(x="Q-CpG",y=output_cpg$means)
q.chh.means<-data.frame(x="Q-CHH",y=output_chh$means)
q.chg.means<-data.frame(x="Q-CHG",y=output_chg$means)
baboon<-data.frame(x="baboon",y=output_baboon$means)
species.means<-rbind(q.cpg.means,q.chh.means,q.chg.means,baboon)
ggplot(species.means,aes(x=x,y=y)) + geom_violin()


prune_data<-function(m,NAr=10,bmeanL=0.01,bmeanU=0.99,bsd=0.01,violin=T,metadata=F,metaCols,diagnostics=T){
  if(metadata){
    b<-as.matrix(m[,-metaCols])
    b.meta<-m[,metaCols]
  } 
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
  return(
    if(metadata){
      list(B=b.lowsdRm, # matrix of bvalues x samples
           means=bmeans.lowsdRm, # mean of bvalue across samples
           mean.vec=c(m.b.total.means,m.b.means.NAr,m.b.means.lowVarRm,m.b.means.lowsdRm), # vector of means for each pruning step
           sd.vec=c(sd.b.total.means,sd.b.means.NAr,sd.b.means.lowVarRm,sd.b.means.lowsdRm), # vector of sd for each pruning step
           num.vec=c(b.total.num,b.NAr.num,b.lowVarRm.num,b.lowsdRm.num), #Vector of number of probes after each pruning step
           metadata=b.meta.lowsdRm) #matrix of meta data for each loci after pruning
    }
    else{
      list(B=b.lowsdRm,
           means=bmeans.lowsdRm,
           mean.vec=c(m.b.total.means,m.b.means.NAr,m.b.means.lowVarRm,m.b.means.lowsdRm), # vector of means for each pruning step
           sd.vec=c(sd.b.total.means,sd.b.means.NAr,sd.b.means.lowVarRm,sd.b.means.lowsdRm), # vector of sd for each pruning step
           num.vec=c(b.total.num,b.NAr.num,b.lowVarRm.num,b.lowsdRm.num) #Vector of number of probes after each pruning step
           )
    }
  )
}


write.csv(file="DATA/Empirical_Data/quercus_example/processed_CpG_bvalues.csv",q.cg.lowvarOmit,row.names = F)


#Extra #####
q.cg.mat<-as.matrix(q.cg[,6:63])

#Basic Diagnostics

# Taking probe means
q.cg.means<-rowMeans(q.cg.mat,na.rm=TRUE) # vector of means

# Order matters here be careful
q.cg.mat<-q.cg.mat[!(q.cg.means >= 1 | q.cg.means <= 0),] # removes rows with no variation (1 or 0)
q.cg.means<-q.cg.means[!(q.cg.means >= 1 | q.cg.means <= 0)]

mean(q.cg.means)
hist(abs(q.cg.means-0.5))
hist(q.cg.means)
violin<-data.frame(x=1,y=q.cg.means)
ggplot(violin,aes(x=x,y=y)) + geom_violin() + geom_boxplot(width=0.02)

#Taking probe SD
q.cg.sd<-apply(q.cg.mat,1,function(x){sd(x)})
hist(q.cg.sd)
paste(dim(q.cg.mat)[1]/length(q.cg$chr)*100, "% of CpG sites are variable for a total of", 
      dim(q.cg.mat)[1], "CpG probes with mean of",mean(q.cg.means)) # 261652 probes left

#Removed NAs (at threshold of 10)
q.cg.NAs<-apply(q.cg.mat,1,function(x){sum(is.na(x))})
hist(q.cg.NAs)
NA.thres<-0
q.cg.NAomit<-q.cg.mat[which(q.cg.NAs <= NA.thres),]
hold<-apply(q.cg.NAomit,1,function(x){sum(is.na(x))})
hist(hold)
q.cg.means.lowNA<-rowMeans(q.cg.NAomit,na.rm=TRUE)
hist(q.cg.means.lowNA)
mean(q.cg.means.lowNA)
paste(dim(q.cg.NAomit)[1]," probes remaing after NA pruning") # 73224 probes left
paste((1-dim(q.cg.NAomit)[1]/dim(q.cg.mat)[1])*100, "% of probes removed with NA threshold of", NA.thres)
violin<-rbind(violin,data.frame(x=2,y=q.cg.means.lowNA))
ggplot(violin,aes(x=as.factor(x),y=y)) + geom_violin(width=1) #+ geom_boxplot(width=0.02)


## Removing SNPs with little variation
# removes means less than 0.01 and greater than 0.09
UVT<-0.99
LVT<-0.01 
q.cg.lowvarOmit<-q.cg.NAomit[which(q.cg.means.lowNA >= LVT & q.cg.means.lowNA <= UVT),]
# hist of mean prob vals after very high and low var removed
adj.quercus.means<-q.cg.means.lowNA[which(q.cg.means.lowNA >= LVT & q.cg.means.lowNA <= UVT)]
write.csv(file="DATA/Empirical_Data/quercus_example/quercus_meanbvalues.csv",adj.quercus.means,row.names = F)
hist(adj.quercus.means)
hist(sample(x =adj.quercus.means,size =1000 ,replace=T))
# mean of mean prob value after v. high and low var removed
mean(q.cg.means.lowNA[which(q.cg.means.lowNA >= LVT & q.cg.means.lowNA <= UVT)])
paste(dim(q.cg.lowvarOmit)[1],"probes remaining after probes with low variation removed") #12980 probes remaining
hist(q.cg.lowvarOmit)
violin<-rbind(violin,data.frame(x=3,y=adj.quercus.means))
ggplot(violin,aes(x=as.factor(x),y=y)) + geom_violin(width=1) #+ geom_boxplot(width=0.02)

# removes sd with greater than var of 0.01 (not currently doing this)
q.cg.sd<-apply(q.cg.lowvarOmit,1,function(x){sd(x,na.rm=T)})
#q.cg.lowvarOmit<-q.cg.lowvarOmit[which(q.cg.sd > 0.01),]
hist(apply(q.cg.lowvarOmit,1,function(x){sd(x,na.rm=T)}))
plot(density(q.cg.lowvarOmit[1,],na.rm=T))
write.csv(file="DATA/Empirical_Data/quercus_example/processed_CpG_bvalues.csv",q.cg.lowvarOmit,row.names = F)











#Diagnostics ####
# Ploting histograms for loci between N and M
N<-2
M<-100
for(i in N:M){
  #plot(density(q.cg.lowvarOmit[i,],na.rm=T),main=paste(i))
  hist(q.cg.lowvarOmit[i,],na.rm=T,main=paste(i))
}

library(MASS)
tt<-q.cg.lowvarOmit[30,]
tt[tt == 0.00] <- 0.000001
tt[tt == 1.0] <- 0.9999999
tt.beta<-fitdistr(tt,dbeta,start=list(shape1=1,shape2=2))
tt.gamma<-fitdistr(tt,dgamma,list(shape = 1, rate = 0.1))
tt.weibull<-fitdistr(tt, "weibull")
tt.beta
tt.gamma
tt.weibull

hist(tt,breaks=100,freq = F)
lines(density((tt)))
tt.b<-rbeta(n = 10000, tt.beta$estimate[1],tt.beta$estimate[2])
tt.g<-rbeta(n = 10000, tt.gamma$estimate[1],tt.gamma$estimate[2])
tt.w<-rweibull(n = 10000, tt.weibull$estimate[1],tt.weibull$estimate[2])
lines(density(tt.b),col="red")
lines(density(tt.g),col="green")
lines(density(tt.w),col="blue")

## Suggesting probe beta-values are not normal...
tt.resample<-sample(tt.b,1000,replace = T)
hist(tt.resample,breaks=100)
# libraries

library(nortest)
library(fitdistrplus)


shapiro.test(q.cg.lowvarOmit[8,]) #not normal
pearson.test(q.cg.lowvarOmit[8,]) #not normal

# Creat data fram of probe means and sds
quercus.var.df<-data.frame(mean=adj.quercus.means,sd=q.cg.sd)

## Variance distribution across mean bins w/ quercus dataset
bins<-seq(from=0,to=1,by=(1/20))
sd.means<-NULL
for(i in 1: (length(bins)-1)){
  print(i)
  temp.b<-which(quercus.var.df$mean>bins[i] & quercus.var.df$mean<=bins[i+1])
  sd.means<-c(sd.means,mean(quercus.var.df$sd[temp.b]))
}
plot(sd.means~c(as.factor(paste(bins[1:20],"-",bins[2:21]))))
hist(quercus.var.df$mean,breaks=20)
hist(quercus.var.df$sd)



## Variance distributino across mean bins w/ simulation data
sim1<-ewas_generator(50,2000,5,
                     freq = sample(x = adj.quercus.means,size = 2000 ,replace=T))
sim1<-ewas_generator(50,2000,5,sd.U = 0.1,sd.V=0.1,sigma=.3) #freq = sample(x = adj.quercus.means,size = 2000 ,replace=T)
hist(adj.quercus.means)

bins<-seq(from=0,to=1,by=(1/20))
sd.means<-NULL

N<-2
M<-100
for(i in N:M){
  #plot(density(q.cg.lowvarOmit[i,],na.rm=T),main=paste(i))
  hist(sim1$Y.collapse[i,],na.rm=T,main=paste(i))
}

hist(sim1$Y[5,])
hist(sim1$Y.collapse[5,])
max(sim1$Y.collapse[5,])
mean.sim<-rowMeans(t(sim1$Y))
sd.sim<-apply(t(sim1$Y),1,function(x){sd(x,na.rm=TRUE)})
sim.df<-data.frame(mean=mean.sim,sd=sd.sim)

for(i in 1:(length(bins)-1)){
  print(i)
  temp.b<-which(sim.df$mean>bins[i] & sim.df$mean<=bins[i+1])
  sd.means<-c(sd.means,mean(sim.df[temp.b,2]))
  }
plot(sd.means~as.factor(bins[2:21]))
hist(sim.df$mean,breaks=20)     
hist(sim.df$sd)


# distribution of loci means
#quercus
hist(sample(x = adj.quercus.means,size = 2000 ,replace=T))
#baboon
hist(mean.chrom) # based on all 20 chromosoms
#sims
hist(sim.df$mean,breaks=20)



range<-sample(1:nrow(quercus.var.df),1000,replace = T)
freq<-sample_n(quercus.var.df,10,replace = T)
freq[1,1]
freq$mean

hist(sim1$M)
dim(sim1$Epsilon)
sim1$Epsilon[1:10,1:10]
hist(sim1$freq)

hold<-sim1$freq[1]+sim1$Z[3,]
hist(sim1$Z[2,])
mean(hold)
mean(sim1$Z[2,])

