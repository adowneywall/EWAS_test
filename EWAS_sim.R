
#### libraries and scripts ####
library(ggplot2)
library(cate)
library(dplyr)
library(reshape2)
library(dSVA)
library(sva)
library(fdrtool)
library(RefFreeEWAS)

source("EWAS_generator.R")
source("EWAS_methodTest.R")

#install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock'))

### Reading in process empirical data ####
q.cpg<-readRDS("DATA/Empirical_Data/quercus_example/cpg_pruned.RData")
q.chh<-readRDS("DATA/Empirical_Data/quercus_example/chh_pruned.RData")
q.chg<-readRDS("DATA/Empirical_Data/quercus_example/chg_pruned.RData")
b.cpg <- readRDS("DATA/Empirical_Data/baboon_example/baboon_prune.RData")

q.cpg.means<-data.frame(x="Q-CpG",y=q.cpg$means)
q.chh.means<-data.frame(x="Q-CHH",y=q.chh$means)
q.chg.means<-data.frame(x="Q-CHG",y=q.chg$means)
baboon<-data.frame(x="Baboon",y=b.cpg$means)
species.means<-rbind(q.cpg.means,q.chh.means,q.chg.means,baboon)
ggplot(species.means,aes(x=x,y=y)) + geom_violin() + labs(y="Mean B-values",x="") +
  theme_bw()

#### Single simulation Tests ####

## Simulation Arguments
N.NUM = 50 # number of individuals
P.NUM = 2000 # number of loci
sample<-q.cpg$B
freq.e<-species.means[species.means$x == "Q-CpG",]

#### Single Dataset Simulation
# With empirically derived b-val means
sim_empiricalMean<-ewas_generator(N.NUM,P.NUM,5,sd.U = 1,sd.V=1,sigma=1,freq = sample(x = freq.e,size = P.NUM ,replace=T)) 
# random means from uniform distribution
sim_randMean<-ewas_generator(N.NUM,P.NUM,5,sd.U = 1,sd.V=1,sigma=1,freq = NULL)

## Visualling b-value distribution from single sample simulations
# Ploting histograms for loci between N and M
par(mfrow = c(2,3))
N<-2 # Start
M<-100 # Stop
for(i in N:M){
  locus = sample(which(freq.e %in% sim_empiricalMean$freq[i]),1)
  locus.bvalues = as.numeric(sample[locus,!is.na(sample[locus,])])
  
  hist(sim_empiricalMean$Y[,i],xlim = c(0,1),breaks=20,main=paste("Empirical Gen. Mean. Loci ",i,sep=""))
  hist(locus.bvalues,xlim = c(0,1),breaks=20,main=paste("Empiricial Data. Loci ",i,sep=""))
  hist(sim1_randMean$Y[,i],xlim = c(0,1),breaks=20,main=paste("Random Mean. Loci ",i,sep=""))
  
  hist(sim_empiricalMean$Y.collapse[,i],xlim = c(0,1),breaks=20,main=paste("Empiricial Gen. Mean w/ collapsed beta-value. Loci ",i,sep=""))
  hist(locus.bvalues,xlim = c(0,1),breaks=20,main=paste("Locus",locus,", bval.mean =",sim_empiricalMean$freq[i]))
  hist(sim1_randMean$Y.collapse[,i],xlim = c(0,1),breaks=20,main=paste("Random Mean w/ collapsed beta-value. Loci ",i,sep=""))
}

### Multi simulation Tests ####
freq.multi<-b.cpg$means

## getting real variance to parameterize pcs
#baboon data
b.pca<-prcomp(b.cpg$B)
screeplot(b.pca) # really only one pc here
b.pca.sd<-b.pca$sdev[1] # possible include more if you want to include more pcs in sims

multi.sim.gen(n=c(250),
              p=c(2000),
              K=c(1),
              freq=freq.multi,
              prop.variance = seq(from=0.1,to=0.9,by=.1),
              sigma=0.1,
              mean.B = c(1,2,5), 
              sd.U=b.pca.sd,
              sd.Um=0.2,
              sd.Usd=0.2,
              rep=100, 
              dir.name = "DATA/EWAS_Sims/", 
              sim.folder = "baboon_Sim01")

## Examining the performance of difference methods on simulated baboon data
sim.folder<-"DATA/EWAS_Sims/baboon_Sim01_2018-04-12/"

sMeta<-simMeta(sim.folder) # enter the name of folder w/ directory, returns parameters list [param] and sim folder name

# Options for selecting a subset of your simulation list for testing
sub.list<-subset(sMeta$param,sMeta$param$n >= 100 & sMeta$param$p >= 2000)
#simRange<-seq(from=1,to=max(sim.list$X))
simRange<- sMeta$param$Sim #sub.list$Sim

baboon.Sim01<-simTest(sMeta$sim,simRange,sMeta$param,K.est = T,lfmm.test = T,cate.test = T,RFEM.test = F,SVA.test = T,dSVA.test = T,
        oracle.test = T,glm.test = T,rep.times = 2)
## currently K.est needs to = T, otherwise script won't run

simTestStore(t4,name="baboon_Sim01",dir=sim.folder)

## Initial analysis of simulation data
ss<-t5$full.df
ss.sub <- subset(ss,ss$n == 250)
ss.sub <- subset(ss.sub,ss.sub$p == 5000)
library(ggplot2)
ggplot(ss.sub,aes(x=method,y=as.numeric(et0))) + geom_boxplot()
ggplot(ss.sub,aes(x=method,y=as.numeric(power_fdr))) + geom_boxplot()
ggplot(ss.sub,aes(x=method,y=as.numeric(error_fdr))) + geom_boxplot()


sub.lfmm <- subset(ss,ss$method == "lfmm")
sub.lfmm <- subset(sub.lfmm,sub.lfmm$p == 5000)
sub.lfmm <- subset(sub.lfmm,sub.lfmm$n == 250)
ggplot(sub.lfmm,aes(x=as.factor(prop.variance),y=as.numeric(et0))) + geom_boxplot()
ggplot(sub.lfmm,aes(x=as.factor(prop.variance),y=as.numeric(power_fdr))) + geom_boxplot()

mean(t1$full.df$K.bcv)
mean(t1$full.df$K.be)


P<-1000
N<-500
log.test<-ewas_generator(n = N,p = P,K = 5,freq = sample(q.cpg.means$y,size = P,replace=T),
                         prop.causal = 0.025,prop.variance = 0.5,
                         sd.U = .2,sd.V = .5,sd.B = .2,mean.B = 2,sigma = .2/5 
                        )
hist(log.test$Y[,10])
hist(log.test$Y.logit[,15])

dim(q.cpg$B[,10])
bb<-as.matrix(q.cpg$B)
dim(bb)
hist(bb[2000,])
I<-6000
hist(log(bb[I,]/(1-bb[I,])))

