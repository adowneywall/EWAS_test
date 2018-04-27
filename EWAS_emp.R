#### empirical data analysis

### libraries and scripts ####
source("EWAS_generator.R")
source("EWAS_methodTest.R")
library(lfmm)
library(missMDA)
library(fdrtool)
library(ggplot2)
### Reading in process empirical data ####

## Quercus
q.cpg<-readRDS("DATA/Empirical_Data/quercus_example/cpg_pruned.RData")
q.chh<-readRDS("DATA/Empirical_Data/quercus_example/chh_pruned.RData")
q.chg<-readRDS("DATA/Empirical_Data/quercus_example/chg_pruned.RData")
q.pred<-readRDS("DATA/Empirical_Data/quercus_example/predictors.RData")

## Baboon
b.cpg <- readRDS("DATA/Empirical_Data/baboon_example/baboon_prune.RData")
b.cpg.raw <- readRDS("DATA/Empirical_Data/baboon_example/baboon_raw.RData")
b.X<-read.table(file ="DATA/Empirical_Data/baboon_example/BSSeq_Baboon/predictor_n50.txt")
b.kin<-read.table(file="DATA/Empirical_Data/baboon_example/BSSeq_Baboon/relatedness_n50.txt")

# Multi Species Plot ####
q.cpg.means<-data.frame(x="Q-CpG",y=q.cpg$means)
q.chh.means<-data.frame(x="Q-CHH",y=q.chh$means)
q.chg.means<-data.frame(x="Q-CHG",y=q.chg$means)
baboon<-data.frame(x="Baboon",y=b.cpg$means)
species.means<-rbind(q.cpg.means,q.chh.means,q.chg.means,baboon)
ggplot(species.means,aes(x=x,y=y)) + geom_violin() + labs(y="Mean B-values",x="") +
  theme_bw()

# Perform PCA on different datasets (quercus and baboon so far) ####

## With pruned data
q.cpg.pca<-prcomp(q.cpg$B)
screeplot(q.cpg.pca,type='l')
plot(q.cpg.pca$sdev[2:10]^2,type="l")

q.chg.pca<-prcomp(q.chg$B)
screeplot(q.chg.pca)
plot(q.chg.pca$sdev[2:10]^2,type="l")

q.chh.pca<-prcomp(q.chh$B)
screeplot(q.chh.pca)
plot(q.chh.pca$sdev[2:10]^2,type="l")

b.cpg.pca<-prcomp(b.cpg$B)
screeplot(b.cpg.pca)
plot(b.cpg.pca$sdev[2:10]^2,type="l")

#Raw data
q.cpg.raw<-readRDS("DATA/Empirical_Data/quercus_example/cpg_raw.RData")[,6:63] # doesn't work with missing values
head(q.cpg.raw) # holes prevent pca
#q.cpg.raw.pca<-prcomp(q.cpg.raw)
#screeplot(q.cpg.raw.pca,type='l')


### All data sets appear to only have one pc of major variation


# Examing structure of variation in empirical dataset + data imputation ####
b.beta.raw<-b.cpg.raw$betaVal[,-c(1,2)]
b.mean.raw<-rowMeans(b.beta.raw,na.rm = T)
b.median.raw<-apply(b.beta.raw, 1, function(x){median(x,na.rm = T)})
b.sd.raw<-apply(b.beta.raw, 1,function(x){sd(x,na.rm=T)})
hist(b.mean.raw)
hist(b.median.raw)
hist(b.sd.raw)

b.raw.sum<-cbind.data.frame(mean=b.mean.raw,median=b.median.raw,sd=b.sd.raw)
ggplot(b.raw.sum,aes(x=median,y=sd)) + geom_point()

source("EWAS_generator.R")
multi.sim.gen(n=c(50),
              p=c(2000),
              K=c(1),
              freq=freq.multi,
              prop.variance = 0.4,#seq(from=0.1,to=0.9,by=.1),
              sigma=0.4,
              mean.B = c(1), 
              sd.U=b.pca.sd,
              sd.Um=0.2,
              sd.Usd=0.1,
              sd.V=0.5,
              rep=1,
              nb.t = nb.baboon,
              dir.name = "DATA/EWAS_Sims/", 
              sim.folder = "baboon.Rerun")

## Examining the performance of difference methods on simulated baboon data
sim.folder<-"baboon.test_2018-04-27/"
example<-simRead(SimRun = sMeta$sim,sim = sMeta$param,reps = 1)
sim.dat<-readRDS("DATA/EWAS_Sims/baboon.test2_2018-04-27/Sim1/Rep150_2000_0.4_1.Rdata")
sim.dat.Y<-sim.dat$Y
dim(sim.dat.Y)
lfmm_out<-lfmm::lfmm_ridge(Y = example$Y2.list[[1]],X = example$X.list[[1]],K =1)
pval<-lfmm_test(Y = example$Y.list[[1]],X = example$X.list[[1]],lfmm = lfmm_out,calibrate = "gif")

pca<-prcomp(example$Y.list[[1]])
screeplot(pca)
pval$gif
hist(pval$pvalue)
hist(pval$calibrated.pvalue)

sim.Y<-example$Y2.list[[1]]
sim.Yt<-t(sim.Y)
b.mean.sim<-rowMeans(sim.Yt,na.rm = T)
b.median.sim<-apply(sim.Yt,1, function(x){median(x,na.rm = T)})
b.sd.sim<-apply(sim.Yt,1,function(x){sd(x,na.rm=T)})
hist(b.mean.sim)
hist(b.median.sim)
hist(b.sd.sim,xlim=c(0,1))
b.sim.sum<-cbind.data.frame(mean=b.mean.sim,median=b.median.sim,sd=b.sd.sim)
ggplot(b.sim.sum,aes(x=mean,y=sd)) + geom_point()

# Baboon data with predictor ####
b.b<-b.cpg$B # create Y (methylation data)
b.meta<-b.cpg$metadata
b.b.X<-as.matrix(b.X) # X

### Model fitting (only run if necessary -  it will be slow) ####
library(lfmm)
library(fdrtool)
mod.lfmm <- lfmm::lfmm_ridge(t(b.b),b.b.X, K = 1) #lfmm estimate of latent factors
lfmmTest<-lfmm_test(t(b.b),b.b.X,lfmm=mod.lfmm)
l.fdr.adj<-fdrtool(as.vector(lfmmTest$score),statistic = c("normal"),verbose = F,plot = T)
l.fdr.adj$lfdr[l.fdr.adj$lfdr < 0.05]
# b.b[which(b.b > 1)] <- 1 # converts all methylation values over 1 to 1.
# scores<-binomal.glm.test(t(b.b),b.b.X,mod.lfmm$U) # glm with probit link function (beta distr.)
# fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = T) #for evaluating fdrtools
# lfmm.list<-list(scores=scores,fdrTool=fdr.adj)
# saveRDS(lfmm.list,"DATA/Empirical_Data/baboon_example/lfmm_BinomLogit_raw.Rdata")
lfmm.linear.list<-list(lfmm_Output=lfmmTest,fdrTool=l.fdr.adj)
saveRDS(lfmm.linear.list,"SampleData/Empirical_Data/baboon_example/lfmm_linear_raw.Rdata")

### Figures and Analysis ####
lfmm.sig<-readRDS("SampleData/Empirical_Data/baboon_example/lfmm_BinomLogit_raw.Rdata")
lfmm.linear<-readRDS("SampleData/Empirical_Data/baboon_example/lfmm_linear_raw.Rdata")

b.lfmm<-cbind.data.frame(b.meta,qval=lfmm.sig$fdrTool$qval)
b.lfmm$sig<-"Insignificant"
b.lfmm$sig[lfmm.sig$fdrTool$lfdr <= 0.05] <- "Significant"
ggplot(b.lfmm,aes(x=c(1:nrow(b.lfmm)),y=-log10(qval),colour=sig)) + geom_point()

linear.lfmm <- cbind.data.frame(b.meta,qval=lfmm.linear$qval)
linear.lfmm$sig <-  "Insiginificant"
linear.lfmm$sig[lfmm.linear$fdrTool$lfdrl <= 0.05] <-  "Significant"
ggplot(linear.lfmm,aes(x=c(1:nrow(linear.lfmm)),y=-log10(qval),colour=sig)) + geom_point()

# Quercus data w/ predictor ####
q.cpgM <-q.cpg$B # create Y (methylation data)
q.meta<-q.cpg$metadata
b.b.X<-as.matrix(b.X) #

### Model fitting (only run if necessary -  it will be slow) ####
# q.lfmm <- lfmm::lfmm_ridge(t(q.cpgM),q.pred$tmax, K = 1) #lfmm estimate of latent factors
# q.lfmmTest<-lfmm_test(t(b.b),b.b.X,lfmm=mod.lfmm)
# q.l.fdr.adj<-fdrtool(as.vector(lfmmTest$score),statistic = c("normal"),verbose = F,plot = T)
# #b.b[which(b.b > 1)] <- 1 # converts all methylation values over 1 to 1.
# scores<-binomal.glm.test(t(q.cpgM),q.pred$tmax,q.lfmm$U) # glm with probit link function (beta distr.)
# fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = T) #for evaluating fdrtools
# q.lfmm.list<-list(scores=scores,fdrTool=fdr.adj)
# saveRDS(q.lfmm.list,"DATA/Empirical_Data/quercus_example/lfmm_BinomLogit_raw.Rdata")
# q.lfmm.linear.list<-list(q.lfmm_Output=lfmmTest,fdrTool=q.l.fdr.adj)
# saveRDS(q.lfmm.list,"DATA/Empirical_Data/quercus_example/lfmm_linear_raw.Rdata")


### Figures and Analysis ####
q.lfmm.linear<-readRDS("DATA/Empirical_Data/quercus_example/lfmm_linear_raw.Rdata")
q.lfmm.logit<-readRDS("DATA/Empirical_Data/quercus_example/lfmm_BinomLogit_raw.Rdata")
  

b.lfmm<-cbind.data.frame(q.meta,qval=fdr.adj$qval)
b.lfmm$sig<-"Insignificant"
b.lfmm$sig[fdr.adj$lfdr <= 0.05] <- "Significant"
ggplot(b.lfmm,aes(x=c(1:nrow(b.lfmm)),
                  y=-log10(qval),
                  #group=interaction(chr,pos),
                  colour=sig)) + geom_point()
hist(fdr.adj$qval)
plot(-log10(b.lfmm$qval))
b.lfmm[7423,]
points(y=-log10(b.lfmm$qval[c(2443,2444,3551,7423)]),x=c(2443,2444,3551,7423),
       col="red",pch=16)


sc<-c("scaffold20751")
sc<-c("C2279863")
sc<-c("C2243561") # no smps
sc<-c("scaffold9024")
sc<-c("C2687893")
sc<-c("C21355867")
sc<-c("scaffold20751")
b.chr<-b.lfmm[which(b.lfmm$chr==sc),]

pos<-c("23308")
pos<-c("599")
pos<-c("395") #no smps
pos<-c("1977")
pos<-c("734")
b.chr[which(as.numeric(b.chr$pos)>=as.numeric(pos)-5 & as.numeric(b.chr$pos)<=as.numeric(pos)+5),]
match(q.cpg)

