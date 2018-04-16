#### empirical data analysis

### Reading in process empirical data ####

## Quercus
q.cpg<-readRDS("DATA/Empirical_Data/quercus_example/cpg_pruned.RData")
q.chh<-readRDS("DATA/Empirical_Data/quercus_example/chh_pruned.RData")
q.chg<-readRDS("DATA/Empirical_Data/quercus_example/chg_pruned.RData")

## Baboon
b.cpg <- readRDS("DATA/Empirical_Data/baboon_example/baboon_prune.RData")
b.X<-read.table(file ="DATA/Empirical_Data/baboon_example/BSSeq_Baboon/predictor_n50.txt")
b.kin<-read.table(file="DATA/Empirical_Data/baboon_example/BSSeq_Baboon/relatedness_n50.txt")

## Multi Species Plot
q.cpg.means<-data.frame(x="Q-CpG",y=q.cpg$means)
q.chh.means<-data.frame(x="Q-CHH",y=q.chh$means)
q.chg.means<-data.frame(x="Q-CHG",y=q.chg$means)
baboon<-data.frame(x="Baboon",y=b.cpg$means)
species.means<-rbind(q.cpg.means,q.chh.means,q.chg.means,baboon)
ggplot(species.means,aes(x=x,y=y)) + geom_violin() + labs(y="Mean B-values",x="") +
  theme_bw()

### Perform PCA on different datasets (quercus and baboon so far)

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

#With Raw data
q.cpg.raw<-readRDS("DATA/Empirical_Data/quercus_example/cpg_raw.RData")[,6:63] # doesn't work with missing values
head(q.cpg.raw)
q.cpg.raw.pca<-prcomp(q.cpg.raw)
screeplot(q.cpg.raw.pca,type='l')
plot(q.cpg.pca$sdev[2:10]^2,type="l")


### All data sets appear to only have one pc of major variation


b.b<-b.cpg$B # create Y (methylation data)
b.meta<-b.cpg$metadata
b.b.X<-as.matrix(b.X) # X
library(lfmm)
mod.lfmm <- lfmm::lfmm_ridge(t(b.b),b.b.X, K = 1) #lfmm estimate of latent factors
b.b[which(b.b > 1)] <- 1 # converts all methylation values over 1 to 1.
scores<-binomal.glm.test(t(b.b),b.b.X,mod.lfmm$U) # glm with probit link function (beta distr.)
fdr.adj<-fdrtool(scores[,2],statistic = c("normal"),verbose = F,plot = T) #for evaluating fdrtools
hist(fdr.adj$lfdr)

b.lfmm<-cbind.data.frame(b.meta,qval=fdr.adj$qval)
ggplot(b.lfmm,aes(x=as.numeric(chrom.loc),
                  y=-log10(qval),
                  group=chrom,
                  colour=chrom)) + geom_point()

plot(-log10(fdr.adj$qval))
points(-log10(fdr.adj$qval)[b.ord]~b.ord,
       col="red",pch=16)



