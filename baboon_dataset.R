

#devtools::install_github("bcm-uga/lfmm")
library(cate)
library(lfmm)
library(Rarity)

tot.count<-read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/counts_chr1_n50.txt",header = T,row.names = 1)
m.count<-read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/mcounts_chr1_n50.txt",header=T,row.names = 1)
p.variable<-read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/predictor_n50.txt")
relate.mat<-read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/relatedness_n50.txt")
cov.mat<-read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/covariates_n50.txt")
dim(tot.count)
dim(m.count)
dim(p.variable)
dim(relate.mat)
dim(cov.mat)


## beta-values and trimming data ####
beta<-m.count/tot.count
beta.mat<-as.matrix(beta)

write.table(x = beta.mat, file = "DATA/Empirical_Data/baboon_example/BSSeq_Baboon/betaval_chr1_n50.txt",sep="")
## Useful function for finding rowxcol of nans or infinites
#which(is.infinite(beta.mat.NOz), arr.ind=TRUE)


## Currently removes all CpGs with any NAs - this may need to be revised!@!
proportion.NA<-function(x){sum(is.na(x))/50}
testing.NA<-apply(beta.mat,1,proportion.NA)
hist(testing.NA)
beta.mat.NONA<-beta.mat[testing.NA == 0.00,]
dim(beta.mat.NONA)
hist(beta.mat.NONA,rm.na=T)
# Removing non sensical beta-values (those above 1)
max.list<-apply(beta.mat.NONA,1,function(x){max(x)})
beta.mat.NONA<-beta.mat.NONA[-c(which(max.list > 1)),]
hist(beta.mat.NONA)

N<-1
M<-100
for(i in N:M){
  #plot(density(q.cg.lowvarOmit[i,],na.rm=T),main=paste(i))
  hist(beta.mat.NONA[i,],main=paste(i))
}


## Removes CpGs with high (>0.9) or low (<0.1) mean methylation
testing.mean<-apply(beta.mat.NONA,1,mean)
hist(testing.mean)
length((which(testing.mean <= 0.05 | testing.mean >= 0.95 )))/length(testing.mean) # removing approximate 2.5% of loci

mean(testing.mean)
# mean of chr 1 of baboon data = 0.6064
sd(testing.mean) 
# sd of chr 1 of baboon data = .2624

beta.mat.NONA.var<-beta.mat.NONA[testing.mean > 0.1 & testing.mean < 0.9,]
dim(beta.mat.NONA.var)

## Remove CpGs with minimal variation (=<5% of overall dataset)
testing.sd<-apply(beta.mat.NONA.var,1,sd)
beta.mat.NONA.varlow<-beta.mat.NONA.var[testing.sd>quantile(testing.sd,probs = c(0.05)),]
dim(beta.mat.NONA.varlow)

dim(beta.mat.NONA.varlow)
dim(sim1$Y)

bins<-seq(from=0,to=1,by=(1/20))
sd.means<-NULL
mean.sim<-rowMeans(beta.mat.NONA.varlow)
sd.sim<-apply(beta.mat.NONA.varlow,1,function(x){sd(x,na.rm=TRUE)})
baboon.df<-data.frame(mean=mean.sim,sd=sd.sim)

for(i in 1:(length(bins)-1)){
  print(i)
  temp.b<-which(baboon.df$mean>bins[i] & baboon.df$mean<=bins[i+1])
  sd.means<-c(sd.means,mean(baboon.df[temp.b,2]))
}
plot(sd.means~as.factor(bins[2:21]))
hist(baboon.df$mean,breaks=20)     
hist(baboon.df$sd)

#### Estimate confounders (K value)####

# Basic PCA w/ scree plot 
dim(t(beta.mat.NONA.varlow))
pca_basic<-prcomp(t(beta.mat.NONA.varlow))
screeplot(pca_basic,npcs = 20)

## Based on cate package
est.confounder.num(~ V1, p.variable, t(beta.mat.NONA.varlow), method = "ed")
est.confounder.num(~ V1, p.variable, t(beta.mat.NONA.varlow), method = "bcv")

### Estimating confounders based on sva package
num_sv     <- sva::num.sv(dat = beta.mat.NONA.varlow, mod = p.variable$V1, method = "be")
num_sv_l   <- sva::num.sv(dat = beta.mat.NONA.varlow, mod = p.variable$V1, method = "leek")

# K ####
K<-2
#### lfmm ####
mod <- lfmm::lfmm_ridge(t(beta.mat.NONA.varlow),p.variable$V1,K = K)
## tests based on GLMs -- replaces lfmm_testing() function in lfmm package due to the binomial nature of data


## Covariance examination between predictor and estimated latent factors ####
## Cov matrix based on lfmm latent factors and predictor
covariance.mat<-cov(cbind(mod$U,p.variable$V1))
cov2cor(covariance.mat)
corPlot(cbind(mod$U,p.variable$V1),method = "spearman")

## Total correlation (r2) between predictor and estimated latent factors
summary(lm(p.variable$V1 ~ mod$U))

## Calculating p-values using glm with probit link function - based on lfmm latent factor estimation ####
p.value <- NULL
z.score <- NULL
effect.mat <- matrix(ncol=K+1,nrow=nrow(beta.mat.NONA.varlow)) #
error.loci<-NULL
for (j in 1:nrow(beta.mat.NONA.varlow)){ #nrow(beta.mat.NONA.varlow)
  if(max(t(beta.mat.NONA.varlow)[,j])<=1){
    mod.glm <- glm(t(beta.mat.NONA.varlow)[,j] ~ ., 
                   data = data.frame(p.variable$V1, mod$U),
                   binomial(link = "probit"))
    p.value[j] <- summary(mod.glm)$coeff[2,4]
    z.score[j] <- summary(mod.glm)$coeff[2,3]
    effect.mat[j,] <- summary(mod.glm)$coeff[2:(K+2),1]
  }
  else{
    error.loci<-c(error.loci,j)
  }
}
p.value.corr<-p.value[!is.na(p.value)]
z.score.corr<-z.score[!is.na(z.score)]
effect.mat.corr<-effect.mat[complete.cases(effect.mat),]

### Examining the effect sizes
library(lattice)
library(reshape2)

# create dataframe with effect sizes for predictor + latent factors
e.m<-data.frame(effect.mat.corr)
colnames(e.m)<-c("Predictor",paste("L",seq(from=1,to=K,by=1),sep=""))

##Plot of predictor effect sizes (shown as abs)

# effect size distribution associated with predictor
densityplot(~abs(Predictor),data=e.m,auto.key = T) 

# effect size distribution associated with predictor and all latent factors
densityplot(~abs(e.m$Predictor)+abs(unlist(e.m[,2:K+1])),auto.key = T) 

# Plot of predictors and all latent factors plotted separately
e.melt<-melt(e.m)
densityplot(~ abs(value), groups = variable, data = e.melt, auto.key = TRUE)

## Basic effect size stats
mean.effect<-apply(effect.mat.corr,2,mean)
sd.effect<-apply(effect.mat.corr,2,sd)

### Calculating p-values based on corrected z-scores w/ genomic inflation factor (gif)
gif <- median(z.score.corr^2)/0.456
p.values.calibrated <- pchisq(z.score.corr^2/gif , df = 1, low = F)
plot(-log10(p.values.calibrated)~c(1:length(p.values.calibrated)),xlab="Loci Position")
#order(simu$causal)
#simu$causal[1]
#points(-log10(p.values.calibrated)[simu$causal]~simu$causal,
#       col="red",pch=16)
points(-log10(p.values.calibrated)[which(-log10(p.values.calibrated) > 3)]~which(-log10(p.values.calibrated) > 3),
       col="blue",pch=21,cex=2)
abline(h = 3,col="red")

#### cate ####
output.EWAS.cate <- cate(~ p.variable$V1,p.variable, t(beta.mat.NONA.varlow), r = 5,adj.method = "rr")

## tests based on GLMs
p.value <- NULL
z.score <- NULL

for (j in 1:nrow(beta.mat.NONA.varlow)){

  if(max(t(beta.mat.NONA.varlow)[,j])<=1){
    mod.glm <- glm(t(beta.mat.NONA.varlow)[,j] ~ ., 
                   data = data.frame(p.variable$V1, output.EWAS.cate$Z),
                   binomial(link = "probit"))
    p.value[j] <- summary(mod.glm)$coeff[2,4]
    z.score[j] <- summary(mod.glm)$coeff[2,3]
    effect.mat[j,] <- summary(mod.glm)$coeff[2:(K+2),1]
  }
  else{
    error.loci<-c(error.loci,j)
  }
}

## Using FDR tool
library(fdrtool)
(fdr_cate_output<-fdrtool(z.score,statistic = c("normal"))) 
fdr_cate_output$pval
plot(-log10(fdr_cate_output$pval)~c(1:nrow(beta.mat.NONA.varlow)),xlab="Loci Position")
#points(-log10(p.values.calibrated.cate)[simu$causal]~simu$causal,
#       col="red",pch=16)
points(-log10(fdr_cate_output$pval)[which(-log10(fdr_cate_output$pval) > 4)]~which(-log10(fdr_cate_output$pval) > 4),
       col="blue",pch=21,cex=2)
abline(h = 4,col="red")

#Using standard GIF
gif <- median(z.score^2)/0.456 #genomic inflation factor
p.values.calibrated.cate <- pchisq(z.score^2/gif , df = 1, low = F)
plot(-log10(p.values.calibrated.cate)~c(1:nrow(beta.mat.NONA.varlow)),xlab="Loci Position")
#points(-log10(p.values.calibrated.cate)[simu$causal]~simu$causal,
#       col="red",pch=16)
points(-log10(p.values.calibrated.cate)[which(-log10(p.values.calibrated.cate) > 4)]~which(-log10(p.values.calibrated.cate) > 4),
       col="blue",pch=21,cex=2)
abline(h = 4,col="red")




#### Taking global methylation means on non-trimed data vs trimmed data ####

#Single Read
# tot.count<-read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/counts_chr5_n50.txt",header = T,row.names = 1)
# m.count<-read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/mcounts_chr5_n50.txt",header=T,row.names = 1)
# beta<-m.count/tot.count
# beta.mat<-as.matrix(beta)
# beta.mat.NOz<-beta.mat
# beta.mat.NOz[is.nan(beta.mat.NOz)] = 0
# beta.mat.NOz[is.infinite(beta.mat.NOz)] = 1
# testing.mean.Est<-apply(beta.mat.NOz,1,mean)
# hist(testing.mean.Est)
# global.mean5<-mean(testing.mean.Est)
# length(testing.mean.Est)
# global.sd5<-sd(testing.mean.Est)

#global.methylation.baboon<-data.frame(chrom=c(1:5),mean=rep(NA,5),sd=rep(NA,5))
#global.methylation.baboon[2,2:3]<-c(global.mean2,global.sd2)


## Currently removes all CpGs with any NAs - this may need to be revised!@!
proportion.NA<-function(x){sum(is.na(x))/50}
testing.NA<-apply(beta.mat,1,proportion.NA)
hist(testing.NA)
beta.mat.NONA<-beta.mat[testing.NA == 0.00,]
dim(beta.mat.NONA)

baboon.file.names<-list.files(path = "DATA/Empirical_Data/baboon_example/BSSeq_Baboon/",pattern = "counts*")
chrom<-NULL
chrom.loc<-NULL
mean.chrom<-NULL
sd.chrom<-NULL
Hi<-0.95
Low<-0.05
for (i in 1:20){
  tot.count<-read.table(paste("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/",baboon.file.names[i],sep = ""),header = T,row.names = 1)
  m.count<-read.table(paste("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/",baboon.file.names[(20+i)],sep = ""),header=T,row.names = 1)
  beta<-m.count/tot.count
  beta.mat<-as.matrix(beta)
  proportion.NA<-function(x){sum(is.na(x))/50}
  testing.NA<-apply(beta.mat,1,proportion.NA)
  beta.mat.NONA<-beta.mat[testing.NA == 0.00,]
  
  ## Removes CpGs with high (>0.9) or low (<0.1) mean methylation
  testing.mean<-apply(beta.mat.NONA,1,mean)
  sd.c<-apply(beta.mat.NONA,1,function(x){sd(x,na.rm=TRUE)})
  rm<-(-which(testing.mean <= Low | testing.mean >= Hi ))
  reduce<-testing.mean[rm]
  lstring<-strsplit(names(reduce),"_")
  
  chrom<-c(chrom,sapply(lstring,function(x){paste(x[1])}))
  chrom.loc<-c(chrom.loc,sapply(lstring, function(x){paste(x[2])}))

  mean.chrom<-c(mean.chrom,reduce)
  sd.chrom<-c(sd.chrom,sd.c[rm])
}
chromGlobalMethyl<-data.frame(chrom = chrom,chrom.loc = chrom.loc, mean.chrom = mean.chrom, sd.chrom = sd.chrom)
length(chrom)
hist(mean.chrom)
plot(density(mean.chrom~chrom.loc))
head(chromGlobalMethyl)
### mean methylation is remarkable consistent in both mean and sd between all chromosomes
#EX
# chrom.label mean.chrom  sd.chrom
# 1   counts_chr1_n50.txt  0.6063702 0.2623622
# 2  counts_chr10_n50.txt  0.6146495 0.2561920
# 3  counts_chr11_n50.txt  0.6240210 0.2553796
# 4  counts_chr12_n50.txt  0.5993940 0.2691549



#### Taking single chrom1 data and inputting into EWAS_Generator Sim ####

source(file = "ewas_generator.R")

p.value<-1000
K.value<-7
sample(e.m$Predictor,size=p,replace = F)
dim(beta.mat.NONA.varlow)
mean.loci.beta<-apply(beta.mat.NONA.varlow,1,mean)

simu.rv <- ewas_generator(n = 50, 
                       p = p.value, 
                       K = K.value, 
                       freq=unname(sample(mean.loci.beta,size=p.value,replace = F)),
                       prop.variance = 0.1,
                       sd.U = runif(K.value),
                       mean.B = 5,
                       sd.B = 0.3,
                       sd.V = 0.1,
                       sigma = .1)
                       #setSeed = 5)

#### Simple run (lfmm) of revised ewas_sim ####
library(cate)
(K <- est.confounder.num(~ X, simu.rv,Y = simu.rv$Y, method = "ed"))
head(simu.rv$Y)

data <- gen.sim.data(n = 50, p = 100, r = 5)
X.data <- data.frame(X1 = data$X1)
dim(X.data)
dim(data$Y)
est.confounder.num(~ X1, X.data, data$Y, method = "ed")
head(data$Y)







