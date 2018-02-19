
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

## Useful function for finding rowxcol of nans or infinites
#which(is.infinite(beta.mat.NOz), arr.ind=TRUE)

## Currently removes all CpGs with any NAs - this may need to be revised!@!
proportion.NA<-function(x){sum(is.na(x))/50}
testing.NA<-apply(beta.mat,1,proportion.NA)
plot(density(testing.NA))
beta.mat[testing.NA == 0.00,]
dim(beta.mat.NONA)



## Removes CpGs with high (>0.9) or low (<0.1) mean methylation
testing.mean<-apply(beta.mat.NONA,1,mean)

hist(testing.mean)
beta.mat.NONA.var<-beta.mat.NONA[testing.mean > 0.1 & testing.mean < 0.9,]
dim(beta.mat.NONA.var)

## Remove CpGs with minimal variation (=<5% of overall dataset)
testing.sd<-apply(beta.mat.NONA.var,1,sd)
beta.mat.NONA.varlow<-beta.mat.NONA.var[testing.sd>quantile(testing.sd,probs = c(0.05)),]
dim(beta.mat.NONA.varlow)




## Estimate confounders ####

## Based on cate package
est.confounder.num(~ V1, p.variable, t(beta.mat.NONA.varlow), method = "ed")
est.confounder.num(~ V1, p.variable, t(beta.mat.NONA.varlow), method = "bcv")

### Estimating confounders based on sva package
num_sv     <- sva::num.sv(dat = beta.mat.NONA.varlow, mod = p.variable$V1, method = "be")
num_sv_l   <- sva::num.sv(dat = beta.mat.NONA.varlow, mod = p.variable$V1, method = "leek")


## lfmm ####
mod <- lfmm::lfmm_ridge(t(beta.mat.NONA.varlow),p.variable$V1,K = 5)
??lfmm_ridge()
## tests based on GLMs -- replaces lfmm_testing() function in lfmm package due to the binomial nature of data

p.value <- NULL
z.score <- NULL
for (j in 1:1000){
  mod.glm <- glm(t(beta.mat.NONA.varlow)[,j] ~ ., 
                 data = data.frame(p.variable$V1, mod$U),
                 binomial(link = "probit"))
  p.value[j] <- summary(mod.glm)$coeff[2,4]
  z.score[j] <- summary(mod.glm)$coeff[2,3] 
}



gif <- median(z.score^2)/0.456
p.values.calibrated <- pchisq(z.score^2/gif , df = 1, low = F)
hist(p.values.calibrated)
plot(-log10(p.values.calibrated)~c(1:1000),xlab="Loci Position")
#order(simu$causal)
#simu$causal[1]
#points(-log10(p.values.calibrated)[simu$causal]~simu$causal,
#       col="red",pch=16)
points(-log10(p.values.calibrated)[which(-log10(p.values.calibrated) > 3)]~which(-log10(p.values.calibrated) > 3),
       col="blue",pch=21,cex=2)
abline(h = 3,col="red")

## cate ####
output.EWAS.cate <- cate(~ p.variable$V1,p.variable, t(beta.mat.NONA.varlow), r = 5,adj.method = "rr")

## tests based on GLMs
p.value <- NULL
z.score <- NULL

for (j in 1:1000){
  mod.glm <- glm(t(beta.mat.NONA.varlow)[,j] ~ ., 
                 data = data.frame(p.variable$V1,output.EWAS.cate$Z),
                 binomial(link = "probit"))
  p.value.cate[j] <- summary(mod.glm)$coeff[2,4]
  z.score.cate[j] <- summary(mod.glm)$coeff[2,3] 
}

gif <- median(z.score.cate^2)/0.456 #genomic inflation factor
p.values.calibrated.cate <- pchisq(z.score.cate^2/gif , df = 1, low = F)
hist(p.values.calibrated.cate)
plot(-log10(p.values.calibrated.cate)~c(1:1000),xlab="Loci Position")
#points(-log10(p.values.calibrated.cate)[simu$causal]~simu$causal,
#       col="red",pch=16)
points(-log10(p.values.calibrated.cate)[which(-log10(p.values.calibrated.cate) > 4)]~which(-log10(p.values.calibrated.cate) > 4),
       col="blue",pch=21,cex=2)
abline(h = 4,col="red")



## Taking global methylation means on non-trimed data vs trimmed data

tot.count<-read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/counts_chr5_n50.txt",header = T,row.names = 1)
m.count<-read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/mcounts_chr5_n50.txt",header=T,row.names = 1)
beta<-m.count/tot.count
beta.mat<-as.matrix(beta)
beta.mat.NOz<-beta.mat
beta.mat.NOz[is.nan(beta.mat.NOz)] = 0
beta.mat.NOz[is.infinite(beta.mat.NOz)] = 1
testing.mean.Est<-apply(beta.mat.NOz,1,mean)
hist(testing.mean.Est)
global.mean5<-mean(testing.mean.Est)
global.sd5<-sd(testing.mean.Est)

#global.methylation.baboon<-data.frame(chrom=c(1:5),mean=rep(NA,5),sd=rep(NA,5))
global.methylation.baboon[2,2:3]<-c(global.mean2,global.sd2)


