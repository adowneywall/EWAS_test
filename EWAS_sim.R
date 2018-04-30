
#### libraries and scripts ####
library(ggplot2)
library(cate)
library(dplyr)
library(reshape2)
library(dSVA)
library(sva)
library(fdrtool)
library(RefFreeEWAS)
library(aod)
library(HRQoL)
library(lfmm)

source("EWAS_generator.R")
source("EWAS_methodTest.R")

#install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock'))

### Reading in process empirical data ####
q.cpg<-readRDS("DATA/Empirical_Data/quercus_example/cpg_pruned.RData")
q.chh<-readRDS("DATA/Empirical_Data/quercus_example/chh_pruned.RData")
q.chg<-readRDS("DATA/Empirical_Data/quercus_example/chg_pruned.RData")
b.cpg <- readRDS("DATA/Empirical_Data/baboon_example/baboon_prune.RData")

# used fitdistr function to estimate alpha and beta (shape 1 and 2) for a negative  binomial dist.
# this can be used to create 'reads' from simulated beta values to mimic real data
nb.baboon<-readRDS("DATA/Empirical_Data/baboon_example/negBinomFitParam.Rdata")

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
b.pca.sd<-b.pca$sdev[1:10] # possible include more if you want to include more pcs in sims

multi.sim.gen(n=c(100),
              p=c(2000),
              K=c(1,3),
              freq=freq.multi,
              prop.variance = 0.6, #seq(from=0.1,to=0.9,by=.2),
              sigma=0.4,
              mean.B = c(3), 
              sd.U=b.pca.sd,
              sd.Um=0.2,
              sd.Usd=0.1,
              sd.V=.5,
              rep=1,
              nb.t = nb.baboon,
              dir.name = "DATA/EWAS_Sims/", 
              sim.folder = "baboon.ReRun2")

sim <- readRDS("DATA/EWAS_Sims/baboon.ReRun2_2018-04-30/Sim1/Rep1/100_2000_0.6_3.Rdata")
dim(sim$Y)

screeplot(prcomp(sim$Y))

## transformations
simINVLogit<-invlink(sim$Y,type = "logit")
simQNORM <- qnorm(sim$Y)
sum(is.infinite(simQNORM))
#lfmm 
lfmm_out<-lfmm::lfmm_ridge(Y = simINVLogit,X = sim$X,K =1)
pval<-lfmm_test(Y = simINVLogit,X = sim$X,lfmm = lfmm_out,calibrate = "gif")

lfmm_out<-lfmm::lfmm_ridge(Y = simQNORM,X = sim$X,K =1)
pval<-lfmm_test(Y = simQNORM,X = sim$X,lfmm = lfmm_out,calibrate = "gif")
pval$gif

hist(pval$pvalue)
hist(pval$calibrated.pvalue)



## lm
lm_inv <- lm(simINVlogit~.)


## Examining the performance of difference methods on simulated baboon data
sim.folder<-"baboon.ReRun2_2018-04-30/"

sMeta<-simMeta(sim.folder) # enter the name of folder w/ directory, returns parameters list [param] and sim folder name

# Options for selecting a subset of your simulation list for testing
sub.list<-sMeta$param#[3,] #subset(sMeta$param,sMeta$param$n >= 100 & sMeta$param$p >= 2000)

#simRange<-seq(from=1,to=max(sim.list$X))
simRange<- sub.list$Sim

example<-simTest(sMeta$sim,simRange,sMeta$param,K.est = T,
                        lfmm.test = T,cate.test = T,RFEM.test = F,SVA.test = T,dSVA.test = T,
                        oracle.test = T,glm.test = T,rep.times = 1)

example<-simRead(SimRun = sMeta$sim,sim = sMeta$param,reps = 1)

lfmm_out<-lfmm::lfmm_ridge(Y = example$Y.list[[1]],X = example$X.list[[1]],K =10)
pca<-prcomp(example$Y.list[[1]])
screeplot(pca)
pval<-lfmm_test(Y = example$Y.list[[1]],X = example$X.list[[1]],lfmm = lfmm_out,calibrate = "gif")
pval$gif
hist(pval$pvalue)
hist(pval$calibrated.pvalue)


#dir.create("DATA/EWAS_Sims/quercus_Final_2018-04-19/Output")
#saveRDS(baboon.Final01,"DATA/EWAS_Sims/quercus_Final_2018-04-19/Output/FinalSimMethod.Rdata")
## currently K.est needs to = T, otherwise script won't run

#### Testing alternative regression methods ####
### Baboon Sample rerun with count data
# sim.folder<-"baboon.final2_2018-04-24/"
# 
# sMeta<-simMeta(sim.folder) # enter the name of folder w/ directory, returns parameters list [param] and sim folder name
# 
# # Options for selecting a subset of your simulation list for testing
# sub.list<-sMeta$param[1,] #subset(sMeta$param,sMeta$param$n >= 100 & sMeta$param$p >= 2000)
# 
# #simRange<-seq(from=1,to=max(sim.list$X))
# simRange<- sub.list$Sim
# sim<-sMeta$param[50,]
# reps=20
# sim.data<-NULL
# sim.data<-simRead(sim.folder,sim,reps)
# saveRDS(sim.data,"DATA/EWAS_Sims/baboonfinal2_Sim1Rep20.Rdata")

sim.data <- readRDS("SampleData/Output_baboonSim/baboonfinal2_Sim52Rep20.Rdata")\

sim.data <- readRDS("SampleData/Output_baboonSim/baboonfinal2_Sim52Rep20.Rdata")
bestim<-sim.data$mcount.list[[1]]/sim.data$tcount.list[[1]]
plot(sim.data$Y.list[[1]][,60]~bestim[,60])
summary(lm(sim.data$Y.list[[1]][,60]~bestim[,60]))

mc<-sim.data$mcount.list[[1]]
tc<-sim.data$tcount.list[[1]]
x<-sim.data$X.list[[1]] # phenotype
ll<-sim.data$B.list[[1]] # effect size of phenotype
y<-as.matrix(sim.data$Y.list[[1]])
lf<-sim.data$U.list[[1]]
as.matrix(unlist(sim.data$U.list[[1]],ncol=4))

summary(lm(x$V1~lf))
screeplot(prcomp(y),type="l")
covMat<-cor(lf[,1])
?cor()

p=2000
i=1
m=max(tc[,i])

fact<-cbind.data.frame(x=x$V1,lf)
lo_z<-NULL
lo_p<-NULL
lom_z<-NULL
lom_p<-NULL
qop_z<-NULL
qop_p<-NULL
qoc_z<-NULL
qoc_p<-NULL
bop_z<-NULL
bop_p<-NULL
boc_z<-NULL
boc_p<-NULL
bocm_z<-NULL
bocm_p<-NULL
pvalB1<-NULL
zscoreB1<-NULL
pvalB2<-NULL
zscoreB2<-NULL

for(i in 1:ncol(y)){
  print(i)
  print("1")
  
  linear_min<-summary(glm(y[,i]~x,
                          data=fact))
  lom_z[i]<-linear_min$coeff[2,3]
  lom_p[i]<-linear_min$coeff[2,4]
  b1<-summary(glm(y[,1]~.,
                    data=fact,
                 binomial(link="logit")))
  b2<-summary(glm(y[,1]~.,
                  data=fact,
                  binomial(link="probit")))
  anova(gaus,b1,b2)
  
  lom_z[i]<-gaus$coeff[2,3]
  lom_p[i]<-linear_min$coeff[2,4]
  
  print("2")
  binom_out_cmin<-summary(glm(cbind(mc[,i],tc[,i]-mc[,i])~x,
                           binomial(link="logit"),
                           data=fact))
  
  bocm_z[i]<-binom_out_cmin$coeff[2,3]
  bocm_p[i]<-binom_out_cmin$coeff[2,4]
}
  
  
  linear_out<-summary(glm(y[,i]~.,
                          data=fact))
  lo_z[i]<-linear_out$coeff[2,3]
  lo_p[i]<-linear_out$coeff[2,4]
  
  linear_min<-summary(glm(y[,i]~x,
                          data=fact))
  lom_z[i]<-linear_min$coeff[2,3]
  lom_p[i]<-linear_min$coeff[2,4]
  
  quasi_out_p<-summary(glm(y[,i]~.,
                           quasibinomial(link="logit"),
                           data=fact))
  qop_z[i]<-quasi_out_p$coeff[2,3]
  qop_p[i]<-quasi_out_p$coeff[2,4]
  
  quasi_out_c<-summary(glm(cbind(mc[,i],tc[,i]-mc[,i])~.,
                           quasibinomial(link="logit"),
                           data=fact))
  qoc_z[i]<-quasi_out_c$coeff[2,3]
  qoc_p[i]<-quasi_out_c$coeff[2,4]
  
  #fact<-cbind.data.frame(x$V1,lf)
  binom_out_p<-summary(glm(y[,i]~.,
                           binomial(link="logit"),
                           data=fact))
  
  bop_z[i]<-binom_out_p$coeff[2,3]
  bop_p[i]<-binom_out_p$coeff[2,4]
  
  binom_out_c<-summary(glm(cbind(mc[,i],tc[,i]-mc[,i])~.,
                           binomial(link="logit"),
                           data=fact))
  
  boc_z[i]<-binom_out_c$coeff[2,3]
  boc_p[i]<-binom_out_c$coeff[2,4]
  
  binom_out_cmin<-summary(glm(cbind(mc[,i],tc[,i]-mc[,i])~x,
                              binomial(link="logit"),
                              data=fact))
  
  bocm_z[i]<-binom_out_cmin$coeff[2,3]
  bocm_p[i]<-binom_out_cmin$coeff[2,4]
  
  bbinom_out_f<-betabin(cbind(mc[,i],tc[,i]-mc[,i])~.,~1,link = "logit",data=fact) # full beta binom
  bbinom_out_m<-betabin(cbind(mc[,i],tc[,i]-mc[,i])~x,~1,link = "logit",data=fact) # just x
  
  if(!is.na(bbinom_out_f@varparam)){
    zscoreB1[i]<-bbinom_out_f@fixed.param[2]/sqrt(bbinom_out_f@varparam[2,2])
    pvalB1[i]<-2 * (1 - pnorm(abs(zscoreB1)))
  }else{
    zscoreB1[i]<-999
    pvalB1[i]<-999
  }
  
  if(!is.na(bbinom_out_m@varparam)){
    zscoreB2[i]<-bbinom_out_m@fixed.param[2]/sqrt(bbinom_out_m@varparam[2,2])
    pvalB2[i]<-2 * (1 - pnorm(abs(zscoreB2)))
  }else{
    zscoreB2[i]<-999
    pvalB2[i]<-999
  }
  print(i)
}
score.df<-data.frame(lo_z,qop_z,qoc_z,bop_z,boc_z,zscoreB1,zscoreB2)
scoreFull.df<-data.frame(pos=c(1:2000),lo_z,qop_z,qoc_z,bop_z,boc_z,zscoreB1,zscoreB2,
                         lom_z,bocm_z,x,ll)
scoreFull.df$CL<-0
scoreFull.df$CL[CL]<-1
zscoreB1[zscoreB1 > 100]
plot(abs(lo_z)~lo_p)
plot(abs(lo_z)~abs(qop_z))
plot(abs(lo_z)~abs(qoc_z))
plot(abs(qop_z)~abs(qoc_z))
plot(abs(boc_z)~abs(qoc_z))
plot(abs(boc_z)~abs(bop_z))
plot(abs(boc_z)~abs(zscoreB2))
plot(abs(boc_z)~abs(zscoreB1),xlim=c(0,60))
plot(abs(zscoreB2)~abs(zscoreB1),xlim=c(0,60))

ggplot(scoreFull.df,aes(x=abs(zscoreB1),y=abs(zscoreB2),colour=as.factor(CL))) + geom_point() +
  xlim(0,75)
ggplot(scoreFull.df,aes(x=abs(lo_z),y=abs(lom_z),colour=as.factor(CL))) + geom_point()
ggplot(scoreFull.df,aes(x=abs(boc_z),y=abs(bocm_z),colour=as.factor(CL))) + geom_point()


ggplot(scoreFull.df,aes(x=abs(lo_z),y=abs(ll),colour=as.factor(CL))) + geom_point()
ggplot(scoreFull.df,aes(x=abs(boc_z),y=abs(ll),colour=as.factor(CL))) + geom_point()
ggplot(scoreFull.df,aes(x=abs(zscoreB1),y=abs(ll),colour=as.factor(CL))) + geom_point() +
  xlim(0,75)



library(reshape2)
full<-melt(score.df) # long data version
full$pos<-rep(c(1:2000),7) # add relative positions
full<-subset(full,full$value < 250)

ggplot(scoreFull.df,aes(x=pos,y=abs(zscoreB2))) + geom_point(alpha=0.2) +
  geom_point(aes(x=CL,y=scoreFull.df$zscoreB2[CL],alpha=0.2))
  geom_vline(xintercept = c(CL),alpha=0.2)

### Full melted plot
ggplot(full,aes(x=pos,y=abs(value),group=variable,colour=variable)) + 
  geom_point(aes(alpha=0.2)) +
  geom_vline(xintercept = c(CL),alpha=0.2)

#### Example BBmm ####
?BBmm()
# Defining the parameters
k <- 100
m <- 10
phi <- 0.5
beta <- c(1.5,-1.1)
sigma <- 0.5

# Simulating the covariate and random effects
x <- runif(k,0,10)
X <- model.matrix(~x)
z <- as.factor(rBI(k,4,0.5,2))
Z <- model.matrix(~z-1)
u <- rnorm(5,0,sigma)


# The linear predictor and simulated response variable
eta <- beta[1]+beta[2]*x+crossprod(t(Z),u)
p <- 1/(1+exp(-eta))
y <- rBB(k,m,p,phi)
dat <- data.frame(cbind(y,x,z))
dat$z <- as.factor(dat$z)

# Apply the model
model <- BBmm(fixed.formula = y~x,random.formula = ~z,m=c(y+m),data=dat)
model


#### Short work up of quercus example ####
#Sim info:
# n=c(50,250,500),
# p=c(2000),
# K=c(1,3),
# freq=freq.multi,
# prop.variance = seq(from=0.1,to=0.9,by=.1),
# sigma=0.2,
# mean.B = c(1,2,5), 
# sd.U=b.pca.sd,
# sd.Um=1,
# sd.Usd=0.2,
# rep=20,
# nb.t = NULL,
# dir.name = "DATA/EWAS_Sims/", 
# sim.folder = "quercus_Final"
# only sims sMeta$param[1:108,] 1-108

working.q.df<-readRDS("DATA/EWAS_Sims/quercus_Final_2018-04-19/Output/FinalSimMethod.Rdata")$full.df

sub<-subset(working.q.df,working.q.df$n == 500)
sub
q.df<-readRDS("DATA/EWAS_Sims/quercus_Final_2018-04-19/Output/FinalSimMethod.Rdata")
q.df$full.df
q
values<-q.df[[2]][[3]][[1]]
new.val<-values[[1]]
new.val[,1]
hist(new.val[,1])
colnames(values)
ss<-working.q.df
ss.sub <- subset(ss,ss$n == 250)
ss.sub <- subset(ss.sub,ss.sub$p == 2000)
library(ggplot2)
ggplot(ss.sub,aes(x=method,y=as.numeric(et0))) + geom_boxplot()
ggplot(ss.sub,aes(x=method,y=as.numeric(power_fdr))) + geom_boxplot()
ggplot(ss.sub,aes(x=method,y=as.numeric(error_fdr))) + geom_boxplot()

sub.lfmm <- subset(ss,ss$method == "lfmm")
sub.lfmm <- subset(sub.lfmm,sub.lfmm$p == 2000)
sub.lfmm <- subset(sub.lfmm,sub.lfmm$n == 250)
ggplot(sub.lfmm,aes(x=as.factor(prop.variance),y=as.numeric(et0))) + geom_boxplot()
ggplot(sub.lfmm,aes(x=as.factor(prop.variance),y=as.numeric(power_fdr))) + geom_boxplot()

sub.cate <- subset(ss,ss$method == "cate")
sub.cate <- subset(sub.cate,sub.cate$p == 2000)
sub.cate <- subset(sub.cate,sub.cate$n == 50)
ggplot(sub.cate,aes(x=as.factor(prop.variance),y=as.numeric(et0))) + geom_boxplot()
ggplot(sub.cate,aes(x=as.factor(prop.variance),y=as.numeric(power_fdr))) + geom_boxplot()

sub.oracle <- subset(ss,ss$method == "oracle")
sub.oracle <- subset(sub.oracle,sub.oracle$p == 2000)
sub.oracle <- subset(sub.oracle,sub.oracle$n == 250)
ggplot(sub.oracle,aes(x=as.factor(prop.variance),y=as.numeric(et0))) + geom_boxplot()
ggplot(sub.oracle,aes(x=as.factor(prop.variance),y=as.numeric(power_fdr))) + geom_boxplot()

# basic summary of eta0,power, and error
ss.summary<- ss %>% group_by(n,p,prop.variance,method,mean.B,K) %>%
  summarise(mean.pFDR = mean(as.numeric(power_fdr),na.rm=T),mean.et0 = mean(as.numeric(et0),na.rm=T),
            mean.eFDR = mean(as.numeric(error_fdr),na.rm=T))
#ss.summary$mean.B<-as.factor(ss.summary$mean.B)
ss.summary<-transform(ss.summary,mean.B = factor(mean.B,levels=1:2,labels=c("Effect Size 1","Effect Size 2")))
ss.summary$Kfact<-as.factor(ss.summary$K)
levels(ss.summary$mean.B)

## Summary of power that has been fdr corrected
ggplot(ss.summary,aes(x=prop.variance,y=mean.pFDR,group=interaction(as.factor(n),Kfact),colour=as.factor(n),linetype=Kfact)) + 
  geom_line() + geom_point() +facet_grid(mean.B~as.factor(method)) +
  theme_bw() + labs(x="Proportion of confounding",y="Estimated Power (FDR <= 0.05)",colour="Sample Size")
ggsave(filename = "DATA/EWAS_Sims/quercus_Final_2018-04-19/Output/q_Final01_FDR.png")
## Summary of eto0
ggplot(ss.summary,aes(x=prop.variance,y=mean.et0,colour=as.factor(n))) + 
  geom_line() + geom_point() +facet_grid(mean.B~as.factor(method)) +
  theme_bw() + labs(x="Proportion of confounding",y="eta0",colour="Sample Size")
ggsave(filename = "DATA/EWAS_Sims/quercus_Final_2018-04-19/Output/q_Final01_et0.png")
## Summary of error
ggplot(ss.summary,aes(x=prop.variance,y=mean.eFDR,colour=as.factor(n))) + 
  geom_line() + geom_point() +facet_grid(mean.B~as.factor(method)) +
  theme_bw() + labs(x="Proportion of confounding",y="Estimated Error (FDR <= 0.05)",colour="Sample Size")

### For looking at hitRank as a measure of quality
hitRankDF<-data.frame(method=NA,sim=NA,rep=NA,q=NA,hr=NA,tr=NA)
folder<-"DATA/EWAS_Sims/quercus_Final_2018-04-19/"
num.HIT<-150 # only look at top 150 values

for(i in 1:108){ ## number of sims in simulation run object
  for(j in 1:20){ ## number of reps 
    cs<-read.csv(Sys.glob(paste(folder,"Sim",i,"/Rep",j,"/caus*.csv",sep="")))[,2]
    for(k in c(2,3,4,5,7,8)){ # this refers to the position of each method in list object
      if(k == 8){
        temp<-matrix(unlist(baboon.Final01[[k]][[i]][[j]]),ncol=5)
        ct<-cbind.data.frame(loci=c(1:nrow(temp)),temp)
        names(ct)<-c("loci","p","z","pcal","fdr.p","fdr.q")
        
      }else{
        temp<-matrix(unlist(baboon.Final01[[k]][[i]][[j]]),ncol=6)
        ct<-cbind.data.frame(loci=c(1:nrow(temp)),temp)
        names(ct)<-c("loci","p","z","pcal","gif","fdr.p","fdr.q")
      }
      ct$logp<-(-log(ct$pcal))
      ct$cs<-0
      ct$cs[cs]<-1
      ct<-ct[rev(order(ct$logp)),] # ranking based on log of the calibrate p value
      # Examined the fdrtools pval, will calculat the same number as the gif
      ct$hitRank<-c(1:nrow(ct))
      ct$sumRank<-cumsum(ct$cs)
      mSim<-cbind.data.frame(method=rep(k,num.HIT),sim=rep(i,num.HIT),rep=rep(j,num.HIT),
                             q=ct$fdr.q[1:num.HIT],hr=ct$hitRank[1:num.HIT],tr=ct$sumRank[1:num.HIT])
      
      hitRankDF<-rbind.data.frame(hitRankDF,mSim) 
    }
  }
  print(i)
}

#ggplot(mSim,aes(x=hr,y=tr)) + geom_line()

hitRankDFfix<-hitRankDF[-1,]
saveRDS(object = hitRankDFfix,file = "DATA/EWAS_Sims/quercus_Final_2018-04-19/Output/hitRankOutput.Rdata")
hitRankDFfix<-readRDS("DATA/EWAS_Sims/quercus_Final_2018-04-19/Output/hitRankOutput.Rdata")
hitRankSum <- hitRankDFfix %>% group_by(method,sim,hr) %>%
  summarise(mean.q=mean(q),mean.tr=mean(tr),sd.tr=sd(tr))
library(ggplot2)
singleSim<-subset(ss,ss$method == "lfmm" & ss$rep == 1)
singleSim<-singleSim[,c(1,2,4,7,10)]
names(singleSim)<-c("sim",names(singleSim)[2:5])
hitRankupdate<- inner_join(x=hitRankSum,y=singleSim)
hitRankupdate<-transform(hitRankupdate,method = factor(method,levels=c(2,3,4,5,7,8),labels=c("lfmm","cate","sva","dsva","oracle","glm")))
hitRankupdate<-transform(hitRankupdate,mean.B = factor(mean.B,levels=1:2,labels=c("Mean effect Size 1","Mean effect Size 2")))
hitRankupdate$kfactor<-as.factor(hitRankupdate$K)

hitRanksingleFactor<-subset(hitRankupdate,hitRankupdate$K == "1")
ggplot(hitRanksingleFactor,aes(x=hr,y=mean.tr,interaction(as.factor(n),as.factor(prop.variance)),colour=as.factor(prop.variance),linetype=as.factor(n))) +
  geom_line() + facet_grid(as.factor(mean.B)~method) + xlim(0,150) + ylim(0,55) + geom_hline(aes(yintercept=50,colour="red")) +
  theme_bw() + labs(x = "Hit Rank",y = "Mean Causal Hits (30 simulations)")


hitRanksingleEffectSize<-subset(hitRankupdate,hitRankupdate$mean.B == "Mean effect Size 1")
ggplot(hitRanksingleEffectSize,aes(x=hr,y=mean.tr,interaction(as.factor(n),as.factor(prop.variance)),colour=as.factor(prop.variance),linetype=as.factor(n))) +
  geom_line() + facet_grid(K~method) + xlim(0,150) + ylim(0,55) + geom_hline(aes(yintercept=50,colour="red")) +
  theme_bw() + labs(x = "Hit Rank",y = "Mean Causal Hits (30 simulations)")


#### Short work up for Baboon based data work up ####
# Based on simulations 
# with variable:
# n = 50,100,250,500
# proprotion of confounding = .1-.9 in .1 intervals
# k = 1 (based on pca/scree plot of real baboon data)
# mean.b (effect size) - 1 or 2
# sigma=0.2,
# sd.U=b.pca.sd,
# sd.Um=1,
# sd.Usd=0.2,
# sim data has 200reps (since deleted), I only ran the test on 20 (10 left)
# only 72 sims
working.df<-readRDS("DATA/EWAS_Sims/baboon_Final_2018-04-16/Output/Final01.Rdata")$full.df

ss<-working.df
ss.sub <- subset(ss,ss$n == 250)
ss.sub <- subset(ss.sub,ss.sub$p == 2000)
library(ggplot2)
ggplot(ss.sub,aes(x=method,y=as.numeric(et0))) + geom_boxplot()
ggplot(ss.sub,aes(x=method,y=as.numeric(power_fdr))) + geom_boxplot()
ggplot(ss.sub,aes(x=method,y=as.numeric(error_fdr))) + geom_boxplot()

sub.lfmm <- subset(ss,ss$method == "lfmm")
sub.lfmm <- subset(sub.lfmm,sub.lfmm$p == 2000)
sub.lfmm <- subset(sub.lfmm,sub.lfmm$n == 250)
ggplot(sub.lfmm,aes(x=as.factor(prop.variance),y=as.numeric(et0))) + geom_boxplot()
ggplot(sub.lfmm,aes(x=as.factor(prop.variance),y=as.numeric(power_fdr))) + geom_boxplot()

sub.cate <- subset(ss,ss$method == "cate")
sub.cate <- subset(sub.cate,sub.cate$p == 2000)
sub.cate <- subset(sub.cate,sub.cate$n == 50)
ggplot(sub.cate,aes(x=as.factor(prop.variance),y=as.numeric(et0))) + geom_boxplot()
ggplot(sub.cate,aes(x=as.factor(prop.variance),y=as.numeric(power_fdr))) + geom_boxplot()

sub.oracle <- subset(ss,ss$method == "oracle")
sub.oracle <- subset(sub.oracle,sub.oracle$p == 2000)
sub.oracle <- subset(sub.oracle,sub.oracle$n == 250)
ggplot(sub.oracle,aes(x=as.factor(prop.variance),y=as.numeric(et0))) + geom_boxplot()
ggplot(sub.oracle,aes(x=as.factor(prop.variance),y=as.numeric(power_fdr))) + geom_boxplot()

# basic summary of eta0,power, and error
ss.summary<- ss %>% group_by(n,p,prop.variance,method,mean.B) %>%
  summarise(mean.pFDR = mean(as.numeric(power_fdr),na.rm=T),mean.et0 = mean(as.numeric(et0),na.rm=T),
            mean.eFDR = mean(as.numeric(error_fdr),na.rm=T))
#ss.summary$mean.B<-as.factor(ss.summary$mean.B)
ss.summary<-transform(ss.summary,mean.B = factor(mean.B,levels=1:2,labels=c("Effect Size 1","Effect Size 2")))
levels(ss.summary$mean.B)

## Summary of power that has been fdr corrected
ggplot(ss.summary,aes(x=prop.variance,y=mean.pFDR,colour=as.factor(n))) + 
  geom_line() + geom_point() +facet_grid(mean.B~as.factor(method)) +
  theme_bw() + labs(x="Proportion of confounding",y="Estimated Power (FDR <= 0.05)",colour="Sample Size")
ggsave(filename = "DATA/EWAS_Sims/baboon_Final_2018-04-16/Output/powerFDRcorrMultiMethod.pdf")
## Summary of eto0
ggplot(ss.summary,aes(x=prop.variance,y=mean.et0,colour=as.factor(n))) + 
  geom_line() + geom_point() +facet_grid(mean.B~as.factor(method)) +
  theme_bw() + labs(x="Proportion of confounding",y="eta0",colour="Sample Size")
ggsave(filename = "DATA/EWAS_Sims/baboon_Final_2018-04-16/Output/et0MultiMethod.pdf")
## Summary of error
ggplot(ss.summary,aes(x=prop.variance,y=mean.eFDR,colour=as.factor(n))) + 
  geom_line() + geom_point() +facet_grid(mean.B~as.factor(method)) +
  theme_bw() + labs(x="Proportion of confounding",y="Estimated Error (FDR <= 0.05)",colour="Sample Size")


### For looking at hitRank as a measure of quality
hitRankDF<-data.frame(method=NA,sim=NA,rep=NA,q=NA,hr=NA,tr=NA)
folder<-"DATA/EWAS_Sims/baboon_Final_2018-04-16/"
num.HIT<-150 # only look at top 150 values

for(i in 1:72){ ## number of sims in simulation run object
  for(j in 1:20){ ## number of reps 
   cs<-read.csv(Sys.glob(paste(folder,"Sim",i,"/Rep",j,"/caus*.csv",sep="")))[,2]
   for(k in c(2,3,4,5,7,8)){ # this refers to the position of each method in list object
     if(k == 8){
       temp<-matrix(unlist(baboon.Final01[[k]][[i]][[j]]),ncol=5)
       ct<-cbind.data.frame(loci=c(1:nrow(temp)),temp)
       names(ct)<-c("loci","p","z","pcal","fdr.p","fdr.q")
       
     }else{
       temp<-matrix(unlist(baboon.Final01[[k]][[i]][[j]]),ncol=6)
       ct<-cbind.data.frame(loci=c(1:nrow(temp)),temp)
       names(ct)<-c("loci","p","z","pcal","gif","fdr.p","fdr.q")
     }
     ct$logp<-(-log(ct$pcal))
     ct$cs<-0
     ct$cs[cs]<-1
     ct<-ct[rev(order(ct$logp)),] # ranking based on log of the calibrate p value
     # Examined the fdrtools pval, will calculat the same number as the gif
     ct$hitRank<-c(1:nrow(ct))
     ct$sumRank<-cumsum(ct$cs)
     mSim<-cbind.data.frame(method=rep(k,num.HIT),sim=rep(i,num.HIT),rep=rep(j,num.HIT),
                            q=ct$fdr.q[1:num.HIT],hr=ct$hitRank[1:num.HIT],tr=ct$sumRank[1:num.HIT])
     
     hitRankDF<-rbind.data.frame(hitRankDF,mSim) 
     }
  }
  print(i)
}

#ggplot(mSim,aes(x=hr,y=tr)) + geom_line()

hitRankDFfix<-hitRankDF[-1,]
saveRDS(object = hitRankDFfix,file = "DATA/EWAS_Sims/baboon_Final_2018-04-16/Output/HitRankEst.Rdata")
hitRankDFfix<-readRDS("DATA/EWAS_Sims/baboon_Final_2018-04-16/Output/HitRankEst.Rdata")
hitRankSum <- hitRankDFfix %>% group_by(method,sim,hr) %>%
  summarise(mean.q=mean(q),mean.tr=mean(tr),sd.tr=sd(tr))

singleSim<-subset(ss,ss$method == "lfmm" & ss$rep == 1)
singleSim<-singleSim[,c(1,2,7,10)]
names(singleSim)<-c("sim",names(singleSim)[2:4])
hitRankupdate<- inner_join(x=hitRankSum,y=singleSim)
hitRankupdate<-transform(hitRankupdate,method = factor(method,levels=c(2,3,4,5,7,8),labels=c("lfmm","cate","sva","dsva","oracle","glm")))
hitRankupdate<-transform(hitRankupdate,mean.B = factor(mean.B,levels=1:2,labels=c("Mean effect Size 1","Mean effect Size 2")))
ggplot(hitRankupdate,aes(x=hr,y=mean.tr,interaction(as.factor(n),as.factor(prop.variance)),colour=as.factor(prop.variance),linetype=as.factor(n))) +
  geom_line() + facet_grid(as.factor(mean.B)~method) + xlim(0,150) + ylim(0,55) + geom_hline(aes(yintercept=50,colour="red")) +
  theme_bw() + labs(x = "Hit Rank",y = "Mean Causal Hits (30 simulations)")
ggsave(filename = "DATA/EWAS_Sims/baboon_Final_2018-04-16/Output/HitPlotMultiMethod.pdf")

hitRankDFfixlabel<-transform(hitRankDFfix,method = factor(method,levels=c(2,3,4,5,7,8),labels=c("lfmm","cate","sva","dsva","oracle","glm")))
hitRankDFfix.sub1<-subset(hitRankDFfixlabel,hitRankDFfixlabel$method=="cate")
hitRankDFfix.sub1<- inner_join(x=hitRankDFfix.sub1,y=singleSim)
hitRankDFfix.sub2<-subset(hitRankDFfix.sub1,hitRankDFfix.sub1$mean.B==2)
ggplot(hitRankDFfix.sub2,aes(x=hr,y=tr,group=as.factor(rep),colour=as.factor(rep))) +
  geom_line() + facet_grid(.~as.factor(prop.variance)) + xlim(0,150)
ggsave(filename = "DATA/EWAS_Sims/baboon_Final_2018-04-16/Output/HitPlotCatemultiePropVar.pdf")

ggplot(hitRankDFfix.sub2, aes(x=hr,y=tr, colour=as.factor(prop.variance))) +
  stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill=as.factor(prop.variance)), alpha=0.3) +
  theme_bw()
ggsave(filename = "DATA/EWAS_Sims/baboon_Final_2018-04-16/Output/HitPlotCateSpline.png")

hit.sub2.sum<-hitRankDFfix.sub2 %>% group_by(prop.variance,hr) %>%
  summarize(mean=mean(tr),sd=sd(tr))
ggplot(hit.sub2.sum,aes(x=hr,y=mean,group=as.factor(prop.variance),colour=as.factor(prop.variance))) +
  geom_line() +
  theme_bw() + xlim(0,150) + geom_hline(yintercept=50)
ggsave(filename = "DATA/EWAS_Sims/baboon_Final_2018-04-16/Output/HitPlotCateSum.pdf")

wmean(ss$K.bcv,na.rm=T)
mean(ss$K.leek)

### Examining performance of logit-transformed beta-values. ####
P<-2000
N<-50
log.test<-ewas_generator(n = N,p = P,K = 1,freq = sample(q.cpg.means$y,size = P,replace=T),
                         prop.causal = 0.025,prop.variance = 0.5,
                         sd.U = 1,sd.V = 1,sd.B = .5,mean.B = 2,sigma = .2)
hist(log.test$Y[,30],xlim=c(0,1))
hist(log.test$Y.logit[,30],xlim=c(0,1))


