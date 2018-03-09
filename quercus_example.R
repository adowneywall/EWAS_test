
library(readr)
library(data.table)
library(dplyr)
q.meth<-read_table2("DATA/Empirical_Data/quercus_example/set3.meth.full.tab.gz",col_names = T)

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
#rm(q.meth.s1,q.meth.s2,q.meth.s3)
#rm(q_meth.start,q_meth.add1,q_meth.add2)

q.cg.mat<-as.matrix(q.cg[,6:63])
q.cg.means<-rowMeans(q.cg.mat,na.rm=TRUE)
hist(q.cg.means)

#Basic Diagnostics

# Taking probe means
q.cg.means<-rowMeans(q.cg.mat,na.rm=TRUE)
mean(q.cg.means)
hist(abs(q.cg.means-0.5))
hist(q.cg.means)
#Taking probe SD
q.cg.sd<-apply(q.cg.mat,1,function(x){sd(x)})
hist(q.cg.sd)
paste(dim(q.cg.mat)[1], "CpG probes") # 261652 probes left

#Removed NAs (at threshold of 10)
q.cg.NAs<-apply(q.cg.mat,1,function(x){sum(is.na(x))})
hist(q.cg.NAs)
NA.thres<-10
q.cg.NAomit<-q.cg.mat[which(q.cg.NAs <= NA.thres),]
hold<-apply(q.cg.NAomit,1,function(x){sum(is.na(x))})
hist(hold)
q.cg.means.lowNA<-rowMeans(q.cg.NAomit,na.rm=TRUE)
hist(q.cg.means.lowNA)
mean(q.cg.means.lowNA)
paste(dim(q.cg.NAomit)[1]," probes remaing after NA pruning") # 73224 probes left
paste((1-dim(q.cg.NAomit)[1]/dim(q.cg.mat)[1])*100, "% of probes removed with NA threshold of", NA.thres)

## Removing SNPs with little variation
# removes means less than 0.01 and greater than 0.09
UVT<-0.99
LVT<-0.01 
q.cg.lowvarOmit<-q.cg.NAomit[which(q.cg.means.lowNA >= LVT & q.cg.means.lowNA <= UVT),]
# hist of mean prob vals after very high and low var removed
adj.quercus.means<-q.cg.means.lowNA[which(q.cg.means.lowNA >= LVT & q.cg.means.lowNA <= UVT)]
hist(adj.quercus.means)
hist(sample(x =adj.quercus.means,size =1000 ,replace=T))
# mean of mean prob value after v. high and low var removed
mean(q.cg.means.lowNA[which(q.cg.means.lowNA >= LVT & q.cg.means.lowNA <= UVT)])
paste(dim(q.cg.lowvarOmit)[1],"probes remaining after probes with low variation removed") #12980 probes remaining
hist(q.cg.lowvarOmit)

# removes sd with greater than var of 0.01 (not currently doing this)
q.cg.sd<-apply(q.cg.lowvarOmit,1,function(x){sd(x,na.rm=T)})
#q.cg.lowvarOmit<-q.cg.lowvarOmit[which(q.cg.sd > 0.01),]
hist(apply(q.cg.lowvarOmit,1,function(x){sd(x,na.rm=T)}))
plot(density(q.cg.lowvarOmit[1,],na.rm=T))
N<-2
M<-100
for(i in N:M){
  #plot(density(q.cg.lowvarOmit[i,],na.rm=T),main=paste(i))
  hist(q.cg.lowvarOmit[i,],na.rm=T,main=paste(i))
}
hist(q.cg.lowvarOmit[74,])
## Suggesting probe beta-values are not normal...
shapiro.test(q.cg.lowvarOmit[8,])
library(nortest)
pearson.test(q.cg.lowvarOmit[8,])
ks.test(q.cg.lowvarOmit[8,])


### fitting different distributions to data

x <- rbeta(100,0.1,10)
hist(x,breaks=100)
## now do this directly with more control.
store<-fitdistr(x,densfun = "beta", list(shape1 = 1, shape2 = 10), lower = 0.001)
store2<-fitdist(x,"beta")
store2
gofstat(store2)
new.dist<-rbeta(1000,0.1101489,7.0510089)
hist(new.dist,breaks=100)

mean(q.cg.lowvarOmit[,20],na.rm=T)
# Creat data fram of probe means and sds
quercus.var.df<-data.frame(mean=adj.quercus.means,sd=q.cg.sd)
library(fitdistrplus)
gofstat(q.cg.lowvarOmit[20,])
hold<-rbeta(10000,10,1)
mean(hold)
(10+1)/10
(10*1)/((10+1)^2*(10+1+1))
hist(rbeta(10000,10,1),breaks=1000)

min(rbeta(1000,10,1))
?rbeta()
## Variance distribution across mean bins w/ quercus dataset
bins<-seq(from=0,to=1,by=(1/20))
sd.means<-NULL
for(i in 1: (length(bins)-1)){
  print(i)
  temp.b<-which(quercus.var.df$mean>bins[i] & quercus.var.df$mean<=bins[i+1])
  sd.means<-c(sd.means,mean(quercus.var.df$sd[temp.b]))
}
plot(sd.means~c(as.factor(paste(bins[1:20],"-",bins[2:21]))))
hist(quercus.var.df$mean)
hist(quercus.var.df$sd)



## Variance distributino across mean bins w/ simulation data
sim1<-ewas_generator(50,2000,5,
                     freq = sample(x = adj.quercus.means,size = 2000 ,replace=T))
sim1<-ewas_generator(50,1000,5,mean.B=1,freq = NULL)
bins<-seq(from=0,to=1,by=(1/20))
sd.means<-NULL
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

