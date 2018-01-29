


#### Scenario 1 ####
 # We simulated DNA methylation data at 2000 CpG sites across 
 # 600 samples, of which the first n CpG sites were associated 
 # with covariates of interest (e.g., level of arsenic exposure) 
 # and a set of latent variables, and the remaining CpG sites were
 # only associated with the latent variables. The set of latent 
 # variables represent “cell types”. 

n = 600 # samples
p = 2000 # number of loci

### Phenotype / covariate of interest

# mean.p = 0
# var.p = 1
# coef.p = 0.3
# inter.p = 0.5

### Latent factors / additional covariates
 # Five latent factors.
 # Normal distributions: N(0,5),
 # N(3,1), N(0,1), N(2,4), N(0,3)

mean.1.lat = 0
var.1.lat = 5
mean.2.lat = 3
var.2.lat = 1
mean.3.lat = 0
var.3.lat = 1
mean.4.lat = 2
var.4.lat = 4
mean.5.lat = 0
var.5.lat = 3

 # The association of DNA methylation and the latent
 # variables was assumed linear and the coefficients
 # were generated from N (0.5, 0.01). 

 assoc.mat.lat = matrix(nrows=2,ncol=5)
 for(i in 1:5){
   FINISH
 }
 
 # The distribution of random errors in the linear 
 # regressions was assumed to be Normal with mean 0 
 # and variance 1.2 for the n CpGs, mean 0 and variance 
 # of 1.2 for the next 100 CpGs, and mean 0 and variance 2
 # for the remaining CpGs. The last setting with larger variance 
 # in random errors was for situations where the influence of 
 # cell types on DNA methylation was weaker.
 
 
 kaushalsim1 <- function(
 ){
   
 }