
TOBECOMPLETED

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

mean.p  = 0 # mean value
var.p   = 1 # variance
coef.p  = 0.3 # coefficient
inter.p = 0.5 # intercept

### Latent factors / additional covariates
 # Five latent factors.
 # Normal distributions: N(0,5),
 # N(3,1), N(0,1), N(2,4), N(0,3)

mean.lats <- c(0,3,0,2,0)
var.lats <- c(5,1,1,4,3)

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
 
 
 kaushalsim1 <- function(times = times,
                         n = n,
                         p = p,
                         outliers = outliers,
                         mean.p  = mean.p, 
                         var.p   = var.p, 
                         coef.p  = coef.p,
                         inter.p = inter.p,
                         mean.lats = mean.lats,
                         var.lats = var.lats,
                         
 ){
   
 }