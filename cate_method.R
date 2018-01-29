#### Example of model 'cate' ####
 # Linear mixed effect model

## Libraries and Dependencies
library(cate)
source(file = "ewas_generator.R")
?cate()

#### Simple Example #### 
## simulate a dataset with 100 observations, 1000 variables and 5 confounders
data <- gen.sim.data(n = 100, p = 1000, r = 5,seed = 5)

X.data <- data.frame(X1 = data$X1) # Intercept

## linear regression without any adjustment
output.naive <- cate(~ X1, X.data, Y = data$Y, r = 0, adj.method = "naive")
## confounder adjusted linear regression
output <- cate(~ X1, X.data, Y = data$Y, r = 5)

## plot the histograms of unadjusted and adjusted regression statistics
par(mfrow = c(1, 2))
hist(output.naive$beta.t)
hist(output$beta.t)


est.confounder.num(~ X1, X.data, data$Y, method = "ed")
est.confounder.num(~ X1, X.data, data$Y, method = "bcv")

#### Example using the generated EWAS data ####
simu <- ewas_generator(n = 200, 
                       p = 1000, 
                       K = 7, 
                       prop.variance = 0.5,
                       sd.U = c(.5, .25, .1,.3,.2,.4,.3),
                       mean.B = .2,
                       sd.B = 0.1,
                       sd.V = 0.5,
                       sigma = .1,
                       setSeed = 5) 


### Estimate the number of latent factors
### estimating latent factors
est.confounder.num(~ X, simu, simu$Y, method = "ed")
est.confounder.num(~ X, simu, simu$Y, method = "bcv",nRepeat = 10)

### Estimating confounders based on sva package
num_sv     <- sva::num.sv(dat = t(simu$Y), mod = simu$X, method = "be")
num_sv_l   <- sva::num.sv(dat = t(simu$Y), mod = simu$X, method = "leek")

#output.EWAS <- cate(~ simu$X,simu, simu$Y, r = 0, adj.method = "naive")
output.EWAS.cate <- cate(~ simu$X,simu, simu$Y, r = 20)
output.EWAS.cate <- cate(~ simu$X,simu, simu$Y, r = 6,adj.method = "rr")

## tests based on GLMs
p.value.cate <- NULL
z.score.cate <- NULL
for (j in 1:1000){
  mod.glm <- glm(simu$Y[,j] ~ ., 
                 data = data.frame(simu$X, output.EWAS.cate$Z),
                 binomial(link = "probit"))
  p.value.cate[j] <- summary(mod.glm)$coeff[2,4]
  z.score.cate[j] <- summary(mod.glm)$coeff[2,3] 
}

gif <- median(z.score.cate^2)/0.456 #genomic inflation factor
p.values.calibrated.cate <- pchisq(z.score.cate^2/gif , df = 1, low = F)
hist(p.values.calibrated.cate)
plot(-log10(p.values.calibrated.cate)~c(1:1000),xlab="Loci Position")
points(-log10(p.values.calibrated.cate)[simu$causal]~simu$causal,
       col="red",pch=16)
points(-log10(p.values.calibrated.cate)[which(-log10(p.values.calibrated.cate) > 4)]~which(-log10(p.values.calibrated.cate) > 4),
       col="blue",pch=21,cex=2)
abline(h = 4,col="red")

# power
length(which(-log10(p.values.calibrated.cate) > 4) %in% simu$causal)
# Type one  error
length(simu$causal %in% which(-log10(p.values.calibrated.cate) > 4))
length(simu$causal)
length(which(-log10(p.values.calibrated.cate) > 4))

plot(density(output.EWAS$beta.p.value))
simu$causal
length(output.EWAS$beta.p.value)
range<-c(1:length(output.EWAS$beta.p.value))
plot(-log10(output.EWAS$beta.p.value)~range)

#store<-fa.pc(simu$Y,7)
