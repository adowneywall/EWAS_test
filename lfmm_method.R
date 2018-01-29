##### Packages and Dependencies #####
#library(devtools)
devtools::install_github("bcm-uga/lfmm")
source(file = "ewas_generator.R")
## Make sure to install additional dependencies in the event
# that some of h

## example:
## simulation of 1000 beta values for n = 200 individuals (Y)
## We have 3 hidden confounders that will correlate with the 
## simulated phenotype (R2 = 0.8)
## SD for confounders c(.8, .6, .5)
simu <- ewas_generator(n = 200, 
                       p = 1000, 
                       K = 7, 
                       prop.variance = 0.4,
                       sd.U = c(.8, .6, .5,.3,.2,.4,.3),
                       mean.B = .2,
                       sd.B = 0.01,
                       sd.V = 0.1,
                       sigma = .1,
                       setSeed = 5)

## check that R2 = 0.8
summary( lm( simu$X ~ simu$U) )

# methylation data for probe 4
hist(as.numeric(simu$Y[,4]))

# mean values
mean(simu$Y[,4])
simu$freq[4] #truth

## list of causal loci
simu$causal

## run EWAS (LFMM) with true K
mod <- lfmm::lfmm_ridge(simu$Y, simu$X, K = 6)

## tests based on GLMs
p.value <- NULL
z.score <- NULL
for (j in 1:1000){
  mod.glm <- glm(simu$Y[,j] ~ ., 
                 data = data.frame(simu$X, mod$U),
                 binomial(link = "probit"))
  p.value[j] <- summary(mod.glm)$coeff[2,4]
  z.score[j] <- summary(mod.glm)$coeff[2,3] 
}

# plot pvalues
plot(-log10(p.value))

# power (typically too low for detection without calibration)
#which(-log10(p.value) > 3) %in% simu$causal)

# plot calibrated pvalues
gif <- median(z.score^2)/0.456
p.values.calibrated <- pchisq(z.score^2/gif , df = 1, low = F)
hist(p.values.calibrated)
plot(-log10(p.values.calibrated)~c(1:1000),xlab="Loci Position")
#order(simu$causal)
#simu$causal[1]
points(-log10(p.values.calibrated)[simu$causal]~simu$causal,
       col="red",pch=16)
points(-log10(p.values.calibrated)[which(-log10(p.values.calibrated) > 4)]~which(-log10(p.values.calibrated) > 4),
       col="blue",pch=21,cex=2)
abline(h = 4,col="red")

# power
length(which(-log10(p.values.calibrated) > 4) %in% simu$causal)
# Type one  error
length(simu$causal %in% which(-log10(p.values.calibrated) > 4))

ewas.eval.out<-ewas_test(n = 200, 
                         p = 1000, 
                         K = 3, 
                         prop.variance = 0.4,
                         sd.U = c(.8, .6, .5),
                         mean.B = .2,
                         sd.B = 0.01,
                         sd.V = 0.1,
                         sigma = .1,
                         times = 10
)
ewas.eval.out