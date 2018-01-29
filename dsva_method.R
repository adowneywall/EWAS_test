
library(dSVA)

### EXAMPLE

### Simulated data from dSVA
output.dSVA<-dSVA(Y = simu$Y,X = simu$X,ncomp = 7)

## tests based on GLMs
p.value.dsva <- NULL
z.score.dsva <- NULL
for (j in 1:1000){
  mod.glm <- glm(simu$Y[,j] ~ ., 
                 data = data.frame(simu$X, output.dSVA$Z),
                 binomial(link = "probit"))
  p.value.dsva[j] <- summary(mod.glm)$coeff[2,4]
  z.score.dsva[j] <- summary(mod.glm)$coeff[2,3] 
}

gif <- median(z.score.dsva^2)/0.456 #genomic inflation factor
p.values.calibrated.dsva <- pchisq(z.score.dsva^2/gif , df = 1, low = F)
hist(p.values.calibrated.dsva)
plot(-log10(p.values.calibrated.dsva)~c(1:1000),xlab="Loci Position")
points(-log10(p.values.calibrated.dsva)[simu$causal]~simu$causal,
       col="red",pch=16)
points(-log10(p.values.calibrated.dsva)[which(-log10(p.values.calibrated.dsva) > 4)]~which(-log10(p.values.calibrated.dsva) > 4),
       col="blue",pch=21,cex=2)
abline(h = 4,col="red")




