#install.packages("RefFreeEWAS")
library(RefFreeEWAS)

#### EXAMPLE ####
data(RefFreeEWAS)


## Sim data
RefFree.Output<-RefFreeEwasModel(t(simu$Y),simu$X,K=5)
dim(RefFree.Output$U)

## Data from RefEWAS
RefFree.Output<-RefFreeEwasModel(rfEwasExampleBetaValues,cbind(1,rfEwasExampleCovariate),K=5)



# z.score.rewas <- NULL
# for (j in 1:1000){
#   mod.glm <- glm(simu$Y[,j] ~ ., 
#                  data = data.frame(simu$X, RefFree.Output$U),
#                  binomial(link = "probit"))
#   p.value.rewas[j] <- summary(mod.glm)$coeff[2,4]
#   z.score.rewas[j] <- summary(mod.glm)$coeff[2,3] 
# }
