##### Program Description ####
   ## Created by: Olivier Francois
   ## Last Edit: 25-01-2018

##### Parameter Description ####

#n	number of individuals
#p	number of response variables.
#K	number of latent factors.
#freq (vector) mean methylation values (if NULL, set randomly) 
#prop.causal proportion of causal variables (probes/loci).
#prop.variance proportion of phenotypic variance explained by latent structure (intensity of confounding).
#sigma standard deviation of residual errors.
#sd.B	standard deviation for effect sizes (B).
#mean.B	(vector) mean of effect sizes.
#sd.U (vector) standard deviations for factors.
#sd.V	standard deviations for loadings.

ewas_generator <- function (n, 
                            p, 
                            K, 
                            freq = NULL,
                            prop.causal = 0.025,  # proportion of loci that are causal
                            prop.variance = 0.6,  # intensity of confounding
                            sigma = 0.2, # standard deviation of residual errors
                            sd.B = 1.0,  # standard deviation of effect size
                            mean.B = 5.0, # mean effect size
                            sd.U = 1.0, # standard deviation for factors
                            sd.V = 1.0, # sd.V standard deviations for loadings
                            setSeed=NULL) 
{
  set.seed(setSeed)
  outlier <- sample.int(p, prop.causal * p)
  outlier.nb = length(outlier)

  if (is.null(freq)) freq <- runif(p)
  cs <- runif(K, min = -1, max = 1)
  theta <- sqrt(prop.variance/sum((cs/sd.U)^2))
  
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K) # Identity matrix with standard deviation along diagonal
  Sigma <- rbind(Sigma, matrix(cs*theta, nrow = 1))
  Sigma <- cbind(Sigma, matrix(c(cs*theta, 1), ncol = 1))
  covar.mat<-Sigma
  
  UX <- MASS::mvrnorm(n, mu = rep(0, K + 1), Sigma = Sigma)
  U <- UX[, 1:K, drop = FALSE]
  X <- UX[, K + 1, drop = FALSE]
  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))
  B <- matrix(0, p, 1)
  B[outlier, 1] <- rnorm(outlier.nb, mean.B, sd.B)
  
  Epsilon = MASS::mvrnorm(n, mu = rep(0, p), Sigma = sigma^2 * diag(p))
  
  
  Z = U %*% t(V) + X %*% t(B) + Epsilon
  Y = matrix(rep(qnorm(freq),n), nrow = n, byrow = T) + Z
  Y = pnorm(Y)
  
  return(list(Y = Y, #beta values
              X = X, #phenotype
              causal = outlier, #set of causal loci
              covar.mat = covar.mat, #covariance matrix
              U = U, #simulated confounders
              V = V, #loadings
              B = B,  #effect sizes
              freq = freq #mean methylation values
              )
         )
}

# simu <- ewas_generator(n = 10000, 
#                        p = 100, 
#                        K = 5, 
#                        prop.variance = .99,
#                        sd.U = c(.05,.05,.05,.05,.05),
#                        mean.B = 5,
#                        sd.B = 0.01,
#                        sd.V = 0.1,
#                        sigma = .1)
# 
# 
# cov2cor(simu$covar.mat)
# summary(lm( simu$X ~ simu$U))
# library(Rarity)
# #pairs(cbind(simu$U,simu$X))
# corPlot(cbind(simu$U,simu$X),method = "spearman")
# 
# hold<-seq(from=0,to=1,by=.01)
# 
# theta.confound.var<-sqrt(hold/sum((0.1/0.5)^2))
# plot(theta.confound.var~hold,type="l",xlab="degree of confounding",ylab="theta")
# 
# theta.csvar<-sapply(hold,function(x){sqrt(1/sum((x/0.5)^2))})
# plot(theta.csvar~hold,type="l",xlab="value of casual loci",ylab="theta")
# 
# theta.sdvar<-sapply(hold,function(x){sqrt(1/sum((.5/x)^2))})
# plot(theta.csvar~hold,type="l",xlab="sd of casual loci",ylab="theta")
# 
# hold*theta.confound.var
# 
# simu$B
# simu$causal
# 
# mod.glm <- glm(simu$Y[,201] ~ simu$X + simu$U,model = )
# plot(simu$Y[,201] ~ simu$X + simu$U)
# summary(mod.glm)
# simu$B[803]
# simu$V[803,]
# mod.glm2 <- glm(simu$Y[,803] ~ simu$X + simu$U[,c(1,7)],
#                family=binomial(link = "probit"))
# mod.glm2 <- glm(simu$Y[,803] ~ .,
#                 data = data.frame(simu$X,simu$U),
#                 family=binomial(link = "probit"))
# plot(simu$Y[,803]~simu$X)
# summary(mod.glm2)
# 
# (Sigma <- diag(x = sd.U^2, nrow = K, ncol = K))
# Sigma <- rbind(Sigma, matrix(cs*theta, nrow = 1))
# Sigma <- cbind(Sigma, matrix(c(cs*theta, 1), ncol = 1))
# UX <- MASS::mvrnorm(n, mu = rep(0, K + 1), Sigma = Sigma)
# U <- UX[, 1:K, drop = FALSE]

ewas_test<-function(n = 200, 
          p = 1000, 
          K = 3, 
          prop.variance = 0.4,
          sd.U = c(.8, .6, .5),
          mean.B = .2,
          sd.B = 0.01,
          sd.V = 0.1,
          sigma = .1,
          times = 10
          ){
  
  eval.mat<-matrix(,nrow=times,ncol=2)
  for(i in 1:times){
    simu<-ewas_generator(n = n, 
                   p = p, 
                   K = K, 
                   prop.variance = prop.variance,
                   sd.U = sd.U,
                   mean.B = mean.B,
                   sd.B = sd.B,
                   sd.V = sd.V,
                   sigma = sigma)
    mod <- lfmm::lfmm_ridge(simu$Y, simu$X, K = 3)
    z.score <- NULL
    for (j in 1:1000){
      mod.glm <- glm(simu$Y[,j] ~ ., 
                     data = data.frame(simu$X, mod$U),
                     binomial(link = "probit"))
      #p.value[j] <- summary(mod.glm)$coeff[2,4]
      z.score[j] <- summary(mod.glm)$coeff[2,3] 
    }
    # power
    eval.mat[i,1]<-mean(which(-log10(p.values.calibrated) > 4) %in% simu$causal)
    # Type one  error
    eval.mat[i,2]<-mean(!(simu$causal %in% which(-log10(p.values.calibrated) > 4)))
  }
  return(eval.mat)
}


ewas.eval.out
