#n	number of individuals
#p	number of response variables.
#K	number of latent factors.
#freq (vector) mean methylation values (if NULL, set randomly) 
#prop.causal proportion of causal variables (probes/loci).
#prop.variance proportion of phenotypic variance explained by latent structure (intensity of confounding).
#sigma.error	standard deviation of residual errors.
#sd.B	standard deviation for effect sizes (B).
#mean.B	(vector) mean of effect sizes.
#sd.U (vector) standard deviations for factors.
#sd.V	standard deviations for loadings.

ewas_generator <- function (n, 
                            p, 
                            K, 
                            freq = NULL,
                            prop.causal = 0.025, 
                            prop.variance = 0.6, 
                            sigma = 0.2, 
                            sd.B = 1.0, 
                            mean.B = 5.0, 
                            sd.U = 1.0, 
                            sd.V = 1.0) 
{
  outlier <- sample.int(p, prop.causal * p)
  outlier.nb = length(outlier)

  if (is.null(freq)) freq <- runif(p)
  
  cs <- runif(K, min = -1, max = 1)
  theta <- sqrt(prop.variance/sum((cs/sd.U)^2))
  
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)
  Sigma <- rbind(Sigma, matrix(cs*theta, nrow = 1))
  Sigma <- cbind(Sigma, matrix(c(cs*theta, 1), ncol = 1))
  
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
              U = U, #simulated confounders
              V = V, #loadings
              B = B,  #effect sizes
              freq = freq #mean methylation values
              )
         )
}

## example:
## simulation of 1000 beta values for n = 200 individuals (Y)
## We have 3 hidden confounders that will correlate with the 
## simulated phenotype (R2 = 0.8)
## SD for confounders c(.8, .6, .5)
simu <- ewas_generator(n = 200, 
                       p = 1000, 
                       K = 3, 
                       prop.variance = 0.8,
                       sd.U = c(.8, .6, .5),
                       mean.B = .2,
                       sd.B = 0.01,
                       sd.V = 0.1,
                       sigma = .1)

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
mod <- lfmm::lfmm_ridge(simu$Y, simu$X, K = 3)

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


# power
which(-log10(p.value) > 3) %in% simu$causal

# plot calibrated pvalues
gif <- median(z.score^2)/0.456
p.values.calibrated <- pchisq(z.score^2/gif , df = 1, low = F)
hist(p.values.calibrated)
plot(-log10(p.values.calibrated))

# power
which(-log10(p.values.calibrated) > 4) %in% simu$causal





