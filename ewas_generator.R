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

#### EWAS Generator #### 

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


#### Repetative Simulation Generator ####

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

#### Multi Simulation Generator ####

multi.sim.gen(n=c(20,50,100,200),p=c(1000,2000,5000,10000),K=5,prop.variance = seq(from=0.1,to=0.9,by=.1),mean.B=seq(from=1,to=5,by=1))
multi.sim.gen<-function(n, 
                        p, 
                        K, 
                        freq = NA,
                        prop.causal = 0.025,  # proportion of loci that are causal
                        prop.variance = 0.6,  # intensity of confounding
                        sigma = 0.2, # standard deviation of residual errors
                        sd.B = 1.0,  # standard deviation of effect size
                        mean.B = 5.0, # mean effect size
                        sd.U = 1.0, # standard deviation for factors
                        sd.V = 1.0, # sd.V standard deviations for loadings
                        setSeed=NULL,
                        dir.name = "DATA/EWAS_Sims/",
                        sim.folder = "FirstSimTestRun"){
  
  scenarios <- expand.grid(n,p,K,NA,prop.causal,prop.variance,sigma,sd.B,mean.B,NA,sd.V)
  colnames(scenarios)<-c("n","p","K","freq","prop.causal","prop.variance","sigma","sd.b","mean.B","sd.U","sd.V")
  scenarios$Sim<-c(1:nrow(scenarios))
  dir.create(sim.file.update<-paste(dir.name,sim.folder,"_",Sys.Date(),sep = ""))
  write.csv(x=scenarios,file = paste(sim.file.update,"/SimulationParameterList.csv",sep=""))
  
  for(i in 1:nrow(scenarios)){
    temp.scenario<-ewas_generator(n = scenarios$n[i], 
                                  p = scenarios$p[i], 
                                  K = scenarios$K[i], 
                                  freq = NULL,
                                  prop.causal = scenarios$prop.causal[i], 
                                  prop.variance = scenarios$prop.variance[i], 
                                  sigma = scenarios$sigma[i], 
                                  sd.B = scenarios$sd.b[i],  
                                  mean.B = scenarios$mean.B[i], 
                                  sd.U = rnorm(n = scenarios$K[i],mean = 0.4,sd = 0.3), 
                                  sd.V = scenarios$sd.V[i])
    paste(new.folder.name)
    new.folder.name<-paste(sim.file.update,"/Sim",i,sep="")
    paste(new.folder.name)
    dir.create(new.folder.name)
    sim.description<-paste("List of Simulation Parameters \n n = ",scenarios$n[i],
                           "\n p = ", scenarios$p[i],
                           "\n K = ", scenarios$K[i],
                           "\n Beta-value frequency = ", scenarios$freq[i],
                           "\n prop.causal = ",scenarios$prop.causal[i],
                           "\n prop.variance = ",scenarios$prop.variance[i],
                           "\n Sigma = ",scenarios$sigma[i],
                           "\n Standard Deviation of Beta-values = ",scenarios$sd.b[i],
                           "\n Mean of Beta-values = ",scenarios$mean.B[i],
                           "\n Standard Deviation of latent factors (U) = ",scenarios$sd.U[i],
                           "\n Standard Deviation of factor loadings (V) = ",scenarios$sd.V[i],
                           sep = "")
    write.table(x=sim.description,
                file = paste(new.folder.name,"/ParameterList_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".txt",sep=""),
                col.names = F,row.names = F)
    write.csv(x=temp.scenario$Y,file = paste(new.folder.name,"/Y_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
    write.csv(x=temp.scenario$X,file = paste(new.folder.name,"/X_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
    write.csv(x=temp.scenario$U,file = paste(new.folder.name,"/U_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
    write.csv(x=temp.scenario$V,file = paste(new.folder.name,"/V_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
    write.csv(x=temp.scenario$B,file = paste(new.folder.name,"/B_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
    write.csv(x=temp.scenario$covar.mat,file = paste(new.folder.name,"/covar_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
    write.csv(x=temp.scenario$causal, file = paste(new.folder.name,"/causalLoci_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
    print(paste("Simulation run in process... approximately ", round(i/nrow(scenarios)*100,digits=1)," % complete",sep=""))
    }
}

#### Extra Code ####

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
