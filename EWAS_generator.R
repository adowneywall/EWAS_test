##### Program Description ####
   ## Created by: Olivier Francois
   ## Last Edit: 25-01-2018

##### Parameter Description ####
library(dplyr)
library(MASS)
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
ewas_generator <- function (n, # sample#
                            p, # loci#
                            K, # latentfactor#
                            freq = NULL, # vector of mean beta-values 
                            prop.causal = 0.025,  # proportion of loci that are causal
                            prop.variance = 0.6,  # intensity of confounding
                            sigma = 0.2/10, # standard deviation of residual errors
                            sd.B = 1.0,  # standard deviation of effect size
                            mean.B = 5.0, # mean effect size
                            sd.U = 1.0, # standard deviation for factors
                            sd.V = 1.0, # sd.V standard deviations for loadings
                            nb.t = NULL, # estimated negative binomial parameters for total read count
                            MinRD = 10, # minimum number of reads in sim
                            MaxRD = 100, 
                            setSeed=NULL) 
{
  set.seed(setSeed) # this is to set permanent seed
  outlier <- sample.int(p, prop.causal * p) # vector of the 'location' of the causal SMPs
  outlier.nb = length(outlier) # total number of causal SMPs

  if (is.null(freq)) freq <- runif(p) # fills freq vector is user doesn't provide mean b-values
  
  ### Generates the covariance matrix
  cs <- runif(K, min = -1, max = 1) 
  theta <- sqrt(prop.variance/sum((cs/sd.U)^2))
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K) # Identity matrix with standard deviation along diagonal
  Sigma <- rbind(Sigma, matrix(cs*theta, nrow = 1))
  Sigma <- cbind(Sigma, matrix(c(cs*theta, 1), ncol = 1))
  covar.mat<-Sigma
  
  ## Uses multi-variate NORMAL dist to create values for the variable of interest (X), latent factors (U), and laten factor 
   # loadings (V), the effect sizes have an expected value of zero with causal SMPs having an effect stated by the user
  UX <- MASS::mvrnorm(n, mu = rep(0, K + 1), Sigma = Sigma) # variable of interest and factors modeled together for known correlation
  U <- UX[, 1:K, drop = FALSE]
  X <- UX[, K + 1, drop = FALSE]
  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))
  
  # Actual effect of variable of interest
  B <- matrix(0, p, 1) # variable will have zero effect on most loci
  B[outlier, 1] <- rnorm(outlier.nb, mean.B, sd.B) # some small proportion do have an effect (magnitude of effect controlled by user)
  
  # Random effects 
  Epsilon = matrix(rnorm(n*p,sd=sigma),nrow=n)
  
  # Adding fixed factors (X *%* t(B)), 
  # latent factors (U *%* t(V)), 
  # and random effects (Epsilon) to determine DNA methylation (Y)
  Z = U %*% t(V) + X %*% t(B) + Epsilon #add together all matrices of variation
  
  # measure of methylation mean (generally pulled from a distribution of true means (or medians) from 
  # a real data set.)
  M = matrix(rep(qnorm(freq),n) , nrow = n, byrow = T)

  Y.raw = M + Z
  Y = exp(Y.raw)/(1+exp(Y.raw))
  
  # Generated a matrix of same dimension as Y (beta vals), but with 'read counts' simulating the number of reads for 
  # each locus for each individual. Standardized to minimum 10 reads using negative binomial distributions from empirical data  if(!is.null(nb.t)){
  if(!is.null(nb.t)){ 
    t.r<-matrix(nrow = nrow(Y),ncol = ncol(Y))
    m.r<-matrix(nrow = nrow(Y),ncol = ncol(Y))
    s<-NULL
    for(i in 1:p){
      s<-sample(1:nrow(nb.t),size = 1)
      temp.r<-rnbinom(nrow(Y)*10,size=nb.t[s,1],mu=nb.t[s,2])
      while(length((temp.r[temp.r >= MinRD & temp.r <= MaxRD])) < 20){
         s<-sample(1:nrow(nb.t),size = 1)
         temp.r<-rnbinom(nrow(Y),size=nb.t[s,1],mu=nb.t[s,2])
      }
      t.r[,i] <-sample(temp.r[which(temp.r >= MinRD & temp.r <=MaxRD)],nrow(Y),replace=T)
      m.r[,i] <-round(Y[,i]*t.r[,i])
    }
  }else{
    t.r<-NULL
    m.r<-NULL
  }
  return(list(Y = Y, # unordered beta values
              t.reads=t.r, # Total reads for a given locus
              m.reads=m.r, # Methylated reads (determined as negative binomial)
              freq = freq, # Vector of the selected beta-value means
              Epsilon=Epsilon, # Matrix of error
              M = M, # Matrix of mean (per locus) methylation values
              X = X, # Phenotype (variable of interest)
              causal = outlier, # Set of causal loci
              covar.mat = covar.mat, #covariance matrix
              U = U, #simulated confounders
              V = V, #loadings
              B = B, #effect sizes
              freq = freq #mean methylation values
              )
         )
}

collapse<-function(x){
  if(median(x)>= 0.5){
    x[x<0.5]<-(1-x[x<0.5])
  }
  else{
    x[x>=0.5]<-(1-x[x>=0.5])
  }
  return(x)
}

#### EWAS Generator FULL ####
# Contains steps for alternative DNA methylation arrangement. Including ordering mean methylation with sd 
# (i.e. high mean values correspond with high sd). 
ewas_generator_full <- function (n, # sample#
                            p, # loci#
                            K, # latentfactor#
                            freq = NULL, # vector of mean beta-values 
                            prop.causal = 0.025,  # proportion of loci that are causal
                            prop.variance = 0.6,  # intensity of confounding
                            sigma = 0.2/10, # standard deviation of residual errors
                            sd.B = 1.0,  # standard deviation of effect size
                            mean.B = 5.0, # mean effect size
                            sd.U = 1.0, # standard deviation for factors
                            sd.V = 1.0, # sd.V standard deviations for loadings
                            nb.t = NULL, # estimated negative binomial parameters for total read count
                            MinRD = 10, # minimum number of reads in sim
                            MaxRD = 100, 
                            setSeed=NULL) 
{
  set.seed(setSeed) # this is to set permanent seed
  outlier <- sample.int(p, prop.causal * p) # vector of the 'location' of the causal SMPs
  outlier.nb = length(outlier) # total number of causal SMPs
  
  if (is.null(freq)) freq <- runif(p) # fills freq vector is user doesn't provide mean b-values
  
  ### Generates the covariance matrix
  cs <- runif(K, min = -1, max = 1) 
  theta <- sqrt(prop.variance/sum((cs/sd.U)^2))
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K) # Identity matrix with standard deviation along diagonal
  Sigma <- rbind(Sigma, matrix(cs*theta, nrow = 1))
  Sigma <- cbind(Sigma, matrix(c(cs*theta, 1), ncol = 1))
  covar.mat<-Sigma
  
  ## Uses multi-variate NORMAL dist to create values for the variable of interest (X), latent factors (U), and laten factor 
  # loadings (V), the effect sizes have an expected value of zero with causal SMPs having an effect stated by the user
  UX <- MASS::mvrnorm(n, mu = rep(0, K + 1), Sigma = Sigma) # variable of interest and factors modeled together for known correlation
  U <- UX[, 1:K, drop = FALSE]
  X <- UX[, K + 1, drop = FALSE]
  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))
  
  # Actual effect of variable of interest
  B <- matrix(0, p, 1) # variable will have zero effect on most loci
  B[outlier, 1] <- rnorm(outlier.nb, mean.B, sd.B) # some small proportion do have an effect (magnitude of effect controlled by user)
  
  # Random effects 
  #Epsilon = MASS::mvrnorm(n, mu = rep(0, p), Sigma = sigma^2 * diag(p))
  Epsilon = matrix(rnorm(n*p,sd=sigma),nrow=n)
  #dim(Epsilon2)
  # Adding fixed factors (X *%* t(B)), 
  # latent factors (U *%* t(V)), 
  # and random effects (Epsilon) to determine DNA methylation (Y)
  
  ##o.order<-c(1:p)
  Z = U %*% t(V) + X %*% t(B) + Epsilon #add together all matrices of variation
  
  ##rZ<-rev(order(abs(rowMeans(t(Z))))) # order rows of matrix by level of variation
  ##o.order<-o.order[rZ]
  ##Z.ord<-t(t(Z)[rZ,])
  # measure of methylation mean (generally pulled from a distribution of true means (or medians) from 
  # a real data set.)
  M = matrix(rep(qnorm(freq),n) , nrow = n, byrow = T)
  ##rM<-order((abs(rowMeans(t(M)-0.5))))
  ##M.ord<-(t(t(M)[rM,]))
  
  # Both mean values and variation have been order to assign the most variation to SMPs with mean values closer
  # to 0.5. This follows the behavior seen in most empirical datasets.
  
  ##Y.raw = M.ord + Z.ord
  Y.raw = M + Z
  ##Y.raw = Y.raw[,order(o.order)]
  Y = exp(Y.raw)/(1+exp(Y.raw)) 
  #Y = pnorm(Y.raw)
  #Y2 = pnorm(Y2.raw)
  
  # Generated a matrix of same dimension as Y (beta vals), but with 'read counts' simulating the number of reads for 
  # each locus for each individual. Standardized to minimum 10 reads using negative binomial distributions from empirical data  if(!is.null(nb.t)){
  if(!is.null(nb.t)){ 
    t.r<-matrix(nrow = nrow(Y),ncol = ncol(Y))
    m.r<-matrix(nrow = nrow(Y),ncol = ncol(Y))
    s<-NULL
    for(i in 1:p){
      s<-sample(1:nrow(nb.t),size = 1)
      temp.r<-rnbinom(nrow(Y)*10,size=nb.t[s,1],mu=nb.t[s,2])
      while(length((temp.r[temp.r >= MinRD & temp.r <= MaxRD])) < 20){
        s<-sample(1:nrow(nb.t),size = 1)
        temp.r<-rnbinom(nrow(Y),size=nb.t[s,1],mu=nb.t[s,2])
      }
      t.r[,i] <-sample(temp.r[which(temp.r >= MinRD & temp.r <=MaxRD)],nrow(Y),replace=T)
      m.r[,i] <-round(Y[,i]*t.r[,i])
    }
    #m.rY<-apply((t.r*Y),c(1,2),function(x){as.integer(x)})
    #m.rYlog<-apply((t.r*Y),c(1,2),function(x){as.integer(x)})
  }else{
    t.r<-NULL
    m.r<-NULL
    #m.rYlog<-NULL
  }
  #Y.collapse = apply(Y,2,function(x){collapse(x)}) # outdated 
  return(list(#Y = Y, #beta values
    Y = Y, # unordered beta values
    #Y.raw = Y.raw, # beta-values before pnorm adjustment
    #Y.logit = Y.logit,
    t.reads=t.r,
    m.reads=m.r,
    #Y.collapse = Y.collapse, # beta-values collapsed into a 0-0.5 or 1-0.5 range
    freq = freq, # Vector of the selected beta-value means
    Epsilon=Epsilon, 
    M = M,
    X = X, #phenotype
    causal = outlier, #set of causal loci
    covar.mat = covar.mat, #covariance matrix
    U = U, #simulated confounders
    V = V, #loadings
    B = B,  #effect sizes
    Z = Z, #matrix of variation
    freq = freq #mean methylation values
  )
  )
}


#### EWAS Generator 2 ####
ewas_generator2 <- function (n, 
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
  #scaling<-function(x){return(0 + (x[,2] - mean(x[,2])) * (1/sd(x[,2])))}
  #Z.scale<-0 + (sim1$Z[,2] - mean(Z[,2])) * (1/sd(Z[,2]))
  Z.scale<-
  Y = 
  return(list(Y = Y, #beta values
              freq = freq,
              Epsilon=Epsilon,
              X = X, #phenotype
              causal = outlier, #set of causal loci
              covar.mat = covar.mat, #covariance matrix
              U = U, #simulated confounders
              V = V, #loadings
              B = B,  #effect sizes
              Z= Z, #matrix of variation
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
                        sd.Um = 0.4, # Mean of sd
                        sd.Usd = 0.3, # sd of sd 
                        sd.V = 1.0, # sd.V standard deviations for loadings
                        runs = 1, # number of times each sim
                        nb.t = NULL, # size and mu estimates for simulating read counts for beta values
                        setSeed=NULL,
                        rep = 1,
                        dir.name = "DATA/EWAS_Sims/",
                        sim.folder = "FirstSimTestRun"){
  initial=Sys.time()
  sim.file.update<-paste(dir.name,sim.folder,"_",Sys.Date(),sep = "")
  if(file.exists(sim.file.update)){
    evalkey(readkey())
    unlink(sim.file.update,recursive = T)
  }
  scenarios <- expand.grid(n,p,K,NA,prop.causal,prop.variance,sigma,sd.B,mean.B,NA,sd.V,1:rep)
  colnames(scenarios)<-c("n","p","K","freq","prop.causal","prop.variance","sigma","sd.b","mean.B","sd.U","sd.V","Rep")
  scenarios.print<-cbind(Sim=1:nrow(subset(scenarios,scenarios$Rep == 1)),subset(scenarios,scenarios$Rep == max(rep)))
  dir.create(sim.file.update)
  write.csv(x=scenarios.print,file = paste(sim.file.update,"/SimulationParameterList.csv",sep=""),
            row.names = F)
  
  for(i in 1:((nrow(scenarios))/rep)){
    new.folder.name<-paste(sim.file.update,"/Sim",i,sep="")
    print(paste("..Sim",i," of Sim",nrow(subset(scenarios,scenarios$Rep == 1,sep=""))))
    suppressWarnings(dir.create(new.folder.name))
    sim.description<-paste("List of Simulation Parameters \n n = ",scenarios$n[i],
                           "\n p = ", scenarios$p[i],
                           "\n K = ", scenarios$K[i],
                           "\n Beta-value Means = ", scenarios$freq[i],
                           "\n prop.causal = ",scenarios$prop.causal[i],
                           "\n prop.variance = ",scenarios$prop.variance[i],
                           "\n Sigma = ",scenarios$sigma[i],
                           "\n Standard Deviation of Beta-values = ",scenarios$sd.b[i],
                           "\n Mean of Beta-values = ",scenarios$mean.B[i],
                           "\n Standard Deviation of latent factors (U) = ",scenarios$sd.U[i],
                           "\n Standard Deviation of factor loadings (V) = ",scenarios$sd.V[i],
                           "\n Replicate (rep) = ",scenarios$Rep[i],
                           sep = "")
    write.table(x=sim.description,
                file = paste(new.folder.name,"/ParameterList_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".txt",sep=""),
                col.names = F,row.names = F)
    # Complete each replicate of the simulation
    for(j in 1:rep){
      s.sub<-subset(scenarios,scenarios$Rep==j)
      temp.scenario<-ewas_generator(n = s.sub$n[i],
                                    p = s.sub$p[i],
                                    K = s.sub$K[i],
                                    freq = sample(freq,size=s.sub$p[i],replace=T),#ifelse(is.null(freq),freq,sample(freq,size=s.sub$p[i],replace=T)),
                                    prop.causal = s.sub$prop.causal[i],
                                    prop.variance = s.sub$prop.variance[i],
                                    sigma = s.sub$sigma[i],
                                    sd.B = s.sub$sd.b[i],
                                    mean.B = s.sub$mean.B[i],
                                    sd.U = sd.U[1:s.sub$K[i]],#ifelse(is.null(sd.U),rnorm(2,mean=0.4,sd=0.3),sd.U),
                                    sd.V = s.sub$sd.V[i],
                                    nb.t = nb.t)
      rep.fold<-paste(new.folder.name,"/Rep",s.sub$Rep[i],sep="")
      dir.create(rep.fold)
      saveRDS(temp.scenario,
              file = paste(rep.fold,
                           "/",s.sub$n[i],
                           "_",s.sub$p[i],
                           "_",s.sub$prop.variance[i],
                           "_",s.sub$mean.B[i],
                           ".Rdata",sep=""))
      # write.csv(x=temp.scenario$Y,
      #           file = paste(rep.fold,"/Y_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""),
      #           row.names = F)
      # write.csv(x=temp.scenario$Y.logit,
      #           file = paste(rep.fold,"/Ylogit_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""),
      #           row.names = F)
      # write.csv(x=temp.scenario$Y2,
      #           file = paste(rep.fold,"/Y2_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""),
      #           row.names = F)
      # write.csv(x=temp.scenario$X,
      #           file = paste(rep.fold,"/X_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""),
      #           row.names = F)
      # write.csv(x=temp.scenario$U,
      #           file = paste(rep.fold,"/U_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""),
      #           row.names = F)
      # write.csv(x=temp.scenario$V,
      #           file = paste(rep.fold,"/V_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""),
      #           row.names = F)
      # write.csv(x=temp.scenario$B,
      #           file = paste(rep.fold,"/B_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""),
      #           row.names = F)
      # 
      # write.csv(x=temp.scenario$covar.mat,file = paste(rep.fold,"/covar_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""))
      # write.csv(x=temp.scenario$causal, file = paste(rep.fold,"/causalLoci_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""))
      # write.csv(x=temp.scenario$freq, file = paste(rep.fold,"/simMeans_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""))
      # if(!is.null(nb.t)){
      #   write.csv(x=temp.scenario$t.reads,file = paste(rep.fold,"/tReads_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""))
      #   write.csv(x=temp.scenario$m.reads,file = paste(rep.fold,"/mReads_",s.sub$n[i],"_",s.sub$p[i],"_",s.sub$prop.variance[i],"_",s.sub$mean.B[i],".csv",sep=""))
      # }
      
      print(paste("....Rep ",j," of ",rep,sep=""))
    }
    print(paste("Simulation run in process... approximately ", round((i*rep)/nrow(scenarios)*100,digits=1)," % complete...",sep="")) 
    }
}

readkey <- function()
{
  cat("This directory already exists! Do you wish to overwrite it? \n[Y/N] then press [ENTER]")
  return(response <- readline())
}
evalkey<-function(x){
 if(x == "Y"){
   print("File has been overwritten.")
 } 
  else{
    if(x == "N"){
      print("Run aborted.")
      exit()
    }
    else{print("Key stroke not a Y or N, run aborted.");exit()}
  }
}
exit <- function() {
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}


# multi.sim.gen2<-function(n, 
#                         p, 
#                         K, 
#                         freq = NA,
#                         prop.causal = 0.025,  # proportion of loci that are causal
#                         prop.variance = 0.6,  # intensity of confounding
#                         sigma = 0.2, # standard deviation of residual errors
#                         sd.B = 1.0,  # standard deviation of effect size
#                         mean.B = 5.0, # mean effect size
#                         sd.U = 1.0, # standard deviation for factors
#                         sd.Um = 0.4, # Mean of sd
#                         sd.Usd = 0.3, # sd of sd 
#                         sd.V = 1.0, # sd.V standard deviations for loadings
#                         runs = 1, # number of times each sim
#                         setSeed=NULL,
#                         dir.name = "DATA/EWAS_Sims/",
#                         sim.folder = "FirstSimTestRun"){
#   
#   scenarios <- expand.grid(n,p,K,NA,prop.causal,prop.variance,sigma,sd.B,mean.B,NA,sd.V)
#   colnames(scenarios)<-c("n","p","K","freq","prop.causal","prop.variance","sigma","sd.b","mean.B","sd.U","sd.V")
#   scenarios$Sim<-c(1:nrow(scenarios))
#   dir.create(sim.file.update<-paste(dir.name,sim.folder,"_",Sys.Date(),sep = ""))
#   write.csv(x=scenarios,file = paste(sim.file.update,"/SimulationParameterList.csv",sep=""))
#   
#   for(i in 1:(nrow(scenarios))*runs){
#     temp.scenario<-ewas_generator(n = scenarios$n[i], 
#                                   p = scenarios$p[i], 
#                                   K = scenarios$K[i], 
#                                   freq = NULL,
#                                   prop.causal = scenarios$prop.causal[i], 
#                                   prop.variance = scenarios$prop.variance[i], 
#                                   sigma = scenarios$sigma[i], 
#                                   sd.B = scenarios$sd.b[i],  
#                                   mean.B = scenarios$mean.B[i], 
#                                   sd.U = ifelse(is.null(ob),rnorm(2,mean=0.4,sd=0.3),ob), 
#                                   sd.V = scenarios$sd.V[i])
#     
#     new.folder.name<-paste(sim.file.update,"/Sim",i,sep="")
#     dir.create(new.folder.name)
#     sim.description<-paste("List of Simulation Parameters \n n = ",scenarios$n[i],
#                            "\n p = ", scenarios$p[i],
#                            "\n K = ", scenarios$K[i],
#                            "\n Beta-value frequency = ", scenarios$freq[i],
#                            "\n prop.causal = ",scenarios$prop.causal[i],
#                            "\n prop.variance = ",scenarios$prop.variance[i],
#                            "\n Sigma = ",scenarios$sigma[i],
#                            "\n Standard Deviation of Beta-values = ",scenarios$sd.b[i],
#                            "\n Mean of Beta-values = ",scenarios$mean.B[i],
#                            "\n Standard Deviation of latent factors (U) = ",scenarios$sd.U[i],
#                            "\n Standard Deviation of factor loadings (V) = ",scenarios$sd.V[i],
#                            sep = "")
#     write.table(x=sim.description,
#                 file = paste(new.folder.name,"/ParameterList_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".txt",sep=""),
#                 col.names = F,row.names = F)
#     write.csv(x=temp.scenario$Y,file = paste(new.folder.name,"/Y_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
#     write.csv(x=temp.scenario$X,file = paste(new.folder.name,"/X_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
#     write.csv(x=temp.scenario$U,file = paste(new.folder.name,"/U_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
#     write.csv(x=temp.scenario$V,file = paste(new.folder.name,"/V_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
#     write.csv(x=temp.scenario$B,file = paste(new.folder.name,"/B_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
#     write.csv(x=temp.scenario$covar.mat,file = paste(new.folder.name,"/covar_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
#     write.csv(x=temp.scenario$causal, file = paste(new.folder.name,"/causalLoci_",scenarios$n[i],"_",scenarios$p[i],"_",scenarios$prop.variance[i],"_",scenarios$mean.B[i],".csv",sep=""))
#     print(paste("Simulation run in process... approximately ", round(i/nrow(scenarios)*100,digits=1)," % complete",sep=""))
#     library(data.table)
#     fwrite(x=temp.scenario$Y,file="DATA/EWAS_Sims/Hold",quote='auto')
#     write.csv(x=temp.scenario$Y,file="DATA/EWAS_Sims/Hold")
#   }
# }
