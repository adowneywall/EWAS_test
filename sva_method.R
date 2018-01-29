
TOBECOMPLETED



#### Default SVA (Surrogate variable Analysis) methods ####

#source("https://bioconductor.org/biocLite.R")
#biocLite("sva")
#biocLite("bladderbatch")
#biocLite("limma")

#browseVignettes("sva")
##Tutorial
library(sva)
library(bladderbatch)
library(pamr)
library(limma)
library(pcadapt)
#limmaUsersGuide()
data(bladderdata)

## pheno = dataframe with individuals by row and phenotype data
 # Sample - numeric ID for each individual
 # outcome - the outcome of the cancer test 
 # batch -  the sequencing batch number
 # cancer - the determination of the sample being cancerous (variable of interest)
pheno = pData(bladderEset)

## edata = Expression data for each of the 57 individuals at 22283 loci
edata = exprs(bladderEset)


dim(pheno) # 57 individuals
dim(edata) # 22283 alleles x 57 individuals

# For SVA we are comparing the full model (the model including all variables) both the variable of
# interest and the adjustment variables

## Full model
 # variables of interest + adjustment variables
mod = model.matrix(~as.factor(cancer), data=pheno) #in this case we only have cancer (variable of interest)
# NOTE: cancer is included as.factor() here since it has multiple levels

## Null model
 # adjustment variables
mod0 = model.matrix(~1,data=pheno)

### Using the SVA() function on your models
  # Two steps: (1) identify the number of latent factors that need to be estimated.

# Estimating number of latent factors using num.sv()
n.sv = num.sv(edata,mod,method="leek")
n.sv

#Scree plot
pca.data<-prcomp(edata)
plot(temp$sdev^2/sum(temp$sdev^2)*100,
     ylab="Percent of Variance explained",
     xlab="Principal Component",
     type="b")


#Estimate the surrogate variables using sva
svobj = sva(edata,mod,mod0,n.sv=n.sv)

## Calculate parametric F-test p-values for each row (each locus)
pValues = f.pvalue(edata,mod,mod0)
# This
hist(pValues[],breaks = 100,probability = T,xlim = c(0,0.2))

qValues = p.adjust(pValues,method="BH")
hist(qValues[],breaks = 100,probability = T,xlim = c(0,0.2))

## The number of significant p and q-values is still very high 
 # This can be corrected be adjusting the analysis to include
 # the surrogate variables.

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")
hist(qValuesSv[],breaks = 100,probability = T,xlim = c(0,0.2))

fit = lmFit(edata,modSv)
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)
eb = eBayes(fitContrasts)
topTableF(eb, adjust="BH")

batch = pheno$batch

######## With EWAS generated data ####
simu

## Full model
# variables of interest + adjustment variables
mod = model.matrix(~Y, data=simu) #in this case we only have cancer (variable of interest)
# NOTE: cancer is included as.factor() here since it has multiple levels

## Null model
# adjustment variables
mod0 = model.matrix(~1,data=simu)

output.EWAS.sva<-sva(simu$Y,mod,mod0,n.sv=7)
