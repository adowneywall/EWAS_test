#library(devtools)
#source("reFactor.R")

##### DEMO #####
if(!exists("refactor", mode="function")) source("reFactor.R")

library(RCurl)

re.data<-"https://raw.githubusercontent.com/cozygene/refactor/master/demo_files/demo_datafile.txt"
re.cell<-"https://raw.githubusercontent.com/cozygene/refactor/master/demo_files/demo_cellproportions.txt"
re.pheno<-"https://raw.githubusercontent.com/cozygene/refactor/master/demo_files/demo_phenotype.txt"
re.covar<-"https://raw.githubusercontent.com/cozygene/refactor/master/demo_files/demo_covariates.txt"

refactor.data <-read.csv(text=getURL(re.data), 
                         header=T,sep = "")
refactor.pheno<-read.csv(text=getURL(re.pheno),
                         header=F,sep = "")
refactor.covariate<-read.csv(text=getURL(re.covar),
                             header=F,sep = "")
refactor.cellpro<-read.csv(text=getURL(re.cell),
                           header=F,sep = "")

associations_test <- function(O, y, model_append)
{
  observed_pvalues <- c()
  for (site in 1:ncol(O))
  {
    model <- lm(y ~ O[,site] + model_append)
    
    pvalue <- coef(summary(model))[2,4]
    observed_pvalues[site] = as.numeric(pvalue)
  }
  
  return(observed_pvalues) 
}

draw_qqplot <- function(y, title, xtitle, ytitle, style='.')
{
  #x <- runif(length(y))
  x <- seq(1/length(y),1,by = 1/length(y))
  x <- sort(-log10(x))
  y <- sort(-log10(y)) 
  plot(x, y, main=title, xlab=xtitle, ylab=ytitle, pch=style, xlim=c(0,ceiling(max(x))),  ylim=c(0,ceiling(max(y))))
  
  # add y=x trend
  xasix<-0:10
  yasix<-0:10
  abline(lm(xasix~yasix),col=2,lty=1)
}

K = 5                                                     # the number of assumed cell types
# Simulated data:
#DATA_FILE = '../demo_files/demo_datafile.txt'             # methylation levels file path
#PHENO_FILE = '../demo_files/demo_phenotype.txt'           # phenotype file path
#CELL_COMP_FILE = '../demo_files/demo_cellproportions.txt' # cell composition file path


print("Reading methylation data...")
O = as.matrix(refactor.data)    
sample_id_O <- O[1, -1] # extract samples ID
O <- O[-1,] # remove sample ID from matrix
cpgnames <- O[, 1] ## set rownames
O <- O[, -1] 
O = t(matrix(as.numeric(O),nrow=nrow(O),ncol=ncol(O)))

print("Reading phenotype data...")
phenotype_matrix = as.matrix(refactor.pheno)
y <-  matrix(as.numeric(as.matrix(phenotype_matrix[, -1]) ))

print("Reading cell proportion data...")
cellpro = as.matrix(refactor.cellpro)

as.matrix(refactor.cellpro)
# run refactor
output <- refactor(re.data, K,data_fromGithub = T,covarfile = re.covar, covar_fromGithub = T,  out = "demo_refactor")

# run expirements
# init plot
filename <- 'demo_results.png'
png(filename)
par(mfrow=c(2,2))

# exp1
print("Unadjusted analysis...")
observed_pvalues <- associations_test(O,y, matrix(0, nrow = nrow(O), ncol = 1))
draw_qqplot(observed_pvalues, title='Unadjusted analysis', xtitle='-log10(expected)', ytitle='-log10(observed)')

# exp2
print("Adjusted analysis using cell proportions...")
observed_pvalues <- associations_test(O, y, cellpro)
draw_qqplot(observed_pvalues, title='Adjusted analysis using cell proportions', xtitle='-log10(expected)', ytitle='-log10(observed)')

# exp3
print("Adjusted analysis using ReFACTor...")
observed_pvalues <- associations_test(O, y, output$refactor_components[,1:K])
draw_qqplot(observed_pvalues, title='Adjusted analysis using ReFACTor', xtitle='-log10(expected)', ytitle='-log10(observed)')

# exp4
print("Adjusted analysis using PCA...")
observed_pvalues <- associations_test(O, y, output$standard_pca);
draw_qqplot(observed_pvalues, title='Adjusted analysis using PCA', xtitle='-log10(expected)', ytitle='-log10(observed)')

dev.off()
print(paste("plot saved to", filename))

system(paste("open", filename))



#### Running on simple baboon example ####

# run refactor
est.confounder.num(~ V1, p.variable, t(test.chrom), method = "ed")

chrom_5_raw<-"DATA/Empirical_Data/baboon_example/BSSeq_Baboon/betaval_chr1_n50.csv"
test.chrom<-as.matrix(read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/betaval_chr1_n50.txt"))
head(test.chrom)
test.covar<-as.matrix(read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/covariates_n50.txt"))
test.pheno<-as.matrix(read.table("DATA/Empirical_Data/baboon_example/BSSeq_Baboon/predictor_n50.txt"))
dim(test.chrom)
dim(test.covar)

chrom_5_covar<-"DATA/Empirical_Data/baboon_example/BSSeq_Baboon/covariates_n50.txt"
output.baboon <- refactor(chrom_5_raw, K,data_fromGithub = F,covarfile = NULL, covar_fromGithub = F,  out = "demo_refactor")

observed_pvalues <- associations_test(O, y, output$refactor_components[,1:K])
draw_qqplot(observed_pvalues, title='Adjusted analysis using ReFACTor', xtitle='-log10(expected)', ytitle='-log10(observed)')


#### Running simple EWAS generator sim ####

source("reFactor.R")
source(file = "ewas_generator.R")
library(cate)

simu <- ewas_generator(n = 50, 
                       p = 1000, 
                       K = 5, 
                       freq= NULL,
                       prop.variance = 0.4,
                       sd.U = c(.8, .6, .5,.3,.2),
                       mean.B = 5,
                       sd.B = 0.01,
                       sd.V = 0.01,
                       sigma = .1)

print(paste((K.simu<-est.confounder.num(~ X, simu, simu$Y, method = "ed")), " latent factors",
            " found with the 'ed' method",sep = ""))
dim(simu$Y)
dim(simu$X)
output.baboon <- refactor(simu$Y, K.simu,data_raw = T,covarfile = NULL, out = "demo_refactor")
dim(simu$Y)

out.test<-refactor(simu$Y,5,out = "demo_refactor")
?refactor()
