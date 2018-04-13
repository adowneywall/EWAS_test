
# devtools::install_github("dcgerard/vicar")

utils::vignette("sample_analysis", package = "vicar")
#For a sample analysis on a simulated dataset using the vicar functions and other confounder adjustment packages.

library(vicar)
library(dplyr)
library(ggplot2)

# data(sim_gtex)
# Y <- sim_gtex$Y
# X <- sim_gtex$X
# ctl <- sim_gtex$ctl
# which_null <- sim_gtex$which_null
# beta <- sim_gtex$beta

data <- gen.sim.data(n = 100, p = 1000, r = 5)
X.data <- data.frame(X1 = data$X1)

Y <- data$Y
X <- data$X1
ctl <- sim_gtex$ctl
which_null <- sim_gtex$which_null
beta <- sim_gtex$beta


X<-as.matrix(read.table("DATA/example/RNAseq/covariates.txt"))
Y<-as.matrix(read.table("DATA/example/RNAseq/expression.txt",header = T,row.names = 1))
ex<-read.table("DATA/example/RNAseq/expression.txt",header = T)

num_sv     <- sva::num.sv(dat = t(Y), mod = X, method = "be")
num_sv_l   <- sva::num.sv(dat = t(Y), mod = X, method = "leek")
num_sv_bcv <- cate::est.confounder.num(~ Intercept + Treatment, X.data = as.data.frame(X), Y = Y,
                                       method = "bcv", bcv.plot = FALSE, nRepeat = 50)$r
num_sv
num_sv_l
num_sv_bcv