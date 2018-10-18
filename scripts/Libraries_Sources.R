



## Sources
source("http://bioconductor.org/biocLite.R")
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")
#source("augPairs.R")


#Using biocLite 
biocLite("LEA")
install.packages("LEA_1.4.0_tar.gz", repos = NULL, type ="source")
devtools::install_github("bcm-uga/lfmm")
devtools::install_github('rstudio/shiny')
library(devtools)
library(shiny)

