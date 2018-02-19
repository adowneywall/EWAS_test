
library(readr)

start.time.std <- Sys.time()
s.meth<-read.table("DATA/Empirical_Data/GSE82310_M.tab.bz2",sep = " ",nrows = 10000000,header = T)
end.time.std <- Sys.time()     

s.tot<-read.table("DATA/Empirical_Data/GSE82310_Cov.tab.bz2",sep = " ",nrows = 10000,header = T)
s.coord<-read.table("DATA/Empirical_Data/GSE82310_coord.bed.bz2",sep = ".",nrows = 1000000,header = T)

??readr()
?readr()
s.methS<-read_table2("DATA/Empirical_Data/GSE82310_M.tab.bz2",n_max = 1000,col_names = F,skip = 1)
s.methheader<-c("chrom_ID",unlist(read.table("DATA/Empirical_Data/GSE82310_M.tab.bz2",sep = " ",nrows = 1,stringsAsFactors = F)[1,]))


start.time.readr <- Sys.time()
s.methS<-read_table2("DATA/Empirical_Data/GSE82310_M.tab.bz2",n_max = 10000000,col_names = F,skip = 1)
end.time.readr <- Sys.time()     

end.time.readr - start.time.readr
end.time.std - start.time.std
