
library(readr)

#s.meth<-read.table("DATA/Empirical_Data/GSE82310_M.tab.bz2",sep = " ",nrows = 10000000,header = T)
#s.tot<-read.table("DATA/Empirical_Data/GSE82310_Cov.tab.bz2",sep = " ",nrows = 10000,header = T)
#s.coord<-read.table("DATA/Empirical_Data/GSE82310_coord.bed.bz2",sep = ".",nrows = 1000000,header = T)

s.methheader<-c("chrom_ID",
                unlist(read.table("DATA/Empirical_Data/stickleback_example/GSE82310_M.tab.bz2",sep = " ",
                                  nrows = 1,stringsAsFactors = F)[1,]))
s.meth<-read_table2("DATA/Empirical_Data/stickleback_example/GSE82310_M.tab.bz2",col_names = F,skip = 1)

s.tot<-read_table2("DATA/Empirical_Data/stickleback_example/GSE82310_Cov.tab.bz2",col_names = F,skip = 1)
colnames(s.tot)<-s.methheader

s.coord<-read_table2("DATA/Empirical_Data/stickleback_example/GSE82310_coord.bed.bz2",col_names = F)
colnames(s.coord)<-c("chrom","start","end","p1","num","p2")

chrom.coord<-cbind(s.coord$chrom,s.coord$start)

num.ind<-12
critical.t<-round(num.ind*0.8)
sub.tot<-s.tot[rowSums(s.tot == 0)>=critical.t,]

head(sub.tot)


rowSums(sub.tot == 0)<=critical.t
