setwd("~/Documents/accmass")
source("functions.R")

sink("output")
#load data

print(Sys.time())
print("load data")
data1<-read.csv("raw_data/test2.csv")
#kuja<-read.csv("raw_data/kujawinski2006.csv")


data1[c("C13", "N15", "H2", "O18")]<-0
data1$mz<-H.substract(data1$mz)

#set Signal to Noise limit
SN.lim<-3

#order mz ascending
data1<-data1[order(data1$mz),T]

# remove masses below S:N
data2<-data1[data1$Intensity>SN.lim*data1$Noise,T]

print(Sys.time())
print ("bruteforce")
data2.bf<-bforce(data2$mz, verbose=F)

print(Sys.time())
print ("relations")
data2.rel<-findrelations(data2$mz)
print(Sys.time())
