setwd("~/Documents/accmass")
source("functions.R")


#sink("output")
#load data

print(Sys.time())
print("load data")
data1<-read.csv("raw_data/test2.csv")
kuja<-read.csv("raw_data/kujawinski2006.csv")

data1[c("C13", "N15", "H2", "O18")]<-0
data1$mz<-H.substract(data1$mz)

#set Signal to Noise limit
SN.lim<-3

#order mz ascending
data1<-data1[order(data1$mz),T]

# remove masses below S:N
data2<-data1[data1$Intensity>SN.lim*data1$Noise,T]

dist<-dist(data2$mz[1:50])

merge(kuja.bf, kuja.rel, clean=T)


print(Sys.time())
print ("bruteforce")
data2.bf<-bforce(data2$mz, verbose=F)
#kuja.bf<-bforce(kuja$mz, verbose=F)


kuja.rel<-findrelations(kuja$mz, operations=operations, FE=FE, nmax=5)

warnings()

print(Sys.time())
print ("C13 isotopes")
data2.c13.rel<-findrelations(data2$mz, operations=iso.operations, FE=FE, nmax=1)


print ("relations") b
data2.rel<-findrelations(data2$mz, operations=operations, FE=FE, nmax=5)
print(Sys.time())

data2.rel.merge<-rbind(data2.rel, data2.c13.rel)

data2.merge<-merge(data2.bf, rbind(data2.rel, data2.c13.rel), arg=F, clean=F)
      hplot.results(data2.rel.merge, ylim=c(0,5))
      )

pdf("vank.pdf", width=7, height=14)
par(mfrow=c(2,1))
vanK.plot(data2.merge, cex=0.4, ylim=c(0,1), elements="CNO")
vanK.plot(data2.merge, cex=0.4, ylim=c(0,3))
dev.off()

hist(data2.merge$diff.ppm)
plot.results(data2.merge)
plot(data2.merge$mz, abs(data2.merge$diff.ppm))

?class
data2.rel.bak<-data2.rel
data2.bf.bak<-data2.bf

for (i in 1:length(data2.rel))
  
  data2.rel[[i]]<-asnumberize(data2.rel[[i]])

asnumberize<-function(rows) {
  rows[,elements]<-as.numeric(rows[,elements])  
}


  
  
}

#data2.rel<-data2.rel.bak
#data2.bf<-data2.bf.bak

#order.rows(data2.bf[[606]])


  
nrrel
nrbf
from
data2.rel[[i]]
i<-606
from
j<-1
j

data2.rel[[606]]
data2.merge[[606]]

d
warnings()

data2.bf
 <-data2.bf[[2]]
data2.merge
<- data2.bf
data2.bf
  rm(
    data2.merge[[606]]
    )
data2.rel[[1]]

plotresults(data2.merge, data2$mz)
  data2.rel
            )

print(Sys.time())


kuja.bf.col<-collect(kuja.bf)
vanK.plot(kuja.bf.col, cex=.5, xlim=c(0,1), ylim=c(0,3), elements="CHO")




dist<-dist(log(data2$mz[1:50]))
clust<-as.hclust(agnes(dist))
cutree(clust, h=log(1E-1))