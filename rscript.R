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

print(Sys.time())
print ("bruteforce")
data2.bf<-bforce(data2$mz, verbose=F)
#kuja.bf<-bforce(kuja$mz, verbose=F)

print(Sys.time())
print ("relations")
data2.rel<-findrelations(data2$mz)
print(Sys.time())


?class
data2.rel.bak<-data2.rel
data2.bf.bak<-data2.bf

for (i in 1:length(data2.rel))
  
  data2.rel[[i]]<-asnumberize(data2.rel[[i]])

asnumberize<-function(rows) {
  rows[,elements]<-as.numeric(rows[,elements])  
}

cleanup<- function(rows) {
  kill=F
  for (i in 1:nrows(rows)) {
     k1<-min(rows[i,elements])<0
     k2<-rows$C+rows$C13<1
     k3<-rows$H>(2+2*rows$C+2*rows$C13+rows$N)
     k4<-rows$H>3*rows$C
     k5<-rows$O>rows$C
     k6<-rows$N>rows$C
     if(k1|k2|k3|k4|k5|k6) {kill[i]<-T} else {kill[i]<-F}
  }
  return(rows[kill=F,T])
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

data2.merge<-data2.bf
for (i in 1:length(data2.rel)) {
  nrrel<-nrow(data2.rel[[i]])
  print(paste(i, "/", length(data2.rel)))
  if(nrrel>0) {
    for(j in 1:nrrel) {
      nrbf<-nrow(data2.merge[[i]])
      from<-which(data2$mz==data2.rel[[i]][1,"mz.from"])
      if (from < i & nrow(data2.merge[[from]])>0)
      {
        data2.merge[[i]][nrbf+1,T]<-  data2.bf[[from]][1,T]
        #data2.merge[[i]][nrbf+1,T]
        data2.merge[[i]][nrbf+1,c("C", "C13", "H", "O", "N", "P", "S", "Na")]<-
          as.numeric(data2.merge[[i]][nrbf+1,c("C", "C13", "H", "O", "N", "P", "S", "Na")]) +
          as.numeric(data2.rel[[i]][j,c("C", "C13", "H", "O", "N", "P", "S", "Na")])
        data2.merge[[i]][nrbf+1,T]<-calc.mass(data2.merge[[i]][nrbf+1,T])
        data2.merge[[i]][nrbf+1,"comment"]<-data2.rel[[i]][j,"comment"]
        data2.merge[[i]][nrbf+1,"mz"]<-data2$mz[i]
        data2.merge[[i]][nrbf+1,"difference"]<-
          as.numeric(data2.merge[[i]][nrbf+1,"calc.mass"])-
          as.numeric(data2.merge[[i]][nrbf+1,"mz"])
        data2.merge[[i]][nrbf+1,"diff.ppm"]<-
          1E6*as.numeric(data2.merge[[i]][nrbf+1,"difference"])/
          as.numeric(data2.merge[[i]][nrbf+1,"mz"])
 
      }
      
    }
  data2.merge[[i]]<-order.rows(data2.merge[[i]])
  }
}

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
