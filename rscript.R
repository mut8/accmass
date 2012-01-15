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

data2

dist<-dist(data2$mz[1:50])

merge(kuja.bf, kuja.rel, clean=T)

head(data2.bf)
print(Sys.time())
print ("bruteforce")
data2.bf<-bforce(data2$mz, verbose=F)
kuja.bf<-bforce(kuja$mz  , verbose=F)

findformula(110)

kuja.rel<-findrelations(kuja$mz, operations=operations, FE=FE, nmax=5)

kuja.rel[,"id"]

warnings()

length(unique(data2.rel$mz))
length(unique(data2.bf$mz))
length(unique(c(data2.rel$mz, data2.bf$mz))
)
print(Sys.time())
print ("C13 isotopes")
data2.c13.rel<-findrelations(data2$mz, operations=iso.operations, FE=FE, nmax=1)


print ("relations") 
data2.rel<-findrelations(data2$mz, operations=operations, FE=FE, nmax=5)
print(Sys.time())

plot.results(data2.rel.merge)
<-rbind(data2.rel, data2.c13.rel)

data2.merge<-merge(data2.bf, rbind(data2.rel, data2.c13.rel), arg=F, clean=T)
      hplot.results(data2.rel.merge, ylim=c(0,5))
      )

pdf("vank.pdf", width=7, height=14)
par(mfrow=c(2,1))
vanK.plot(data2.merge, cex=0.4, ylim=c(0,1), elements="CNO")
vanK.plot(data2.merge, cex=0.4, ylim=c(0,3))
dev.off()


head(data2.rel)
mz<-data2$mz
rel<-data2.rel.merge
bf<-data2.bf
rel
results
length(tmp2)
clean=F

length(unique(rel$mz))
length(unique(bf$mz))

length(unique(c(rel$mz, bf$mz)))
clean=T

which(is.element(unique(data2.rel$mz), data2$mz)!=T)
which(is.element(
merge(data2$mz, data2.bf, data2.c13.rel, arg=F, clean=T)
  , 
  unique(data2$mz)
  )==F)

merge

output<-merge(mz, bf, rel, arg=F, clean=T)
output
vanK.plot
  results[results$C13==NA,T]
  )
ylim=c()
rel[rel$mz>500,T]

hist(data2.merge$diff.ppm)
plot.results(
  data2.merge[3000:3100,T]
             )
plot(
  output$mz
  , 
  abs(output$diff.ppm)
  )
plot.results(output)

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

bf<-  data2.bf
rel<- rbind(data2.rel, data2.c13.rel)
i<-2
results

mz<-data2$mz
head(mz)

head(
  bf$mz[1]==results$mz[1:6]
  )
head(mz)
i<-1601
cleanup(rows)
rows
results
$calc.mass
kuja.
kuja.
rel$id
i<-83
rel
bf<-kuja.bf
rel<-kuja.rel
sink("dump")
merge<-function(bf, rel, arg=F, clean=T) {
  rel<-rel[rel$id>rel$from.id,T]
  results<-row.template
  nr.unique<-unique(c(bf$id, rel$id))
  results[1:length(nr.unique), "id"]<-nr.unique
  results$id<-sort(nr.unique)
  nr<-length(results$id)
  for (i in 1:nr){
  print(paste(i, "/", nr))
      #print(rel[rel$mz==results$mz[i],T])
      #print(add.rel(results, rel[rel$mz==results$mz[i],T]))
      #print(bf[bf$mz==results$mz[i],T])
      tmp1<-rel[rel$id==results$id[i],T]
      tmp1a<-tmp1[which(is.element(tmp1[,"from.id"], results$id)),T]
        if (nrow(tmp1a)>0)
              for (j in 1:nrow(tmp1a)) {
                  tmp1a[j, elements]<-
                    tmp1a[j,elements]+
                    results[results$id==tmp1a[j,"from.id"], elements]
                  tmp1a[j,T]<-calc.diffs(tmp1a[j,T])
                }
      tmp2<-bf[bf$id==results$id[i],T]
      tmp3<-rbind(tmp1a, tmp2)
      tmp4<-tmp3[is.na(tmp3$calc.mass)!=T,T]
      print(tmp4)
      if (nrow(tmp4)>0){
        print(tmp4)
        tmp5<-points.manipulation(tmp4, arg=arg, clean=clean)
        print(tmp5)
        if (sum(is.na(tmp5$calc.mass)!=T)>0 ) {
          tmp6<-tmp5[which(tmp5$points==max(tmp3$points)),T]
          print(tmp6)
          results[i,2:ncol(results)]<-tmp6[1,2:ncol(tmp6)]
        }
      }
  }
  return(results)
}

i<-1
results
[results$id==79,elements]+tmp1[j,elements]
  tmp1
tmp1a
<-rel[rel$id==results$id[i],T]

i<-84
  

j<-1
results$id
merge(kuja.bf, kuja.rel)
i<-121
results
tmp1

tmp4
tmp5
tmp6
sink()
max(results$C13)

results[2150,T]
results$mz
results$calc.mass
results

h