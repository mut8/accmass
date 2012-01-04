setwd("~/Documents/accmass")
source("functions.R")

#load data
data1<-read.csv("raw_data/test2.csv")
kuja<-read.csv("raw_data/kujawinski2006.csv")


data1[c("C13", "N15", "H2", "O18")]<-0
data1$mz<-H.substract(data1$mz)

SN.lim<-3

data1<-data1[order(data1$mz),T]
data2<-data1[data1$Intensity>SN.lim*data1$Noise,T]

data2.res<-bforce(data2$mz)


res<-bforce(kuja$mz)

operations<-  row.template
  operations[1,T]<-c(NA,-1, 1, 0, 0, 0,0,0,0, NA, NA, NA, "+C13 -C")
  operations[2,T]<-c(NA, 1, 0, 2, 0, 0,0,0,0, NA, NA, NA, "+CH2")
  operations[3,T]<-c(NA, 1, 0, 4, 0,-1,0,0,0, NA, NA, NA, "+CH4 -O")
  operations[4,T]<-c(NA, 0, 0, 2, 0, 0,0,0,0, NA, NA, NA, "+H2")
  operations[5,T]<-c(NA, 2, 0, 4, 0, 1,0,0,0, NA, NA, NA, "+C2H4O")
  operations[6,T]<-c(NA, 1, 0, 0, 0, 2,0,0,0, NA, NA, NA, "+CO2")
  operations[7,T]<-c(NA, 2, 0, 2, 0, 1,0,0,0, NA, NA, NA, "+C2H2O")
  operations[8,T]<-c(NA, 0, 0, 0, 0, 1,0,0,0, NA, NA, NA, "+O")
  operations[9,T]<-c(NA, 1, 0, 0, 1, 0,0,0,0, NA, NA, NA, "+C2H2O")

for (i in 1:nrow(operations))
  operations[i,T]<-calc.mass(operations[i,T])

diffs<-diffs.matrix(data2$mz)

relations.type<-diffs
relations.type[T,T]<-NA

relations.mult<-diffs
relations.mult[T,T]<-NA

nmax=5
mz<-data2$mz

for (k in 1:nrow(operations))
  for (i in 1:nrow(diffs))
    for (j in (i+1):ncol(diffs))
      for (n in 1:nmax) {
      print(paste("k=", k, "/",nrow(operations) , " i=", i,"/",nrow(diffs),  " j=", j,"/",ncol(diffs),  " n=", n, "/", nmax))
      if (
        operations$calc.mass[k] < diffs[i,j] + FE*sqrt(2)*max(mz[i], mz[j]) &
          operations$calc.mass[k] > diffs[i,j] - FE*sqrt(2)*max(mz[i], mz[j]))
      {
        relations.type[i,j]<-operations$comment[k]
        relations.mult[i,j]<-n
      }
    }
    
sum(is.na(relations)==F)


which(
  diffs > (as.numeric(operations$calc.mass[1])*(1-C13E)) & 
  diffs < (as.numeric(operations$calc.mass[1])*(1+C13E)))
      






?rbind

calc.mass(operations)
*RE

calc.mass(operations)

mass.C12

findformula(284.294, Na.max=0, p.max=0, s.max=0, c13.max=0)
kuja[1,T]

lst<-res[[2]]
mz<-res[[1]]

res2<-relationloop(lst[[5]][1,T], res)
findrelation(lst[[2]][1,T], res, mult=1:5, c=1, h=2)

lst[[3]]

for (i in 1:length(lst)) {
  tmp<-lst[[i]]
  if (nrow(tmp)>0) {
  lst[[i]]<-tmp[order(abs(tmp$diff.ppm)),T]
  print(i)
  res<-relationloop(lst[[i]][1,T], res)
  }
}




#data import
############



#files<-


head(data1)

#crop data
#data2<- data1[8:nrow(data1),1:7]
#colnames(data2)<-c("mz", "Intensity", "RelInt", "Resolution", "Charge", "Baseline", "Noise")
  
#order data by m/z
#data3<-data2[order(data2$mz),]

#attach data3

#attach(data3)


formulas2<-1
i<-1

for (i in 1:nrow(kuja)){
  tmp<-reslist[[i]]
  formulas2[i]<-printformula(tmp$C[1],tmp$H[1],tmp$O[1],tmp$N[1])
}

sum(formulas==kuja[,1])/nrow(kuja)

reslist[[1]]

res<-findformula(kuja$mz[1], Na.max=0, s.max=0, p.max=0)
res
reslist<-list(res)

for (i in 1: nrow(kuja)) {
res<-findformula(kuja$mz[i], Na.max=0, s.max=0, p.max=0)
reslist[[i]] <- res1
}

rel<-mass.C12+2*mass.H1
relationloop(res$calc.mass[1], kuja$mz)
reslist
kuja$mz

findformula(100.1)

i

reslist

head(data1[order(data1$relint, decreasing=T),T])

plot(1:nrow(data1), data1$int[order(data1$int, decreasing=T)],  type="l", ylab="intensity", xlab="rank"
     , log="xy"
     )

?log

dev.off()

##loop to analyse all masses

for(i in 1:1000) {

res <- findformula(data1$mz[i])
print(res)

if(nrow(res)!=0)
{
switch<-F
  if(is.na(data1$difference[i])== T) {switch<-T} else 
    if(abs(data1$difference[i]) > abs(res$difference[1])) {switch<-T}
if(switch)
{
data1[i, 4:9]<-res[1,T]
data1[i, c("C13", "H2", "N15", "O18")]<-0
data1$comment[i]<-"direct hit"
}
}
data1$hits[i]<-nrow(res)


 if(nrow(res)>0)
 {
   
   #test for C13 

   testmass<-data1$calc.mass[1]+mass.C13-mass.C12
 for(j in which(data1$mz>testmass*(1-error)&data1$mz<testmass*(1+error)))
 {
  #print("a")
  if(is.na(data1$difference[j])==T) {replace<-T} else 
    if (abs(data1$difference[j]) > abs(data1$mz[j]-testmass)) {replace<-T}
  #print("b")
  if(replace==T)
  {
    #print("c")
    data1[j, 4:9]<-res[1,T]
    data1$C13[j]<-data1$C13[i]+1
  }
 }

   #test for H2   
   
 testmass<-data1$calc.mass[1]+mass.H2-mass.H1
 for(j in which(data1$mz>testmass*(1-error)&data1$mz<testmass*(1+error)))
 {
  #print("a")
  if(is.na(data1$difference[j])==T) {replace<-T} else 
    if (abs(data1$difference[j]) > abs(data1$mz[j]-testmass)) {replace<-T}
  #print("b")
  if(replace==T)
  {
    #print("c")
    data1[j, 4:9]<-res[1,T]
    data1$H2[j]<-data1$H2[i]+1
  }
 } 
 
 }
}

which(data1$mz>testmass*(1-error)&data1$mz<testmass*(1+error))

findformula(data1$mz[7069])
findformula(data1$mz[7069]-mass.C13+mass.C12)



[]

testmass<-data1$calc.mass[7069]+mass.C13-mass.C12

testmass<-data1[which(data1$mz>testmass*(1-error) & data1$mz<testmass*(1+error)),"mz"]

data1[9000:9050, T]
[data1$mz>(testmass*(1-error))&data1$mz<(testmass*(1+error)),]
  
tmp<-findformula(data1$mz[7002])
tmp[order(abs(tmp$difference)), T]

H.substract(data1$mz[20001]))

data1
[20000:20200,T]

