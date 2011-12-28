options(java.parameters = "-Xmx1000m") 

install.packages("rJava")
install.packages("xlsx")

require("xlsx")

setwd("~/R/angelika")

#files<-

#load data
setwd("~/R/angelika/raw_data")
data1<-read.csv("test.csv")
setwd("~/R/angelika")

data1[c("C13", "N15", "H2", "O18")]<-0

head(data1)

#crop data
data2<- data1[8:nrow(data1),1:7]
colnames(data2)<-c("mz", "Intensity", "RelInt", "Resolution", "Charge", "Baseline", "Noise")
  
#order data by m/z
data3<-data2[order(data2$mz),]

#attach data3

attach(data3)

#accurate masses
mass.e<-0.0005485
mass.H1<-1.007825037
mass.H2<-2.014101787
mass.C12<-12
mass.C13<-13.00335484
mass.N14<-14.00307401   
mass.N15<-15.00010898
mass.O16<-15.99491464
mass.O18<-17.99915939


#MS error margin in
error=3E-6

x<-240.17595
x<-275.01668

#data1$mz<-data1$mz-x-mass.H1+mass.e
H.substract<-function(x) return(x-mass.H1+mass.e)

i

for(i in 1:70) {

res <- findformula(data1$mz[i])
print(res)



if(is.na(data1$difference[i])== F)
  if(data1$difference[i] > abs(res$difference))

data1[i, 4:9]<-res[1,T]
data1$hits[i]<-nrow(res)
data1$comment[i]<-"direct hit"


if(nrow(res)>0)
{
testmass<-data1$calc.mass[1]+mass.C13-mass.C12
data1[data1$mz>(testmass*(1-error))&data1$mz<(testmass*(1+error)), 4:9]<-res[1,T]
data1[data1$mz>(testmass*(1-error))&data1$mz<(testmass*(1+error)), 10]<-(data1$C13[i]+1)
}
}

data1
[data1$mz>(testmass*(1-error))&data1$mz<(testmass*(1+error)),]
  
findformula(H.substract(data1$mz[20001]))

data1[20000:20200,T]

findformula<-function(x) {

results<-data.frame(matrix(ncol=6, nrow=0))
colnames(results)<-c("C","H","N","O","calc.mass", "difference")

Cmax<-as.integer(x / (mass.C12 + mass.H1/3))
Cmin<-as.integer((x - 2*mass.H1) / (mass.C12 + 2* mass.H1 + mass.O16))
for(c in Cmin:Cmax)
for(n in 1:min((Cmax-c+1), c))
for(o in as.integer((Cmax-c-n+1)/2):min((Cmax-c-n+1),c))
for(h in as.integer((c/3)):(2*c+2))  
{

#check for N-rule (n+h %% 2 == 0)
if (n+h %% 2 == 0){}
{
  
#calculate mass of compound  
calc.mass<-c*mass.C12+n*mass.N14+o*mass.O16+h*mass.H1

#check if calculated mass is within error margin
if(calc.mass>x-(x*error) & calc.mass<x+(x*error)){
{
  #print result
  print(paste("C", c, "H", h, "N", n, "O", o, "m=", calc.mass, "difference", calc.mass-x))
  results<-rbind(results, c(c, h, n, o, calc.mass, calc.mass-x))

}
}
} 
}  
colnames(results)<-c("C","H","N","O","calc.mass", "difference")
return(results[order(results$difference),])
}