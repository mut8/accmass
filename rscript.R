options(java.parameters = "-Xmx1000m") 

install.packages("rJava")
install.packages("xlsx")

#functions
##########

##
##H.substract: substracts mass of H+
##
H.substract<-function(x) return(x-mass.H1+mass.e)


##
##findformula: algorithm to calculate sum formula from accurate mass
##

findformula<-function(x) {
  
results<-data.frame(matrix(ncol=6, nrow=0))
colnames(results)<-c("C","H","N","O","calc.mass", "difference")

#C max and C min calculate max/min number of C atoms (and sum of C+N+O atoms)
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
return(results[order(abs(results$difference)),T])
}

#data import
############

setwd("~/Documents/accmass")

#files<-

#load data
setwd("~/Documents/accmass/raw_data")
data1<-read.csv("test.csv")
setwd("~/Documents/accmass")

data1[c("C13", "N15", "H2", "O18")]<-0
data1$mz<-H.substract(data1$mz)

head(data1)

#crop data
#data2<- data1[8:nrow(data1),1:7]
#colnames(data2)<-c("mz", "Intensity", "RelInt", "Resolution", "Charge", "Baseline", "Noise")
  
#order data by m/z
#data3<-data2[order(data2$mz),]

#attach data3

#attach(data3)

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

#x<-240.17595
#x<-275.01668

#data1$mz<-data1$mz-x-mass.H1+mass.e

i

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

testmass<-data1$calc.mass[7069]+mass.C13-mass.C12

testmass<-data1[which(data1$mz>testmass*(1-error) & data1$mz<testmass*(1+error)),"mz"]

data1[9000:9050, T]
[data1$mz>(testmass*(1-error))&data1$mz<(testmass*(1+error)),]
  
tmp<-findformula(data1$mz[7002])
tmp[order(abs(tmp$difference)), T]

H.substract(data1$mz[20001]))

data1
[20000:20200,T]

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
return(results[order(abs(results$difference)),T])
}