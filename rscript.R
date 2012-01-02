#options(java.parameters = "-Xmx1000m") 

#install.packages("rJava")
#install.packages("xlsx")

#functions
##########

##
##H.substract: substracts mass of H+
##
H.substract<-function(x) return(x-mass.H1+mass.e)


##
##findformula: algorithm to calculate sum formula from accurate mass
##

calc.mass<-function(data) {
  print(data$C*mass.C12+data$N*mass.N14+data$O*mass.O16+data$H*mass.H1+data$C13*mass.C13+data$Na*mass.Na23+data$P*mass.P31+data$S*mass.S32)
}

findformula<-function(x, c.h.min=1/3, c.o.max=1, c13.max=2, c.n.max=1, s.max=3, p.max=3, Na.max=0, screen.error=3E-2, FE=3E-6) {
  
results<-data.frame(matrix(ncol=12, nrow=0))
colnames(results)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm")

#C max and C min calculate max/min number of C atoms (and sum of C+N+O atoms)
Cmax<-as.integer(x / (mass.C12 + mass.H1*c.h.min))
Cmin<-as.integer((x - 2*mass.H1) / (mass.C12 + 2* mass.H1 + mass.O16))
for(c in Cmin:Cmax)
for(c13 in 0:c13.max)
for(n in 0:max(0,min ((Cmax-c+1), c*c.n.max)))
for(o in 0:(c*c.o.max))
for(h in as.integer(((c13+c)/3)):(2*(c13+c)+2))
for(p in 0:p.max)
for(Na in 0:Na.max)
for(s in 0:s.max)
  
  #check for N-rule (n+h %% 2 == 0)
if (n+h %% 2 == 0){


#calculate appromimate mass of compound
#scr.mass<-c*12+h+n*14+o*16+s*32+p*31+Na*23+c13*13
#if(scr.mass>x*(1-screen.error) & scr.mass>x*(1+screen.error)) {
  
#calculate exact mass of compound  
calc.mass<-c*mass.C12+n*mass.N14+o*mass.O16+h*mass.H1+c13*mass.C13+Na*mass.Na23+p*mass.P31+s*mass.S32

#check if calculated mass is within error margin
if(calc.mass>x-(x*FE) & calc.mass<x+(x*FE))
{
  #print result
  print(paste("C", c, "C13", c13, "H", h, "N", n, "O", o,"P",p,"S",s,"Na",Na, "m=", calc.mass, "difference", calc.mass-x, "diff.ppm", 1E6*(calc.mass-x)/x))
  results<-rbind(results, c(x, c, c13, h, n, o, s, p, Na,calc.mass, calc.mass-x, 1E6*(calc.mass-x)/x))

}
#}
}  
colnames(results)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm")
return(results[order(abs(results$difference)),T])
}

findrelation<-function(mass, rel.mass, data, RE=2E-5) {
  target.mass=mass+rel.mass
  rel.error=rel.mass*RE
  print(target.mass)
  print(which(data$mz>(target.mass-rel.error)&data$mz<(target.mass+rel.error)))
  return(which(data$mz>(target.mass-rel.error)&data$mz<(target.mass+rel.error)))
}

relationloop<-function(row, data, list, RE=2E-5) {

print(data)  
  
#13c
  print("13c")
dif=mass.C13-mass.C12
for (i in 1:2){
tmp<-findrelation(data$mz[row], i*dif, data)
if(length(tmp)>0) {
  print(data[row,T])
  row<-data[row,T]

  row$C13[1]<-data$C13[row]+i
  row$C[1]<-data$C[row]-i
  row$calc.mass[1]<-calc.mass(row)
  print(row)
  
for(j in tmp) {
  row1<-row
  row1$mz[1]<-data$mz[j]
  row1$difference[1]<-row1$calc.mass-row1$mz
  row1$dif.ppm[1]<-10E6*(row1$calc.mass-row1$mz)/row1$mz
  list[[j]]<-rbind(list[[j]], row1)

}
}
}

#CH2
  print("+CH2")
dif=mass.C12+2*mass.H1
for (i in 1:5){
tmp<-findrelation(data$mz[row], i*dif, data)
if(length(tmp)>0) {
  row<-data[row,T]
  print(row)
  row$C[1]<-data[row,"C"]+i
  row$H[1]<-data[row,"H"]+1*i
  row[1,10]<-calc.mass(row)
    print(row)
  
for(j in tmp) {
  row1<-row
  row1$mz[1]<-data$mz[j]
  row1$difference[1]<-row1$calc.mass-row1$mz
  row1$diff.ppm[1]<-10E6*(row1$calc.mass-row1$mz)/row1$mz
  list[[j]]<-rbind(list[[j]], row1)
}
}
}


#CH4-O
  print("+CH4 -O")
dif=mass.C12+4*mass.H1-mass.O16
for (i in 1:5){
tmp<-findrelation(data$mz[row], i*dif, data)
if(length(tmp)>0) {
  row<-data[row,T]
  row$C[1]<-data$C[row]+i
  row$H[1]<-data$H[row]+2*i
  row$O[1]<-data$O[row]-i
  row$calc.mass[1]<-calc.mass(row)
  
for(j in tmp) {
  row1<-row
  row1$mz[1]<-data$mz[j]
  row1$difference[1]<-row1$calc.mass-row1$mz
  row1$diff.ppm[1]<-10E6*(row1$calc.mass-row1$mz)/row1$mz
  list[[j]]<-rbind(list[[j]], row1)

}
}
}



#H2
  print("H2")
dif=2*mass.H1
for (i in 1:5){
tmp<-findrelation(data$mz[row], i*dif, data)
if(length(tmp)>0) {
  row<-data[row,T]
  row$H[1]<-data$H[row]+2*i
  row$calc.mass[1]<-calc.mass(row)
  
for(j in tmp) {
  row1<-row
  row1$mz[1]<-data$mz[j]
  row1$difference[1]<-row1$calc.mass[1]-row1$mz[1]
  row1$diff.ppm[1]<-10E6*(row1$calc.mass[1]-row1$mz[1])/row1$mz[1]
  list[[j]]<-rbind(list[[j]], row1)
}
}
}

#C2H4O
    print("+C2H4O")
dif=2*mass.C12+4*mass.H1+mass.O16
for (i in 1:5){
tmp<-findrelation(data$mz[row], i*dif, data)
if(length(tmp)>0) {
  row<-data[row,T]
  row$C[1]<-data$C[row]+2*i
  row$H[1]<-data$H[row]+4*i
  row$O[1]<-data$O[row]+i
  row$calc.mass[1]<-calc.mass(row)
  
for(j in tmp) {
  row1<-row
  row1$mz[1]<-data$mz[j]
  row1$difference[1]<-row1$calc.mass[1]-row1$mz[1]
  row1$diff.ppm[1]<-10E6*(row1$calc.mass[1]-row1$mz[1])/row1$mz[1]
  list[[j]]<-rbind(list[[j]], row1)
}
}
}
  
#CO2
    print("+CO2")
dif=mass.C12+2*mass.O16
for (i in 1:5){
tmp<-findrelation(data$mz[row], i*dif, data)
if(length(tmp)>0) {
  row<-data[row,T]
  row$C[1]<-data$C[row]+i
  row$O[1]<-data$O[row]+2*i
  row$calc.mass[1]<-calc.mass(row)
  
for(j in tmp) {
  row1<-row
  row1$mz[1]<-data$mz[j]
  row1$difference[1]<-row1$calc.mass[1]-row1$mz[1]
  row1$diff.ppm[1]<-10E6*(row1$calc.mass[1]-row1$mz[1])/row1$mz[1]
  list[[j]]<-rbind(list[[j]], row1)
}
}  
}  

#C2H2O
    print("+C2H2O")
dif=2*mass.C12+2*mass.H1+mass.O16
for (i in 1:5){
tmp<-findrelation(data$mz[row], i*dif, data)
if(length(tmp)>0) {
  row<-data[row,T]
  row$C[1]<-data$C[row]+2*i
  row$H[1]<-data$H[row]+4*i
  row$O[1]<-data$O[row]+i
  row$calc.mass[1]<-calc.mass(row)
  
for(j in tmp) {
  row1<-row
  row1$mz[1]<-data$mz[j]
  row1$difference[1]<-row1$calc.mass[1]-row1$mz[1]
  row1$diff.ppm[1]<-10E6*(row1$calc.mass[1]-row1$mz[1])/row1$mz[1]
  list[[j]]<-rbind(list[[j]], row1)

}
}
}

#O
    print("+O")
dif=mass.O16
for (i in 1:5){
tmp<-findrelation(data$mz[row], i*dif, data)
if(length(tmp)>0) {
  row<-data[row,T]
  row$O[1]<-data$O[row]+2*i
  row$calc.mass[1]<-calc.mass(row)
  
for(j in tmp) {
  row1<-row
  row1$mz[1]<-data$mz[j]
  row1$difference[1]<-row1$calc.mass[1]-row1$mz[1]
  row1$diff.ppm[1]<-10E6*(row1$calc.mass[1]-row1$mz[1])/row1$mz[1]
  list[[j]]<-rbind(list[[j]], row1)

}
}
}

#NH
    print("+NH")
dif=mass.N14+mass.H1
for (i in 1:5){
tmp<-findrelation(data$mz[row], i*dif, data)
if(length(tmp)>0) {
  row<-data[row,T]
  row$N[1]<-data$N[row]+i
  row$H[1]<-data$H[row]+i
  row$calc.mass[1]<-calc.mass(row)
  
for(j in tmp) {
  row1<-row
  row1$mz[1]<-data$mz[j]
  row1$difference[1]<-row1$calc.mass[1]-row1$mz[1]
  row1$diff.ppm[1]<-10E6*(row1$calc.mass[1]-row1$mz[1])/row1$mz[1]
  list[[j]]<-rbind(list[[j]], row1)

}
}
}
print(list)
return(list)
}

mainloop<-function(mz, FE=3E-6, RE=2E-5, bruteforce.lim=500) {

  data<-data.frame(matrix(nrow=length(mz), ncol=12))
  colnames(data)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm")
  data$mz<-sort(mz)
  tmp<-data.frame(matrix(nrow=0, ncol=12))
  colnames(tmp)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm")
  plot(1:nrow(data), 1:nrow(data), type="n", ylim=c(0,20), ylab="hits", xlab="mass nr")
  for(i in 1:length(mz))
    list[[i]]<-tmp

  for(i in 1:length(mz)) {
    print(paste("check-check",i))
    if(data$mz[i]<bruteforce.lim) {
      print(paste("check-check",i,"a"))
      tmp<-findformula(data$mz[i], s.max=0, p.max=0, Na.max=0)
      print(paste("check-check",i,"b"))
      if(nrow(tmp)==0) {
        print("check long run")
        tmp<-findformula(data$mz[i])
      }
    list[[i]]<-rbind(tmp, list[[i]])
    }
    tmp<-list[[i]]
    print(paste("check-check",i,"c"))
    print(tmp[1,1:12])
    tmp<-tmp[order(abs(tmp$diff.ppm))]
    print(class(tmp))
    print(tmp)
    data[i,1:12]<-tmp[1,1:12]
    list<-relationloop(i, data, list)
    print("check")
    points(i, nrow(list[[i]]))    
  }
                  
                   
  return(list(data,list))                                    
}

findformula(  kuja$mz[2])
            , s.max=0, p.max=0, Na.max=0)
mainloop(kuja$mz)

#data import
############

setwd("~/Documents/accmass")

#files<-

#load data
setwd("~/Documents/accmass/raw_data")
data1<-read.csv("test.csv")
kuja<-read.csv("kujawinski2006.csv")

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
mass.Na23<-22.9897697
mass.P31<-30.9737634
mass.S32<-31.9720718

#MS error margin in
error=3E-6

#x<-240.17595
#x<-275.01668

#data1$mz<-data1$mz-x-mass.H1+mass.e

printformula<-function(c,h,o,n) { return (paste("C",c,"H",h,"O",o,"N",n, sep=""))}

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

