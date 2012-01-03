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

findrelation<-function(row, lst, mult=1:5, c=0, h=0, o=0, n=0, s=0, p=0, Na=0, c13=0, RE=2E-5) {
  mz<-lst[[1]]
  lst<-lst[[2]]
  rel.mass= c*mass.C12+c13*mass.C13+h*mass.H1+o*mass.O16+n*mass.N14+s*mass.S32+p*mass.P31+Na*mass.Na23
  
  for (i in mult) {
  target.mass=row$mass+rel.mass
  rel.error=rel.mass*RE

  new.row<-row
  new.row$C[1]<-row$C[1]
  new.row$C13[1]<-row$C13[1]
  new.row$H[1]<-row$H[1]
  new.row$O[1]<-row$O[1]
  new.row$N[1]<-row$N[1]
  new.row$S[1]<-row$S[1]
  new.row$P[1]<-row$P[1]
  new.row$Na[1]<-row$Na[1]
  
  new.row$calc.mass<-calc.mass(new.row)
  
  hits<-which(mz>(target.mass-rel.error)&mz<(target.mass+rel.error))
  
  for(j in hits) {
    new.row2<-new.row
    new.row2$mz[1]<-mz[j]
    new.row2$difference<-new.row2$calc.mass-new.row2$mz
    new.row2$diff.ppm<-new.row2$difference/new.row2$mz
    lst[[j]]<-rbind(list[[j]],new.row2)
  }
    
  }
  
  return(list(mz, lst))
}

relationloop<-function(row,  lst, RE=2E-5) {
  
#13c
  print("13c")
  lst<-findrelation(row, lst, n=1:2, c=-1, c13=1)

#CH2
  print("+CH2")
  lst<-findrelation(row, lst, c=1, h=2)

#CH4-O
  print("+CH4 -O")
  lst<-findrelation(row, lst, c=1, h=4, o=-1)

#H2
  print("H2")
  lst<-findrelation(row, lst, h=2)

#C2H4O
  print("+C2H4O")
  lst<-findrelation(row, lst, c=2, h=4, o=1)

#CO2
  print("+CO2")
  lst<-findrelation(row, lst, c=1, o=2)

#C2H2O
  print("+C2H2O")
  lst<-findrelation(row, lst, c=2, h=2, o=1)

#O
  print("+O")
  lst<-findrelation(row, lst, o=1)

#NH
  print("+NH")
  lst<-findrelation(row, lst, n=1, h=1)

return(lst)
}

bforce<-function(mz, FE=3E-6, RE=2E-5, bruteforce.lim=500, verbose=T) {

#  #create result dataframe, set colnames
#  data<-data.frame(matrix(nrow=length(mz), ncol=12))
#  colnames(data)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm")
  
  #add ascending m/z values to results dataframe
  mz<-sort(mz)
  
  #create result list, set colnames
  tmp<-data.frame(matrix(nrow=0, ncol=12))
  colnames(tmp)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm")
  
  lst<-vector("list", length(mz))
  
  for(i in 1:length(mz))
    lst[[i]]<-tmp
  
  #subset of peaks with mz < mz.lim, nr of peaks analysed, initialize counter 
  sel<-which(mz<bruteforce.lim)
  nr<-length(sel)
  count<-1
  
  # for each m/z found
  for(i in sel) {
    print(paste(count,"/",nr,sep=""))
    count<-count+1
        
      #attempt to find the formula using CHON
      tmp<-findformula(mz[i], s.max=0, p.max=0, Na.max=0, c13.max=0)
    
      print("brief BF")
    
      #if no result is found...
      if(nrow(tmp)==0) {
        print("extensive BF")
        
        #try with CHONPSNa
        tmp<-findformula(mz[i], c13.max=0)
      }
    
    #if results were found..
    if(nrow(tmp)>0) {
    #add all results found to the results list
    tmp<-rbind(tmp, lst[[i]])
    #order results and add them to lst
    lst[[i]]<-tmp[order(abs(tmp$diff.ppm)),T]
    }
  }
  # create plot
  if (verbose==T) {
      rm(nr)
      nr<-1
      for (i in 1:length(mz)) {
        nr[i]<-nrow(lst[[i]])
      }
    plot(mz, nr, cex=.5, ylim=c(0,20), ylab="hits", xlab="mass nr")
  }

  return(list(mz, lst))
}

plotresults<-function(lst) {
  mz<-lst[[1]]
  lst<-lst[[2]]
  nr<-1
  for (i in 1:length(mz)) {
    nr[i]<-nrow(lst[[i]])    
  }
  plot(mz, nr, cex=.5, ylim=c(0,20), ylab="hits", xlab="mass nr")
}



res<-bforce(kuja$mz)
lst<-res[[2]]
mz<-res[[1]]
for (i in 1:length(lst)) {
  tmp<-lst[[i]]
  lst[[i]]<-tmp[order(abs(tmp$diff.ppm))]
  print(lst[[i]])
  relationloop(lst[[i]][1,T], res)
}



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

