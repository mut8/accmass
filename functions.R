#we will need the "cluster" package to align masses
require("cluster")

#predefined error margins
error=3E-6
FE=error
RE=2E-5
C13E=1E-4

#exact masses of isotopes (not all currently used)
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


# define row.template and subsets thereoff
row.template<-data.frame(matrix(ncol=18, nrow=0))
colnames(row.template)<-c("id","mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm", "points", "from.id", "relation.type","multiplier", "comment")
elements<-c("C", "H", "O", "N", "P", "S", "Na", "C13")
numerics<-c(elements, "mz", "calc.mass", "difference", "diff.ppm", "points", "from.id", "relation.type", "multiplier")

numberize<-function(rows){
  for (i in numerics) {
    rows[,i]<-as.numeric(rows[,i])
  }
  return(rows)
}

iso.operations<-  row.template
  iso.operations[1,T]<-c(NA, NA,-1, 1, 0, 0, 0,0,0,0, NA, NA, NA, NA, NA,"+C13 -C",1, NA)


operations<-  row.template
  operations[1,T]<-c(NA, NA, 1, 0, 2, 0, 0,0,0,0, NA, NA, NA, NA, NA,"+CH2",1, NA)
  operations[2,T]<-c(NA, NA, 1, 0, 4, 0,-1,0,0,0, NA, NA, NA, NA, NA,"+CH4 -O",1, NA)
  operations[3,T]<-c(NA, NA, 0, 0, 2, 0, 0,0,0,0, NA, NA, NA, NA, NA,"+H2",1, NA)
  operations[4,T]<-c(NA, NA, 2, 0, 4, 0, 1,0,0,0, NA, NA, NA, NA, NA,"+C2H4O",1, NA)
  operations[5,T]<-c(NA, NA, 1, 0, 0, 0, 2,0,0,0, NA, NA, NA, NA, NA,"+CO2",1, NA)
  operations[6,T]<-c(NA, NA, 2, 0, 2, 0, 1,0,0,0, NA, NA, NA, NA, NA,"+C2H2O",1, NA)
  operations[7,T]<-c(NA, NA, 0, 0, 0, 0, 1,0,0,0, NA, NA, NA, NA, NA,"+O",1, NA)
  operations[8,T]<-c(NA, NA, 1, 0, 0, 1, 0,0,0,0, NA, NA, NA, NA, NA,"+C2H2O",1, NA)



#H.substract: substract proton from M, add electron (to correct for positive ionisation)
H.substract<-function(x) return(x-mass.H1+mass.e)


calc.mass<-function(data) {
  data[,c("C", "N", "C13", "O","H", "S",  "P", "Na")]<-as.numeric(data[,c("C", "N", "C13", "O","H", "S",  "P", "Na")])
  data$calc.mass<-as.numeric(data$C*mass.C12+data$N*mass.N14+data$O*mass.O16+data$H*mass.H1+data$C13*mass.C13+data$Na*mass.Na23+data$P*mass.P31+data$S*mass.S32)
  return(data)
}

for (i in 1:nrow(operations))
  operations[i,T]<-calc.mass(operations[i,T])

for (i in 1:nrow(iso.operations))
  iso.operations[i,T]<-calc.mass(iso.operations[i,T])


##
##findformula: algorithm to calculate sum formula from accurate mass
##

x<-210

findformula<-function(x, c.h.min=1/3, c.o.max=1, c13.max=2, c.n.max=1, s.max=3, p.max=3, Na.max=0, screen.error=3E-2, FE=3E-6) {
  
  results<-row.template
  target<-c(x*(1-FE), x*(1+FE))
  
  #C max and C min calculate max/min number of C atoms (and sum of C+N+O atoms)
  
  for(Na in 0:Na.max){
    massleftNa<-x - mass.Na23*Na
    
    for(p in 0:p.max){
      massleftNaP<-massleftNa - mass.P31*p
    
      for(s in 0:s.max){
        massleftNaPS<-massleftNaP - mass.S32 *s
        Cmax<-as.integer(massleftNaPS/12)
        Cmin<-as.integer(massleftNaPS/50)
  
        #print(paste( as.integer(massleft/50), as.integer(massleft/12)))
        #Cmin<-1
        for(c13 in 0:c13.max) { 
        massleftNaPSC13<-massleftNaPS - mass.C13*c13
          
          for(c in Cmin:(Cmax-c13)){
            ctot <- c13 + c
            massleftNaPSC13C<-massleftNaPSC13 - c * mass.C12
            
            for(n in 0:(min(ctot, ceiling(massleftNaPSC13C/12)))){
              massleftNaPSC13CN<-massleftNaPSC13C - mass.N14 * n
              
              for(o in 0:min(ctot, ceiling(massleftNaPSC13CN/14))){
                massleftNaPSC13CNO<-massleftNaPSC13CN - mass.O16*o
                #print(massleftNaPSC13CNO)
                
                for(h in max(floor(ctot/3), floor(massleftNaPSC13CNO*0.8), 0) : max(0, min(((2*ctot+n)+2), ceiling(massleftNaPSC13CNO*1.2)))) {
                  #print(paste("C", c, "H", h, "O", o, "N", n, "S", s, "P", p, "Na", Na))
    
                  #check for N-rule (n+h %% 2 == 0)n
                  if ((n+h) %% 2 == 0){
  
                    #calculate appromimate mass of compound
                    #scr.mass<-c*12+h+n*14+o*16+s*32+p*31+Na*23+c13*13
                    #if(scr.mass>x*(1-screen.error) & scr.mass>x*(1+screen.error)) {
                      
                    #calculate exact mass of compound  
                    calc.mass<-c*mass.C12+n*mass.N14+o*mass.O16+h*mass.H1+c13*mass.C13+Na*mass.Na23+p*mass.P31+s*mass.S32
  
                    #check if calculated mass is within error margin
                    if(calc.mass>target[1] & calc.mass<target[2])
                    {
                      #print result
                      #print(paste("C", c, "C13", c13, "H", h, "N", n, "O", o,"P",p,"S",s,"Na",Na, "m=", calc.mass, "difference", calc.mass-x, "diff.ppm", 1E6*(calc.mass-x)/x))
                      tmp<-nrow(results)
                      results[tmp+1, T]<-c(NA, x, c, c13, h, n, o, s, p, Na, calc.mass, calc.mass-x, 1E6*(calc.mass-x)/x, 1E-6*abs(x/(calc.mass-x)), NA, NA,NA, "bruteforce")
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  results<-numberize(results)
  
  results<-results[order(results$points, decreasing=T),T]
  
  return(results)
}

bforce<-function(mz, FE=3E-6, RE=2E-5, bruteforce.lim=500, verbose=T) {

#  #create result dataframe, set colnames
#  data<-data.frame(matrix(nrow=length(mz), ncol=12))
#  colnames(data)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm")
  
  #add ascending m/z values to results dataframe
  mz<-sort(mz)
  imax<-length(mz)
  
  #create result list, set colnames
  results<-row.template

  #subset of peaks with mz < mz.lim, nr of peaks analysed, initialize counter 
  i<-1
  # for each m/z found
  for(i in 1:imax) {
    print(paste(i,"/",imax,sep=""))
        if(mz[i]<bruteforce.lim){
            
      #print(mz[i])
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
      if (nrow(tmp)>0){
      #add all results found to the results list
      tmp$id<-i
      results<-rbind(results, tmp)
      }
    }
  }
  
  # create plot
  plot.results(results)
  return(results[order(results$mz), T])
}

plot.results<-function(res, ...) {
  tmp<-as.data.frame(table(res$mz))
  plot(tmp[,1], tmp[,2], ylab="hits", xlab="mass nr")
}


findrelations<-function(inp, operations=operations, nmax=5, FE=FE){
inp<-sort(inp)
  
kmax<-nrow(operations)
imax<-length(inp)
jmax<-imax
hits=0  
  
print("calculate diffs")
diffs<-matrix(nrow=length(inp), ncol=length(inp))
for (i in 1:length(inp))
for (j in 1:length(inp)) {
  diffs[i,j]<-inp[i]-inp[j]
}

print("calculate diff error margins")

diffsminerror<-diffs
diffspluerror<-diffs

sr2<-sqrt(2)

for (i in 1:imax) {
print(paste("calculate diff error margins", i, "/", imax))
  for (j in 1:jmax) {

diffsminerror[i,j]<-diffs[i,j]-FE*sr2*max(inp[c(i,j)])
diffspluerror[i,j]<-diffs[i,j]+FE*sr2*max(inp[c(i,j)])
}}

results<-row.template

print("find relations loop")
for(k in 1:nrow(operations))
{
  print(paste("looking for relations:", operations$relation.type[k], "hits:", hits))
  for(n in c(-nmax:-1, 1:nmax)) 
  {
    print(paste("multiplicator", n))
    x<- as.numeric(operations$calc.mass[k])* n
    mat<-diffspluerror > x & diffsminerror < x 
      for (i in 2:(imax))
        for (j in 1:(i-1))
          if (mat[i,j])
          {
            nr<-nrow(results)
            results[(nr+1),c("C","C13","H","N","O","S","P","Na")]<-as.numeric(operations[k,c("C","C13","H","N","O","S","P","Na")])*n
            results[(nr+1),"from.id"]<-j
            results[(nr+1),"mz"]<-inp[i]
            results[(nr+1),"id"]<-i
            results[(nr+1),"multiplier"]<-n*-1
            results[(nr+1),"relation.type"]<-operations[k, "relation.type"]
            results[(nr+1),"comment"]<-paste("relations, m/z difference =", diffs[i,j])
            results[(nr+1),T]<-calc.mass(results[nr+1,T])
            results[(nr+1),"difference"]<-results[nr+1,"calc.mass"]+results[nr+1,"from.id"]-results[nr+1,"mz"]
            results[(nr+1),"diff.ppm"]<-1E6*results[(nr+1),"difference"]/results[(nr+1),"calc.mass"]
            results[nr+1,T]<-numberize(results[nr+1,T])

            results[(nr+2),c("C","C13","H","N","O","S","P","Na")]<-as.numeric(operations[k,c("C","C13","H","N","O","S","P","Na")])*n*-1
            results[(nr+2),"from.id"]<-i
            results[(nr+2),"mz"]<-inp[j]
            results[(nr+2),"id"]<-j
            results[(nr+2),"multiplier"]<-n*-1
            results[(nr+2),"relation.type"]<-operations[k, "relation.type"]
            results[(nr+2),"comment"]<-paste("relations (reverse), difference =", diffs[i,j])
            results[(nr+2),T]<-calc.mass(results[nr+2,T])
            results[(nr+2),"difference"]<-results[nr+2,"calc.mass"]+results[nr+2,"from.id"]-results[nr+2,"mz"]
            results[(nr+2),"diff.ppm"]<-1E6*results[(nr+2),"difference"]/results[(nr+2),"calc.mass"]
            hits<-hits+1
            results[nr+2,T]<-numberize(results[nr+2,T])
          }
  }
}

row.template[,"relation.type"]
"relation.type"

return(results[order(results$mz),T])
}

calc.diffs<-function(row) {
row<-calc.mass(row)
row[,"difference"]<-row[,"calc.mass"]-row[,"mz"]
row[,"diff.ppm"]<-1E6*row[,"difference"]/row[,"mz"]
row[,"points"]<-1/row[,"diff.ppm"]
return(row)
}

points.manipulation<-function(rows, arg=F, clean=F)  {
  if(clean==T) {rows<-cleanup(rows)}
  return(rows)  
}

# bf<-data2.bf
# rel<-data2.c13.rel
# clean=T


vanK.plot<-function(df, ylab="H:C", xlab="O:C", elements="CHO", ...) {
  x <- as.numeric(df$O) / (as.numeric(df$C) + as.numeric(df$C13))
if (elements == "CHO")
  y <- as.numeric(df$H) / (as.numeric(df$C) + as.numeric(df$C13))
if (elements == "CNO")
  y <- as.numeric(df$N) / (as.numeric(df$C) + as.numeric(df$C13))
if (elements == "CNO" & ylab=="H:C")
  ylab<-"N:C"
  
  plot(x, y, ylab=ylab, xlab=xlab, ...)
}
  
  
dbe<-function(row) {
  return(1+row$C+row$C13+row$N/2-row$H/2)
}

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
          tmp6<-tmp5[which(tmp5$points==max(tmp5$points)),T]
          print(tmp6)
          results[i,2:ncol(results)]<-tmp6[1,2:ncol(tmp6)]
        }
      }
  }
  return(results)
}


(results$H-2)/(results$C + results$C13)

cleanup<- function(rows) {
#rows<-rows[-is.na(rows$calc.mass),T]
  kill=F
  for (i in which(is.na(rows$calc.mass)!=T)) {
     k1<-min(rows[i,elements])<0
     k2<-rows$C[i]+rows$C13[i]<1
#     k3<-rows$H[i]>(2+2*rows$C[i]+2*rows$C13[i]+rows$N[i])
     k3<-rows$H[i]>(2+2*rows$C[i]+2*rows$C13[i])
     k4<-rows$H[i]>3*rows$C[i]
     k5<-rows$O[i]>rows$C[i]
     k6<-rows$N[i]>rows$C[i]
     if(k1|k2|k3|k4|k5|k6) {kill[i]<-T} else {kill[i]<-F}
  }
  return(rows[kill==F,T])
}


