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
numerics<-c(elements, "mz", "calc.mass", "difference", "diff.ppm", "points", "from.id",  "multiplier")

numberize<-function(rows){
  for (i in numerics) {
    rows[,i]<-as.numeric(rows[,i])
  }
  return(rows)
}  

operations<-  row.template
  operations[1,T]<-c(NA, NA, 1, 0, 2, 0, 0,0,0,0, NA, NA, NA, NA, NA,"+CH2",5, NA)
  operations[2,T]<-c(NA, NA, 1, 0, 4, 0,-1,0,0,0, NA, NA, NA, NA, NA,"+CH4 -O",5, NA)
  operations[3,T]<-c(NA, NA, 0, 0, 2, 0, 0,0,0,0, NA, NA, NA, NA, NA,"+H2",5, NA)
  operations[4,T]<-c(NA, NA, 2, 0, 4, 0, 1,0,0,0, NA, NA, NA, NA, NA,"+C2H4O",5, NA)
  operations[5,T]<-c(NA, NA, 1, 0, 0, 0, 2,0,0,0, NA, NA, NA, NA, NA,"+CO2",5, NA)
  operations[6,T]<-c(NA, NA, 2, 0, 2, 0, 1,0,0,0, NA, NA, NA, NA, NA,"+C2H2O",5, NA)
  operations[7,T]<-c(NA, NA, 0, 0, 0, 0, 1,0,0,0, NA, NA, NA, NA, NA,"+O",5, NA)
  operations[8,T]<-c(NA, NA, 1, 0, 0, 1, 0,0,0,0, NA, NA, NA, NA, NA,"+C2H2O",5, NA)
  operations[9,T]<-c(NA, NA,-1, 1, 0, 0, 0,0,0,0, NA, NA, NA, NA, NA,"+C13 -C",1, NA)



#H.substract: substract proton from M, add electron (to correct for positive ionisation)
H.substract<-function(x) return(x-mass.H1+mass.e)


calc.mass<-function(data) {
  data[,c("C", "N", "C13", "O","H", "S",  "P", "Na")]<-as.numeric(data[,c("C", "N", "C13", "O","H", "S",  "P", "Na")])
  data$calc.mass<-as.numeric(data$C*mass.C12+data$N*mass.N14+data$O*mass.O16+data$H*mass.H1+data$C13*mass.C13+data$Na*mass.Na23+data$P*mass.P31+data$S*mass.S32)
  return(data)
}

for (i in 1:nrow(operations))
  operations[i,T]<-calc.mass(operations[i,T])




##
##findformula: algorithm to calculate sum formula from accurate mass
##

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

bforce<-function(data, FE=3E-6, RE=2E-5, bruteforce.lim=500, verbose=T, clean=T, arg=F) {

#  #create result dataframe, set colnames
#  data<-data.frame(matrix(nrow=length(mz), ncol=12))
#  colnames(data)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm")
  
  #add ascending m/z values to results dataframe
  imax<-length(data$mz)
  
  #create result list, set colnames
  results<-row.template

  #subset of peaks with mz < mz.lim, nr of peaks analysed, initialize counter 
  # for each m/z found
  #for(i in 1:179) {
  for(i in 1:imax) {
    print(paste(i,"/",imax,sep=""))

      tmp<-row.template      
    
        #ad all results found to the results list    print(paste(i,"/",imax,sep=""))
        if(data$mz[i]<bruteforce.lim){
      
        #print(mz[i])
          #attempt to find the formula using CHON
          print("brief BF")
          tmp<-findformula(data$mz[i], s.max=0, p.max=0, Na.max=0, c13.max=0)

          #if no result is found...
          if(nrow(tmp)==0) {
            print("extensive BF")
            
            #try with CHONPSNa
            tmp<-findformula(data$mz[i], c13.max=0, Na.max=0)
          }
        
        #ad all results found to the results list
      
            if(nrow(tmp)>0)
        tmp$id<-data$id[i]
      }
    if (verbose==T) {print(results)}
    tmp2 <- 
      findrelations(data[i, T], results,operations=operations, FE=FE)
    if (verbose==T) {print(tmp2)}
    tmp3<-
      points.manipulation(rbind(tmp2, tmp), arg=arg, clean=clean)
    if (nrow(tmp3) > 0) {
    tmp4<-tmp3[which(tmp3$points==max(tmp3$points)),T]
      if (nrow(tmp4) > 0) {
        nr<-nrow(results) 
        results[nr+1,T]<-tmp4[1,T] 
        rm(tmp4)
      }
      rm(tmp3)
    }      
  rm(tmp)
  rm(tmp2)
  }
  
  return(results)
}


findrelations<-function(row, results, operations=operations,FE=FE){

  if(nrow(results)==0) {return(row.template)}
  
kmax<-nrow(operations)
imax<-length(results)
  
diffs<-0
print("calculate diffs")
for (i in 1:nrow(results))
  diffs[i] <- row$mz - results$mz[i]

print("calculate diff error margins")

margin<-sqrt(2)*row$mz*FE

diffsminerror<-diffs-margin
diffspluerror<-diffs+margin

ret<-row.template

print("find relations loop")
for(k in 1:nrow(operations))
{
  print(paste("looking for relations:", operations$relation.type[k]))
  for(n in c(-1*operations$multiplier[k]:1, 1:operations$multiplier[k])) {
#    print(paste("multiplicator", n))
    x <- as.numeric(operations$calc.mass[k])* n
    mat<- results$id[which(diffspluerror > x & diffsminerror < x)]
#    print(length(mat))
      for (i in mat)
        {
            nr<-nrow(ret)
            print("hit")
            new.row<-row.template
#            print(new.row)
            new.row[1,c("id", "mz")]<-row[,c("id", "mz")]
#            print(new.row)
            new.row[1,elements]<-
              as.numeric(operations[k,elements])*n+ results[results$id==i, elements]
#            print(new.row)
            new.row[1,"from.id"]<-i
#            print(new.row)
            new.row[1,"multiplier"]<-n
#            print(new.row)
            new.row[1,"relation.type"]<-operations[k, "relation.type"]
#            print(new.row)
            new.row[1,"comment"]<-paste("relations, m/z difference =", diffs[i])
#            print(new.row)
            new.row<-calc.diffs(new.row)
#            print(new.row)
            new.row<-numberize(new.row)
#            print(new.row)
                  ret<-rbind(ret, new.row)
        }
  }
}

return(ret)
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
  if(arg=="minNS") {rows$points<-1/(rows$N+rows$S)}
  if(arg=="minNS.diff") {rows$points<-1000/(rows$N+rows$S)+1/rows$diff.ppm}
  return(rows)  
}
# bf<-data2.bf
# rel<-data2.c13.rel
# clean=T


vanK.plot<-function(df, ylab="H:C", xlab="O:C", elements="CHO", add=F, ...) {
  x <- as.numeric(df$O) / (as.numeric(df$C) + as.numeric(df$C13))
if (elements == "CHO")
  y <- as.numeric(df$H) / (as.numeric(df$C) + as.numeric(df$C13))
if (elements == "CNO")
  y <- as.numeric(df$N) / (as.numeric(df$C) + as.numeric(df$C13))
if (elements == "CNO" & ylab=="H:C")
  ylab<-"N:C"
  if (add==F) {plot(x, y, ylab=ylab, xlab=xlab, ...)} else {points(x, y, ylab=ylab, xlab=xlab, ...)}
}
  
  
dbe<-function(row) {
  return(1+row$C+row$C13+row$N/2-row$H/2)
}

merge.results<-function(bf, rel, arg=F, clean=T) {
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


cleanup<- function(rows) {
#rows<-rows[-is.na(rows$calc.mass),T]
  kill=F
  for (i in which(is.na(rows$calc.mass)!=T)) {
     k1<-min(rows[i,elements])<0
     k2<-rows$C[i]+rows$C13[i]<1
     k3<-rows$H[i]>(2+2*rows$C[i]+2*rows$C13[i]+rows$N[i])
#     k3<-rows$H[i]>(2+2*rows$C[i]+2*rows$C13[i])
     k4<-rows$H[i]>3*rows$C[i]
     k5<-rows$O[i]>rows$C[i]
     k6<-rows$N[i]>rows$C[i]
     if(k1|k2|k3|k4|k5|k6) {kill[i]<-T} else {kill[i]<-F}
  }
  return(rows[kill==F,T])
}



run.all<-function(data, mz) {
  
data<-data[order(data[,mz]),T]
id<-1:nrow(data)
ret<-cbind(id, data )
data.bf<-bforce(ret[,"mz"], ret$id)
data.rel<-  findrelations(ret[,"mz"], ret$id, operations=operations, FE=FE)
data.merge<-merge.results(data.bf, data.rel)

ret<-  merge(ret, data.merge, all.x=T, by="id")

return(ret)

}


test.el<-function(sample, results){
res<-1
for (i in sample$id) {
  if (length(which(results$id==sample$id[i]))==0)
    {res[i]=NA} else {
    if (sum(sample[i, elements]!=results[results$id==sample$id[i], elements])>0) {res[i]<-F} else {res[i]<-T}
  }
}
res2<-data.frame(not.det=sum(is.na(res)), correct.det=sum(res[is.na(res)==F]==T), false.found=sum(res[is.na(res)==F]==F))
return(res2)
}



create.synth<-function(dat,IE=1E-6, n.sel=60, n.synth=1000) {

synth<-row.template
synth [1:n.sel,c("C", "H", "N", "S", "O", "P")]<-
  dat[sample(1:nrow(dat),n.sel, replace=F),c("C", "H", "N", "S", "O", "P")]

for (el in elements)
synth[is.na(synth[,el]),el]<-0
for (i in 1: nrow(synth))
synth[i,T]<-calc.mass(synth[i,T])

operations<-operations[1:8,T]

for(i in 1:(n.synth-nrow(synth))) {
  print(paste(i, "/", (n.synth-n.sel)))
  new.row<-row.template
  new.row<-synth[sample(1:nrow(synth),1),T]
  new.row[,elements]<-as.numeric(operations[sample(1:nrow(operations), 1),elements])*sample(1:5, 1)+as.numeric(new.row[,elements])
  
  new.row<-calc.mass(new.row)
  synth<-rbind(synth, new.row)
}

synth$mz<-synth$calc.mass*(1+rnorm(nrow(synth))*IE)
synth<-synth[order(synth$mz),T]
synth$id<-1:nrow(synth)
return(synth)
}

test.alg<-function(IE=1E-6, n=10, FE=1E-6, n.sel=60, n.synth=1000, verbose=F) {
ret<-data.frame(matrix(ncol=3, nrow=0))
colnames(ret)<-c("not.det", "correct.det", "false.det")

for (i in 1:n) {

  sample<-create.synth(kuja, IE=IE, n.sel=n.sel, n.synth=n.synth)
  results<-bforce(sample, FE=FE, verbose=verbose)
  ret[i,T]<-test.el(sample, results)
  vanK.plot(results, pch=16, cex=0.5)
  vanK.plot(sample, add=T, pch=17, col="red", cex=0.5)
}
return(ret)  
}


