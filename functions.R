#MS error margins


error=3E-6
FE=error
RE=2E-5
C13E=1E-4

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


H.substract<-function(x) return(x-mass.H1+mass.e)

row.template<-data.frame(matrix(ncol=13, nrow=0))
colnames(row.template)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm", "comment")


calc.mass<-function(data) {
  data[,c("C", "N", "C13", "O","H", "S",  "P", "Na")]<-as.numeric(data[,c("C", "N", "C13", "O","H", "S",  "P", "Na")])
  data$calc.mass<-as.numeric(data$C*mass.C12+data$N*mass.N14+data$O*mass.O16+data$H*mass.H1+data$C13*mass.C13+data$Na*mass.Na23+data$P*mass.P31+data$S*mass.S32)
  return(data)
}

##
##findformula: algorithm to calculate sum formula from accurate mass
##


findformula<-function(x, c.h.min=1/3, c.o.max=1, c13.max=2, c.n.max=1, s.max=3, p.max=3, Na.max=0, screen.error=3E-2, FE=3E-6) {
  
  results<-row.template
  
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
                    if(calc.mass>x-(x*FE) & calc.mass<x+(x*FE))
                    {
                      #print result
                      #print(paste("C", c, "C13", c13, "H", h, "N", n, "O", o,"P",p,"S",s,"Na",Na, "m=", calc.mass, "difference", calc.mass-x, "diff.ppm", 1E6*(calc.mass-x)/x))
                      tmp<-nrow(results)
                      results[tmp+1, T]<-c(x, c, c13, h, n, o, s, p, Na, calc.mass, calc.mass-x, 1E6*(calc.mass-x)/x, "bruteforce")
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

  results<-results[order(abs(as.numeric(results$difference))),T]
  return(results)
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
    if(nrow(tmp)>0) {
    #add all results found to the results list
    tmp<-rbind(tmp, lst[[i]])
    #order results and add them to lst
    lst[[i]]<-tmp[order(abs(as.numeric(tmp$diff.ppm))),T]
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

  return(lst)
}

plotresults<-function(lst, mz, xlab="mass nr", ylab="hits", cex=.5, ...) {
  nr<-1
  for (i in 1:length(mz)) {
    nr[i]<-nrow(lst[[i]])    
  }
  plot(mz, nr, cex=cex, ylab=ylab, xlab=xlab, ...)
}



findrelations<-function(inp){

kmax<-nrow(operations)
imax<-length(inp)
jmax<-imax
nmax=5
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

for (i in 1:imax) {
print(paste("calculate diff error margins", i, "/", imax))
  for (j in 1:jmax) {

diffsminerror[i,j]<-diffs[i,j]-FE*1.414214*max(inp[c(i,j)])
diffspluerror[i,j]<-diffs[i,j]+FE*1.414214*max(inp[c(i,j)])
}}

relationlist<-list(row.template)
for (i in 1:imax) {
  relationlist[[i]]<-row.template
  colnames(relationlist[[i]])<-c("mz.from", "C", "C13", "H", "N", "O", "S", "P", "Na", "n", "difference", "difference.ppm", "comment")
}


print("find relations loop")
for(k in 1:nrow(operations))
{
  print(paste("looking for relations:", operations$comment[k]))
  for(n in c(-nmax:-1, 1:nmax)) 
  {
    print(paste("multiplicator", n))
    x<- as.numeric(operations$calc.mass[k])* n
    mat<-diffspluerror > x & diffsminerror < x 
      for (i in 2:(imax))
        for (j in 1:(i-1))
          if (mat[i,j])
          {
            nr<-nrow(relationlist[[i]])
            relationlist[[i]][(nr+1),c("C","C13","H","N","O","S","P","Na")]<-as.numeric(operations[k,c("C","C13","H","N","O","S","P","Na")])*n
            relationlist[[i]][(nr+1),"mz.from"]<-inp[j]
            relationlist[[i]][(nr+1),"n"]<-n
            relationlist[[i]][(nr+1),"comment"]<-paste(operations[k, "comment"], diffs[i,j])
            nr<-nrow(relationlist[[j]])
            relationlist[[j]][(nr+1),c("C","C13","H","N","O","S","P","Na")]<-as.numeric(operations[k,c("C","C13","H","N","O","S","P","Na")])*n-1
            relationlist[[j]][(nr+1),"mz.from"]<-inp[i]
            relationlist[[j]][(nr+1),"n"]<-n
            relationlist[[j]][(nr+1),"comment"]<-paste(operations[k,"comment"], -1*diffs[i,j])
            hits<-hits+1
          }
  }
}

return(relationlist)
}

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
  
  
order.rows<-function(df, order="difference") {
  if(order=="difference"|order=="dif"|order=="diff"|order=="d")
  res<-df[order(abs(as.numeric(df$diff.ppm))),T]
  #if(order=="NPSNa")
  #res<-df[order(rowSums(as.numeric(df[,c("N", "P", "S", "Na")])))]
  return(res)
}

dbe<-function(row) {
  return(1+row$C+row$C13+row$N/2-row$H/2)
}


collect<-function(lst){
  mat<-row.template
  for(i in 1:length(lst))
  mat[i,T]<-lst[[i]][1,T]
  return(mat)
}
 
