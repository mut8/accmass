#MS error margins
error=3E-6
FE=error
RE=2E-5
C13E=1E-4

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
    
                  #check for N-rule (n+h %% 2 == 0)
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
  #colnames(results)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm")
  #colnames(results)<-
  #  colnames(row.template)
  #if(nrow(results)>0)
  results<-results[order(abs(as.numeric(results$difference))),T]
  return(results)
}

findrelation<-function(row, res, mult=1:5, c=0, h=0, o=0, n=0, s=0, p=0, Na=0, c13=0, RE=2E-5) {
  mz<-res[[1]]
  lst<-res[[2]]
  rel.mass= c*mass.C12+c13*mass.C13+h*mass.H1+o*mass.O16+n*mass.N14+s*mass.S32+p*mass.P31+Na*mass.Na23
  
  for (i in mult) {
  target.mass=row$mass+rel.mass*i
  rel.error=rel.mass*i*RE

  new.row<-row
  new.row$C[1]<-row$C[1]+c*i
  new.row$H[1]<-row$H[1]+h*i
  new.row$O[1]<-row$O[1]+o*i
  new.row$N[1]<-row$N[1]+n*i
  new.row$S[1]<-row$S[1]+s*i
  new.row$P[1]<-row$P[1]+p*i
  new.row$Na[1]<-row$Na[1]+Na*i
  
  new.row<-calc.mass(new.row)
  
  hits<-which(mz>(target.mass-rel.error)&mz<(target.mass+rel.error))
  print(hits)
  
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

bforce<-function(mz, FE=3E-6, RE=2E-5, bruteforce.lim=500, verbose=T) {

#  #create result dataframe, set colnames
#  data<-data.frame(matrix(nrow=length(mz), ncol=12))
#  colnames(data)<-c("mz","C","C13","H","N","O","S","P","Na","calc.mass", "difference", "diff.ppm")
  
  #add ascending m/z values to results dataframe
  mz<-sort(mz)
  
  #create result list, set colnames
  tmp<-data.frame(matrix(nrow=0, ncol=12))#CH4-O
  print("+CH4 -O")
  res<-findrelation(row, res, c=1, h=4, o=-1)

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

  return(list(mz, lst))
}

plotresults<-function(res) {
  mz<-res[[1]]
  lst<-res[[2]]
  nr<-1
  for (i in 1:length(mz)) {
    nr[i]<-nrow(lst[[i]])    
  }
  plot(mz, nr, cex=.5, ylim=c(0,20), ylab="hits", xlab="mass nr")
}

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


diffs.matrix<-function(inp){
diffs<-matrix(nrow=length(inp), ncol=length(inp))
relations<-matrix(nrow=length(inp), ncol=length(inp))
for (i in 1:length(inp))
for (j in 1:length(inp)) {
  diffs[i,j]<-inp[i]-inp[j]
}
return(diffs)
}
