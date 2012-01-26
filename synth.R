setwd("~/Documents/accmass")
source("functions.R")
#sink("synth_echo.txt")


kuja<-read.csv("raw_data/kujawinski2006.csv", sep="\t")

# sample<-create.synth(kuja, IE=0, n.sel=50, n.synth=100)
# bforce(sample)
# data<-
#   sample
# i<-1
# clean<-T
# arg<-F
# results
# i<-2
# row<-data[2,T]
# results
# bruteforce.lim=500
# 
# n<- -1
# k<-1
# diffspluerror
# minerror
# x
# mat
# ret

results.table<-data.frame(matrix(nrow=0,ncol=5))
colnames(results.table)<-c("not.det", "correct.det","false.det","IE", "FE")

for (FE in 1E-6*c(0.1, 0.5, 1.5, 5)){
for (RE in 1E-6*c(20)){
for (IE in 1E-6*c(0, 0.1, 0.5, 1.5, 3)){
  for(i in 1:3){
print(  Sys.time()  )
print(paste("FE =", FE, "IE =", IE, "repeat", i))

res<-test.alg(n.sel=50, n.synth=400, n=1, IE=IE, FE=FE, RE=RE, error.type=RE)

res$IE<-IE
res$FE<-FE
res$RE<-RE
results.table<-rbind(results.table, res)
write.csv(results.table, file="synth_output.csv")
}
}
}
}
#sink()

# outp<-read.csv("synth_output.csv")
# 
# source("~/R/functions.R")
# 
# timeseries(outp$correct.det, outp$IE, outp$FE)
# timeseries(outp$false.det, outp$IE, outp$FE, col="red", add=T)
# 
# # create.synth(kuja, IE=0, n.sel=30, n.synth=100)
# test.alg(n.sel=70, n.synth=1000, n=1, IE=.5E-6, FE=1E-6)
#synth<-create.synth(kuja, n.sel=50, n.synth=400, IE=0)
# for (i in 1:1000) {
#   print(paste(i, "/ 1000"))
#   tmp<-findformula(synth[i,"mz"], c13.max=0)
#   if(nrow(tmp)>0) {results[i,T]<-tmp[1]}
# }
#   