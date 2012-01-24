setwd("~/Documents/accmass")
source("functions.R")
sink("synth_echo.txt")

kuja<-read.csv("raw_data/kujawinski2006.csv", sep="\t")


results.table<-data.frame(matrix(nrow=0,ncol=5))
colnames(results.table)<-c("not.det", "correct.det","false.det","IE", "FE")

for (FE in 1E-6*c(0.1, 0.3, 0.5, 1, 1.5, 2.5, 5)){
for (IE in 1E-6*c(0.1, 0.3, 0.5, 1, 1.5, 2.5, 5)){
  for(i in 1:10){
print(  Sys.time()  )
print(paste("FE =", FE, "IE =", IE, "repeat", i))

res<-test.alg(n.sel=70, n.synth=1000, n=1, IE=0, FE=3)

res$IE<-IE
res$FE<-FE
results.table<-rbind(results.table, res)
write.csv(results.table, file="synth_output.csv")
}
}
}
sink()

# outp<-read.csv("synth_output.csv")
# 
# source("~/R/functions.R")
# 
# timeseries(outp$correct.det, outp$IE, outp$FE)
# timeseries(outp$false.det, outp$IE, outp$FE, col="red", add=T)
# 
# create.synth(kuja, IE=0, n.sel=30, n.synth=100)