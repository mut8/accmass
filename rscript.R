setwd("~/Documents/accmass")
source("functions.R")

files<-c(
"10G_Laubextrakt.csv",
"11G_Laubextrakt.csv",
"12G_Laubextrakt.csv",
"1BFE0_Laubextrakt.csv",
"2BFE1_Laubextrakt.csv",
"3BFE_Laubextrakt.csv",
"4B_Laubextrakt.csv",
"5B_Laubextrakt.csv",
"6B_Laubextrakt.csv",
"7GFE0_1_Laubextrakt.csv",
"8GFE_Laubextrakt.csv",
"9GFE_Laubextrakt.csv",
"W29112010_Laubextrakt.csv",
"WPETACMEOH_G_Laubextrakt.csv",
"WPETACMEOH_MQ_Laubextrakt.csv"
)


dir<-"data_litter/"
for (i in 1:length(files)) {
  data1<-read.csv(paste(dir, files[i], sep=""))
  data1<-data1[8:nrow(data1), 1:3]
  colnames(data1)<-c("mz", "Intensity", "rel.int")
  sample<-rep(files[i], nrow(data1))
  data1<-cbind(data1, sample)
  if (i == 1) {data.merge<-data1} else {data.merge<-rbind(data.merge, data1)}
}

data.merge2<-data.merge
data.merge2$mz<-as.numeric(as.character(data.merge$mz))
data.merge2$Intensity<-as.numeric(as.character(data.merge$Intensity))
data.merge2$rel.int<-as.numeric(as.character(data.merge$rel.int))

data1[8:, T]
data.red<-data.merge2[data.merge2$rel.int>1, T]

data.red<-data.red[order(data.red$mz),T]

dist<-dist(data.red$mz)

mz.min<-min(data.red$mz)
mz.max<-max(data.red$mz)

(mz.max-mz.min)/20

data.red2<-data.red

data.red2$groups<-0
groups.max<-0
for (i in 80+(0:40)*40) {
  print(paste(i, i+40, nrow(data.red2[data.red2$mz>i & data.red2$mz < i + 40,T])))
dist<-dist(data.red2$mz[data.red2$mz>i & data.red2$mz < i + 40])  
clust<-hclust(dist)
data.red2$groups[data.red2$mz>i & data.red2$mz < i + 40]<-cutree(clust, h=(i+20)*3E-6)+groups.max
groups.max<-max(data.red2$groups)
}

dist<-dist(data.red2$mz[data.red2$mz>i])  
clust<-hclust(dist)
data.red2$groups[data.red2$mz>i]<-cutree(clust, h=(i+20)*1E-6)+groups.max
groups2.max<-max(data.red2$groups)

data.3a<-data.frame(unique(data.red2$groups))
data.4a<-cbind(data.3a,  tapply(data.red2$mz, data.red2$groups, mean), tapply(data.red2$mz, data.red2$groups, sd))

colnames(data.4a)<-c("id", "mz", "mz.sd")

head(data.4a)
$data.3a
nrow(data.4)

results.data.4a.3ppm<-results.data.4a
results.data.4a.3ppm.minNS<-bforce(H.substract(data.4a), arg="minNS")

results.data.4a$id<-ceiling(results.data.4a$id)+1
results.data.4a.3ppm.minNS$id<-ceiling(results.data.4a.3ppm.minNS$id)


merge()

pdf("vank2_minppm_3ppm.pdf", height=10, width=6)

par(mfrow=c(2,1), tck=0.01)

for(i in 1:length(files)) {
results.file.merge<-
  merge(data.red2[data.red2$sample==files[i],T], results.data.4a.3ppm, all.x=T, by.x="groups", by.y="id")
cex=(log(results.file.merge$Intensity)-min(log(results.file.merge$Intensity)))/(max(log(results.file.merge$Intensity))-min(log(results.file.merge$Intensity)))

par(mar=c(0,4.1,4.1,2.1))
vanK.plot(results.file.merge, cex=cex, main=files[i], elements="CNO", xaxt="n")
axis(1, labels=F)
axis(3, labels=F)
axis(4, labels=F)
par(mar=c(4.1,4.1,0,2.1))
vanK.plot(results.file.merge, cex=cex*0.7)  
axis(3, labels=F)
axis(4, labels=F)
}
dev.off()


pdf("vank2_minNS_3ppm.pdf", height=10, width=6)

par(mfrow=c(2,1), tck=0.01)

for(i in 1:length(files)) {
results.file.merge<-
  merge(data.red2[data.red2$sample==files[i],T], 
results.data.4a.3ppm.minNS, all.x=T, by.x="groups", by.y="id")
cex=(log(results.file.merge$Intensity)-min(log(results.file.merge$Intensity)))/(max(log(results.file.merge$Intensity))-min(log(results.file.merge$Intensity)))

par(mar=c(0,4.1,4.1,2.1))
vanK.plot(results.file.merge, cex=cex, main=files[i], elements="CNO", xaxt="n")
axis(1, labels=F)
axis(3, labels=F)
axis(4, labels=F)
par(mar=c(4.1,4.1,0,2.1))
vanK.plot(results.file.merge, cex=cex*0.7)  
axis(3, labels=F)
axis(4, labels=F)
}
dev.off()


results.file.merge
i<-1

files

i<-1
hist(cex)


, bruteforce.lim=400)

, tapply(data.red$groups, data.red$groups, length))

data.red$mz[data.red$mz>i]

length(unique(data.red$groups))

groups.max

data.red$mz[20000]

90+(0:80)*20
  
10:1600, step=20  
  
i
<-90

bforce(uja)
data2<-data2[order(data2$mz),T]
data2$mz<-H.substract(data2$mz)

outp<-bforce(data2)

data2$km<-data2$mz*14/(mass.H1*2+mass.C12)

data2$kmd<-ceiling(data2$km)-data2$km
data2$z<-as.integer(data2$mz) %% 14 - 14

plot(data2$kmd, data2$z, cex=0.2)

data2$id<-1:nrow(data2)
head(data2)
files<-list.files("/home/lkohl/Documents/accmass/raw_data/Litter")
files<-files[grep(".csv", files)]
files[i]
#
read.csv("data_litter/3BFE.csv")


data2
data.merge<-merge(data2, results, by.x="id", by.y="id", all.x=T)
vanK.plot(data.merge, cex=cex)
max(cex)

cex=1+log10(data.merge$Intensity/sum(data.merge$Intensity))/5

head(data2)


files<-files[1]

file<-1
for (i in 1:length(files))
{
  print(paste(files[i], ": file", i, "/", length(files)))
data<-read.csv(paste("raw_data/SoilWater/", files[i], sep=""))
data2<-data.frame(data[-(1:7), 1:7])
colnames(data2)<-c("mz", "Intensity", "Relative Intensity", "Resolution", "Charge", "Baselin", "Noise")
length(file)<-nrow(data2)
file<-rep(data[1,1], nrow(data2))
data3<-data.frame(data2, file)
if (i == 1) {
  data.collected<-data3
  } else { data.collected<-rbind(data.collected, data3)}

}

data.collected2<-data.collected
data.collected2$Noise<-as.numeric(data.collected$Noise)
data.collected2$Intensity<-as.numeric(data.collected$Intensity)

data.reduced<-
  data.collected2[data.collected2$Intensity>3*data.collected2$Noise, T]
data.reduced$mz<-as.numeric(as.character(data.reduced$mz))

data.reduced<-data.reduced[order(data.reduced$mz),T]

data2<-read.csv("raw_data/test2.csv")
data2<-data2[data2$Intensity > 3* data2$Noise,T ]

kuja<-read.csv("raw_data/kujawinski2006.csv")    
  
cluster.R	
	
rscr

head(data2$mz)


head(
  data2
  )

data2<-data2[order(data2$mz),T]

data.reduced$mz

dist<-dist(log(data.reduced$mz))
clust<-hclust(dist)

data.reduced$mz<-H.substract(data.reduced$mz)
Sys.time()
results<-run.all(data.reduced, "mz")
Sys.time()

write.csv("results.csv", results)
sink()
