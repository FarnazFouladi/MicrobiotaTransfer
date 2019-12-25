#Farnaz Fouladi
#08-18-2019

#This R code compares Firmicutes:Bacteroidetes ratio between human fecal pellets,
#slurries, and mouse fecal pellets

rm(list = ls())

library(ggplot2)
library(ggsignif)
library(plyr)

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)

data<-read.table(paste0(input,"PhylumTableDada2MetaAdded.txt"),sep="\t",header = TRUE,comment.char = "")
data1<-data[,2:18]
#Table is log-normalized
data1<-(10^data1)-1
data1$ratio<-log10(data1$Firmicutes/data1$Bacteroidetes)
data2<-cbind(data1,data[,19:ncol(data)])
#Comparing ratio between human donors, slurry, and mice:
data2<-data2[data2$Sample.type=="Human.donor"|data2$Sample.type=="Fecal.slurry"|data2$Sample.type=="Mouse.feces",]
data2$variable2<-factor(data2$Sample.type,levels = c("Human.donor","Fecal.slurry","Mouse.feces"))


pdf("FirmicutestoBacteroidetes.pdf",width=5, height=5)
theme_set(theme_gray(base_size = 12))
ggplot(data=data2,aes(x=variable2,,y=ratio))+geom_boxplot()+labs(x="",y="log10(Firmicutes:Bacteroidetes)")+
  geom_signif(y_position=c(3.1,3.3,3.5), xmin=c(1,1,2), xmax=c(2,3,3),annotation=c("***","***","***"), tip_length=0,vjust = 0.5,textsize =4)+
  theme(axis.title.x = element_blank())+scale_x_discrete(labels=c("Human fecal samples","Slurries","Mouse fecal pellets"))

dev.off()

result<-aov(lm(ratio~as.factor(variable2),data=data2))
TukeyHSD(result)

ddply(data2,.(variable2),colwise(mean))[,c("variable2","ratio")]
ddply(data2,.(variable2),colwise(sd))[,c("variable2","ratio")]
