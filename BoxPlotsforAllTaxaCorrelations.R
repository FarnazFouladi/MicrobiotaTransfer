#Farnaz Fouladi
#08-18-2019

#This R code generates Supplementary Figure 5 and Figure 6A. 

rm(list = ls())

library(ggsignif)
library(ggplot2)
library(cowplot)
library(plyr)

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)
myT<-read.table(paste0(output,"SlurryVsMouseCor.txt"),header = TRUE,sep = "\t")
# Attaching taxanomy
taxa<-read.table("/Users/farnazfouladi/CarrollMouseTransfer/taxForwardReads.txt",header = TRUE,sep = "\t",na.strings = "NA")
rownames(taxa)<-sapply(seq(1:nrow(taxa)),function(x){paste0("SV",x)})
taxa<-cbind(rownames(taxa),taxa)
myT1<-merge(myT,taxa,by.x="bugnames",by.y="rownames(taxa)",sort=FALSE)

#Functions for plotting
tukeyResult<-function(dataframe,taxa,time){
  dataframe1<-dataframe[dataframe$time==time,]
  result<-aov(lm(rho~dataframe1[,which(colnames(dataframe1)==taxa)],dataframe1))
  return(sort(TukeyHSD(result)$`dataframe1[, which(colnames(dataframe1) == taxa)]`[,"p adj"]))
}

myPlot<-function(dataframe,taxa,time){
  dataframe1<-dataframe[dataframe$time==time,]
  ggplot(data=dataframe1,aes(x=dataframe1[,which(colnames(dataframe1)==taxa)],y=rho))+
    geom_boxplot(outlier.size=0.5,aes(fill=factor(dataframe1[,which(colnames(dataframe1)==taxa)])))+
    theme(axis.text.x = element_blank())+labs(y="Rho",x="")+
    labs(title=paste0("Week= ",time))+theme(legend.position = "none")
  
}


theme_set(theme_classic(base_size = 8))
myList<-list()
index<-1

#Week1, Phylum
tukeyResult(myT1,"Phylum",1)
plot<-myPlot(myT1,"Phylum",1)
plot1<-plot+geom_signif(y_position=c(1), xmin=c(2), xmax=c(3),annotation=c("***"), tip_length=0,vjust = 0.5,textsize =2)+theme(legend.position = "none")
myList[[index]]<-plot1
index<-index+1
#Week1, Class
tukeyResult(myT1,"Class",1)
plot<-myPlot(myT1,"Class",1)
plot1<-plot+geom_signif(y_position=c(1,1.05), xmin=c(3,3), xmax=c(5,8),annotation=c("***","**"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week1, Order
tukeyResult(myT1,"Order",1)
plot<-myPlot(myT1,"Order",1)
plot1<-plot+geom_signif(y_position=c(1,1.05), xmin=c(1,1), xmax=c(4,8),annotation=c("***","**"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week1, Family
tukeyResult(myT1,"Family",1)
plot<-myPlot(myT1,"Family",1)
plot1<-plot+geom_signif(y_position=c(1,1.05,1.1,1.15), xmin=c(3,3,14,14), xmax=c(14,21,18,19),annotation=c("***","**","**","*"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week2, Phylum
tukeyResult(myT1,"Phylum",2)
plot<-myPlot(myT1,"Phylum",2)
plot1<-plot+geom_signif(y_position=c(1,1), xmin=c(1,2.1), xmax=c(1.9,3),annotation=c("*","***"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week2, Class
tukeyResult(myT1,"Class",2)
plot<-myPlot(myT1,"Class",2)
plot1<-plot+geom_signif(y_position=c(1,1.05,1.1), xmin=c(2,3,3), xmax=c(3,5,8),annotation=c("*","***","***"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week2, Order
tukeyResult(myT1,"Order",2)
plot<-myPlot(myT1,"Order",2)
plot1<-plot+geom_signif(y_position=c(1,1.05,1.1), xmin=c(1,1,1), xmax=c(4,8,9),annotation=c("***","***","*"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week2, Family
tukeyResult(myT1,"Family",2)
plot<-myPlot(myT1,"Family",2)
plot1<-plot+geom_signif(y_position=c(1,1.05,1.1,1.15,1.2,1.25), xmin=c(3,3,3,14,14,19), xmax=c(11,14,21,18,19,21),annotation=c("*","***","**","**","**","*"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week3, Phylum
tukeyResult(myT1,"Phylum",3)
plot<-myPlot(myT1,"Phylum",3)
plot1<-plot+geom_signif(y_position=c(1,1), xmin=c(1,2.1), xmax=c(1.9,3),annotation=c("**","***"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week3, Class
tukeyResult(myT1,"Class",3)
plot<-myPlot(myT1,"Class",3)
plot1<-plot+geom_signif(y_position=c(1.05,1.1), xmin=c(3,3), xmax=c(5,8),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week3, Order
tukeyResult(myT1,"Order",3)
plot<-myPlot(myT1,"Order",3)
plot1<-plot+geom_signif(y_position=c(1,1.05), xmin=c(1,1), xmax=c(4,8),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week3, Family
tukeyResult(myT1,"Family",3)
plot<-myPlot(myT1,"Family",3)
plot1<-plot+geom_signif(y_position=c(1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4), xmin=c(3,3,3,11,11,14,14,18,19), xmax=c(11,14,21,18,19,18,19,21,21),annotation=c("**","***","**","*","*","***","**","**","*"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week4, Phylum
tukeyResult(myT1,"Phylum",4)
plot<-myPlot(myT1,"Phylum",4)
plot1<-plot+geom_signif(y_position=c(1,1), xmin=c(1,2.1), xmax=c(1.9,3),annotation=c("**","***"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week4, Class
tukeyResult(myT1,"Class",4)
plot<-myPlot(myT1,"Class",4)
plot1<-plot+geom_signif(y_position=c(1.05,1.1), xmin=c(3,3), xmax=c(5,8),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1

#Week4, Order
tukeyResult(myT1,"Order",4)
plot<-myPlot(myT1,"Order",4)
plot1<-plot+geom_signif(y_position=c(1,1.05), xmin=c(1,1), xmax=c(4,8),annotation=c("***","***"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1
#Week4, Family
tukeyResult(myT1,"Family",4)
plot<-myPlot(myT1,"Family",4)
plot1<-plot+geom_signif(y_position=c(1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4), xmin=c(3,3,3,11,11,14,14,19,18), xmax=c(11,14,21,18,19,18,19,21,21),annotation=c("**","***","***","*","*","**","**","*","*"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1

#Find average abd sd for rho
data4<-myT1[myT1$time==4,]
ddply(data4,.(Class),colwise(mean))[,c(1,4)]
ddply(data4,.(Class),colwise(sd))[,c(1,4)]
ddply(data4,.(Family),colwise(mean))[,c(1,4)]
ddply(data4,.(Family),colwise(sd))[,c(1,4)]

png("boxPlotsforALLtaxaCorrelationPanel.png", units="in", width=10, height=10,res=300)
plot_grid(myList[[1]],myList[[5]],myList[[9]],myList[[13]],
          myList[[2]],myList[[6]],myList[[10]],myList[[14]],
          myList[[3]],myList[[7]],myList[[11]],myList[[15]],
          myList[[4]],myList[[8]],myList[[12]],myList[[16]],nrow=4,ncol=4)

dev.off()
