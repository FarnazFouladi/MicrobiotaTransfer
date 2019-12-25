#Farnaz Fouladi
#08-18-2019

#This R code compares colonization efficiency between taxa.

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

#Extract Legend
plot<-ggplot(data=myT1,aes(x=Phylum,y=rho))+geom_boxplot(outlier.size=0.5,aes(fill=factor(Phylum,labels=c("Actinobacteria (12)","Bacteroidetes (39)","Firmicutes (211)","Proteobacteria (7)","Verrucomicrobia (2)"))))+
  theme(axis.text.x = element_blank(),legend.text = element_text(size=7),legend.title = element_text(size = 8),legend.key.size = unit(0.7,"line"))+labs(y="rho",x="")+
  labs(title="",fill="Phylum")

pdf("legendPhylumRho.pdf",height = 1,width = 2)
legendp <- cowplot::get_legend(plot)
grid.newpage()
grid.draw(legendp)
dev.off()

#Week1, Class
tukeyResult(myT1,"Class",1)
plot<-myPlot(myT1,"Class",1)
plot1<-plot+geom_signif(y_position=c(1,1.05), xmin=c(3,3), xmax=c(5,8),annotation=c("***","**"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1

#ExtractLegend
plot<-ggplot(data=myT1,aes(x=Class,y=rho))+geom_boxplot(outlier.size=0.5,aes(fill=factor(Class,labels=c("Actinobacteria (1)","Bacilli(3)","Bacteroidia (39)","Betaproteobacteria (4)","Clostridia (188)","Coriobacteriia (11)","Deltaproteobacteria (2)","Erysipelotrichia (18)","Gammaproteobacteria (1)","Negativicutes (2)","Verrucomicrobiae (2)"))))+
  theme(axis.text.x = element_blank(),legend.text = element_text(size=7),legend.title = element_text(size = 8),legend.key.size = unit(0.7,"line"))+labs(y="rho",x="")+
  labs(title="",fill="Class")

pdf("legendClassRho.pdf",height = 2,width = 3)
legendp <- cowplot::get_legend(plot)
grid.newpage()
grid.draw(legendp)
dev.off()

#Week1, Order
tukeyResult(myT1,"Order",1)
plot<-myPlot(myT1,"Order",1)
plot1<-plot+geom_signif(y_position=c(1,1.05), xmin=c(1,1), xmax=c(4,8),annotation=c("***","**"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1

#Extract legend
plot<-ggplot(data=myT1,aes(x=Order,y=rho))+geom_boxplot(outlier.size=0.5,aes(fill=factor(Order,labels=c("Bacteroidales (39)","Bifidobacteriales (1)","Burkholderiales (4)","Clostridiales (187)","Coriobacteriales (11)","Desulfovibrionales (2)","Enterobacteriales (1)","Erysipelotrichales (18)","Lactobacillales (3)","Selenomonadales (2)","Verrucomicrobiales (2)"))))+
  theme(axis.text.x = element_blank(),legend.text = element_text(size=7),legend.title = element_text(size = 8),legend.key.size = unit(0.7,"line"))+labs(y="rho",x="")+
  labs(title="",fill="Order")

pdf("legendOrderRho.pdf",height = 2,width = 3)
legendp <- cowplot::get_legend(plot)
grid.newpage()
grid.draw(legendp)
dev.off()

#Week1, Family
tukeyResult(myT1,"Family",1)
plot<-myPlot(myT1,"Family",1)
plot1<-plot+geom_signif(y_position=c(1,1.05,1.1,1.15), xmin=c(3,3,14,14), xmax=c(14,21,18,19),annotation=c("***","**","**","*"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot1
index<-index+1

#Extract legend
plot<-ggplot(data=myT1,aes(x=Family,y=rho))+geom_boxplot(outlier.size=0.5,aes(fill=factor(Family,labels=c("Acidaminococcaceae (1)","Alcaligenaceae (3)","Bacteroidaceae (17)","Bifidobacteriaceae (1)","Christensenellaceae (5)","Coriobacteriaceae (11)","Defluviitaleaceae (2)","Desulfovibrionaceae (2)","Enterobacteriaceae (1)","Enterococcaceae (1)","Erysipelotrichaceae (18)","Eubacteriaceae  (3)","Family_XIII (10)","Lachnospiraceae (95)","Lactobacillaceae (1)","Oxalobacteraceae (1)","Peptostreptococcaceae (2)","Porphyromonadaceae (8)","Prevotellaceae (6)","Rikenellaceae (8)","Ruminococcaceae (67)","Streptococcaceae (1)","Veillonellaceae (1)","Verrucomicrobiaceae (2)"))))+
  theme(axis.text.x = element_blank(),legend.text = element_text(size=7) ,legend.title = element_text(size = 8),legend.key.size = unit(0.7,"line"))+
  labs(title="",fill="Family")

pdf("legendFamilyRho.pdf",height = 2,width = 4)
legendp <- cowplot::get_legend(plot)
grid.newpage()
grid.draw(legendp)
dev.off()

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

#Find average and sd for rho
data4<-myT1[myT1$time==4,]
ddply(data4,.(Class),colwise(mean))[,c(1,4)]
ddply(data4,.(Class),colwise(sd))[,c(1,4)]
ddply(data4,.(Family),colwise(mean))[,c(1,4)]
ddply(data4,.(Family),colwise(sd))[,c(1,4)]

pdf("boxPlotsforALLtaxaCorrelationPanel.pdf",height = 10,width = 10)
plot_grid(myList[[1]],myList[[5]],myList[[9]],myList[[13]],
          myList[[2]],myList[[6]],myList[[10]],myList[[14]],
          myList[[3]],myList[[7]],myList[[11]],myList[[15]],
          myList[[4]],myList[[8]],myList[[12]],myList[[16]],nrow=4,ncol=4)

dev.off()
