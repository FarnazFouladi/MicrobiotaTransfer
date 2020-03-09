#Farnaz Fouladi
#08-18-2019

#This R code compares correlations between SVs and phenotypes between taxa. 

rm(list = ls())

library(ggsignif)
library(ggplot2)
library(cowplot)

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)

myT<-read.table(paste0(output,"interactionPhenotypeBug.text"),header = TRUE,sep = "\t")
# Attaching taxanomy
taxa<-read.table("/Users/farnazfouladi/CarrollMouseTransfer/taxForwardReads.txt",header = TRUE,sep = "\t",na.strings = "NA")
rownames(taxa)<-sapply(seq(1:nrow(taxa)),function(x){paste0("SV",x)})
taxa<-cbind(rownames(taxa),taxa)
myT1<-merge(myT,taxa,by.x="bugName",by.y="rownames(taxa)",sort=FALSE)

#Change the name of variables
myT1$variableName<-sapply(as.character(myT1$variableName),function(x){
  if (x=="Body.wt") return("Body weight")
  else if (x=="Body.wt.pct.change") return("Body weight percentage change")
  else if (x=="Brown.fat.wt") return("Brown fat weight")
  else if (x=="Cecum.wt") return("Cecum weight")
  else if (x=="Cum.food.consumption") return("Cumulative food intake")
  else if (x=="Fat.mass") return("Fat mass")
  else if (x=="Fat.mass.pct.change") return("Fat mass percentage change")
  else if (x=="Gonadal.fat.wt") return("Gonadal fat weight")
  else if (x=="Rel.brown.fat.wt") return("Relative brown fat weight")
  else if (x=="Rel.cecum.wt") return("Relative cecum weight")
  else if (x=="Rel.gonadal.fat.wt") return("Relative gonadal fat weight")
  else if (x=="Rel.SI.wt") return("Relative small intestinal weight")
  else if (x=="SI.wt") return("Small intestinal weight")
  else if (x=="Wk.food.consumption") return("Weekly food intake")
  else return(print("error"))
})

# Plots for taxa that were significantly associate with a phenotype and also there was a significant difference between taxa
#*** p<0.001
#** p<0.01
#* p<0.05

#Functions for plotting

tukeyResult<-function(dataframe,taxa,time,variable){
  dataframe1<-dataframe[dataframe$time==time & dataframe$variableName==variable,]
  result<-aov(lm(log10(adjustedSpearman)~dataframe1[,which(colnames(dataframe1)==taxa)],dataframe1))
  return(sort(TukeyHSD(result)$`dataframe1[, which(colnames(dataframe1) == taxa)]`[,"p adj"]))
}


myPlot<-function(dataframe,taxa,time,variable){
  dataframe1<-dataframe[dataframe$time==time & dataframe$variableName==variable,]
  ggplot(data=dataframe1,aes(x=dataframe1[,which(colnames(dataframe1)==taxa)],y=-log10(adjustedSpearman)))+
    geom_boxplot(outlier.size=0.5,aes(fill=factor(dataframe1[,which(colnames(dataframe1)==taxa)])))+
    theme(axis.text.x = element_blank())+labs(y="-log10(adjusted p-value)",x="")+
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black", size=1)+
    labs(title=paste0(variable,"\nWeek= ",time))+theme(legend.position = "none")
  
}

myList<-list()
index<-1
theme_set(theme_classic(base_size = 9.5))

###Phylum
tukeyResult(myT1,"Phylum",time=2,"Fat mass percentage change")
plot<-myPlot(myT1,"Phylum",time=2,"Fat mass percentage change")
  
myList[[index]]<-plot
index<-index+1


#Extract Legend
plot<-ggplot(data=myT1,aes(x=Phylum,y=-log10(adjustedSpearman)))+geom_boxplot(outlier.size=0.5,aes(fill=factor(Phylum,labels=c("Actinobacteria (12)","Bacteroidetes(39)","Firmicutes (198)","Proteobacteria (7)","Verrucomicrobia (2)"))))+
  theme(axis.text.x = element_blank(),legend.text = element_text(size=7),legend.title = element_text(size = 8),
        legend.key.size = unit(0.7,"line"))+geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black", size=1)+labs(y="-log10(adjusted p-value)",x="")+
  labs(title="Fat mass percentage change\nWeek=2",fill="Phylum")+geom_signif(y_position=c(3.5, 3.5), xmin=c(2), xmax=c(3),annotation=c("*"), tip_length=0,vjust = 0.5,textsize =4)
pdf("legendPhylumIntercation.pdf",height = 1,width = 2)
legendp <- cowplot::get_legend(plot)
grid.newpage()
grid.draw(legendp)
dev.off()

tukeyResult(myT1,"Phylum",time=3,"Fat mass percentage change")
plot<-myPlot(myT1,"Phylum",time=3,"Fat mass percentage change")+
  geom_signif(y_position=c(1.5), xmin=c(2), xmax=c(3),annotation=c("**"), tip_length=0,vjust = 0.5,textsize =4)
myList[[index]]<-plot
index<-index+1

tukeyResult(myT1,"Phylum",time=4,"Fat mass percentage change")
plot<-myPlot(myT1,"Phylum",time=4,"Fat mass percentage change")+
  geom_signif(y_position=c(1.5,1.55), xmin=c(2,2), xmax=c(3,4),annotation=c("***", "**"), tip_length=0,vjust = 0.5,textsize =4)
myList[[index]]<-plot
index<-index+1

###Class
tukeyResult(myT1,"Class",time=3,"Fat mass percentage change")
plot<-myPlot(myT1,"Class",time=3,"Fat mass percentage change")+
  geom_signif(y_position=c(1.5), xmin=c(3), xmax=c(5),annotation=c("*"), tip_length=0,vjust = 0.5,textsize =4)
myList[[index]]<-plot
index<-index+1

#Extract Legend

plot<-ggplot(data=myT1,aes(x=Class,y=-log10(adjustedSpearman)))+geom_boxplot(outlier.size=0.5,aes(fill=factor(Class,labels=c("Actinobacteria (1)","Bacilli(3)","Bacteroidia (39)","Betaproteobacteria (4)","Clostridia (177)","Coriobacteriia (11)","Deltaproteobacteria (2)","Erysipelotrichia (17)","Gammaproteobacteria (1)","Negativicutes (1)","Verrucomicrobiae (2)"))))+
  theme(axis.text.x = element_blank(),legend.text = element_text(size=7),legend.title = element_text(size = 8),
        legend.key.size = unit(0.7,"line"))+geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black", size=1)+labs(y="-log10(adjusted p-value)",x="")+
  labs(title="Fat mass percentage change\nWeek=3",fill="Class")+geom_signif(y_position=c(3), xmin=c(3), xmax=c(5),annotation=c("**"), tip_length=0,vjust = 0.5,textsize =4)
pdf("legendClassInteraction.pdf",height = 2,width = 3)
legendp <- cowplot::get_legend(plot)
grid.newpage()
grid.draw(legendp)
dev.off()

tukeyResult(myT1,"Class",time=4,"Fat mass percentage change")
plot<-myPlot(myT1,"Class",time=4,"Fat mass percentage change")+
  geom_signif(y_position=c(1.5), xmin=c(3), xmax=c(5),annotation=c("***"), tip_length=0,vjust = 0.5,textsize =4)
myList[[index]]<-plot
index<-index+1

### Order

tukeyResult(myT1,"Order",time=3,"Fat mass percentage change")
plot<-myPlot(myT1,"Order",time=3,"Fat mass percentage change")+
  geom_signif(y_position=c(1.5), xmin=c(1), xmax=c(4),annotation=c("*"), tip_length=0,vjust = 0.5,textsize =3)
myList[[index]]<-plot
index<-index+1

#Extract legend
plot<-ggplot(data=myT1,aes(x=Order,y=-log10(adjustedSpearman)))+geom_boxplot(outlier.size=0.5,aes(fill=factor(Order,labels=c("Bacteroidales (39)","Bifidobacteriales (1)","Burkholderiales (4)","Clostridiales (177)","Coriobacteriales (11)","Desulfovibrionales (2)","Enterobacteriales (1)","Erysipelotrichales (17)","Lactobacillales (3)","Selenomonadales (1)","Verrucomicrobiales (2)"))))+
  theme(axis.text.x = element_blank(),legend.title = element_text(size = 8),
        legend.text = element_text(size=7),legend.key.size = unit(0.7,"line"))+geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black", size=1)+labs(y="-log10(adjusted p-value)",x="")+
  labs(title="Fat mass percentage change\nWeek=3",fill="Order")+geom_signif(y_position=c(3), xmin=c(1), xmax=c(4),annotation=c("**"), tip_length=0,vjust = 0.5,textsize =3)
pdf("legendOrderInteraction.pdf",height = 2,width = 3)
legendp <- cowplot::get_legend(plot)
grid.newpage()
grid.draw(legendp)
dev.off()

tukeyResult(myT1,"Order",time=4,"Fat mass percentage change")
plot<-myPlot(myT1,"Order",time=4,"Fat mass percentage change")+
  geom_signif(y_position=c(1.5), xmin=c(1), xmax=c(4),annotation=c("***"), tip_length=0,vjust = 0.5,textsize =3)
myList[[index]]<-plot
index<-index+1

###Family

tukeyResult(myT1,"Family",time=3,"Fat mass percentage change")
plot<-myPlot(myT1,"Family",time=3,"Fat mass percentage change")
myList[[index]]<-plot
index<-index+1

#Extract legend
plot<-ggplot(data=myT1,aes(x=Family,y=-log10(adjustedSpearman)))+geom_boxplot(outlier.size=0.5,aes(fill=factor(Family,labels=c("Acidaminococcaceae (1)","Alcaligenaceae (3)","Bacteroidaceae (17)","Bifidobacteriaceae (1)","Christensenellaceae (5)","Coriobacteriaceae (11)","Defluviitaleaceae (2)","Desulfovibrionaceae (2)","Enterobacteriaceae (1)","Enterococcaceae (1)","Erysipelotrichaceae (17)","Eubacteriaceae  (3)","Family_XIII (9)","Lachnospiraceae (90)","Lactobacillaceae (1)","Oxalobacteraceae (1)","Peptostreptococcaceae (1)","Porphyromonadaceae (8)","Prevotellaceae (6)","Rikenellaceae (8)","Ruminococcaceae (64)","Streptococcaceae (1)","Verrucomicrobiaceae (2)"))))+
  theme(axis.text.x = element_blank(),legend.text = element_text(size=7) ,legend.title = element_text(size = 8),
        legend.key.size = unit(0.7,"line"))+geom_hline(yintercept=-log10(0.05), linetype="dashed", color="black", size=1)+labs(y="-log10(adjusted p-value)",x="")+
  labs(title="Fat mass percentage change\nWeek=3",fill="Family")+geom_signif(y_position=c(3.2,3.4,3.6,3.8), xmin=c(2,2,2,2), xmax=c(6,13,14,21),annotation=c("*","*","*","*"), tip_length=0,vjust = 0.5,textsize =2)
pdf("legendFamilyInteraction.pdf",height = 2,width = 4)
legendp <- cowplot::get_legend(plot)
grid.newpage()
grid.draw(legendp)
dev.off()

tukeyResult(myT1,"Family",time=4,"Fat mass percentage change")
plot<-myPlot(myT1,"Family",time=4,"Fat mass percentage change")+
  geom_signif(y_position=c(1.5,1.56,1.6), xmin=c(14,20,11), xmax=c(20,21,20),annotation=c("***","***","***"), tip_length=0,vjust = 0.5,textsize =2)
myList[[index]]<-plot
index<-index+1


png("boxPlotsforALLtaxaCorrelationPhenotypePanel.png", units="in", width=10, height=10,res=300)
pdf("boxPlotsforALLtaxaCorrelationPhenotypePanel_new.pdf",height = 10,width = 10)
plot_grid(myList[[1]],myList[[2]],myList[[3]],NULL,myList[[4]],myList[[5]],
          NULL,myList[[6]],myList[[7]],NULL,myList[[8]],myList[[9]],
          nrow=4,ncol=3,scale=0.9)

dev.off()
