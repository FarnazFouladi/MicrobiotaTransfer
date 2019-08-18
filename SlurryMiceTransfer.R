#Farnaz Fouladi
#08-18-2019

#This R code generates Figure 5 and Supplementary Figure 2 and Supplementary Figure 3B. 

rm(list =ls())

library(ggplot2) 
library(grid)
library(gridExtra) 

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)
data<-read.table(paste0(output,"svMeta.txt"),sep = "\t",header=TRUE,check.names=FALSE,na.strings = "NA" ,comment.char="")
finishAbundanceIndex<-which(colnames(data)=="Sample")-1
#Removing low abundant sv
sv<-data[,1:finishAbundanceIndex]
sv<-sv[,colMeans(sv>0)>=0.1]
myT<-cbind(sv,data[,(finishAbundanceIndex+1):ncol(data)])
finishAbundanceIndex<-which(colnames(myT)=="Sample")-1
Donors<-as.character(myT$Donor[myT$Sample.type=="Human.donor"])
SlurryNames<-(myT[myT$Sample.type=="Fecal.slurry",])$Slurry.ID1
SlurryNames<-SlurryNames[SlurryNames!="T1.slurry.EAN40"]#For this slurry there is no mouse

spearman<-vector()
rho<-vector()
time<-vector()
Slurry<-vector()
Mouse<-vector()
outerindex<-1

for (week in 1:4){
  
  myTweek<-myT[myT$Week==week|is.na(myT$Week),]
  
  for (slurry in SlurryNames){
    
    bugH<-myTweek[myTweek$Slurry.ID1==slurry & myTweek$Sample.type=="Fecal.slurry",]
    myTweekSlurry<-myTweek[myTweek$Slurry.ID1==slurry & myTweek$Sample.type=="Mouse.feces",]
    
    for (j in 1:nrow(myTweekSlurry)){
      
      abundanceinSlurry<-vector()
      AbundanceinMice<-vector()
      index<-1  
      
      bugM<-myTweekSlurry[j,]
      
      for (i in 1:finishAbundanceIndex){
        
        if (!(bugH[,i]==0 & bugM[,i]==0)){
          abundanceinSlurry[index]<-bugH[,i]
          AbundanceinMice[index]<-bugM[,i]
          index<-index+1
        }
      }
      
      spearman[outerindex]<-cor.test(abundanceinSlurry,AbundanceinMice,method = "spearman")$p.value
      rho[outerindex]<-cor.test(abundanceinSlurry,AbundanceinMice,method = "spearman")$estimate
      time[outerindex]<-week
      Slurry[outerindex]<-slurry
      Mouse[outerindex]<-as.character(bugM$Sample.ID)
      outerindex<-outerindex+1
    }
  }
}

aframe<-data.frame(Slurry,Mouse,time,spearman,rho)
aframe$Adjustedspearman<-p.adjust(aframe$spearman,method ="BH")
write.table(aframe,"SlurryMiceTransferSV.txt",sep = "\t",row.names = FALSE)

#Plotting
myList<-list()
index<-1

for (i in 1:nrow(aframe)){
  
  data<-data.frame(abundanceinSlurry=as.numeric(myT[myT$Sample.type=="Fecal.slurry" & myT$Slurry.ID1==as.character(aframe$Slurry[i]),1:finishAbundanceIndex]),
                   AbundanceinMice=as.numeric(myT[myT$Sample.type=="Mouse.feces" & myT$Slurry.ID1==as.character(aframe$Slurry[i]) & 
                                                    myT$Week==aframe$time[i] & !is.na(myT$Week) & myT$Sample.ID==as.character(aframe$Mouse[i]) ,1:finishAbundanceIndex]))
  data1<-data[!(data$abundanceinSlurry==0 & data$AbundanceinMice==0), ]
  text1<-paste0("p-value= ",format(aframe$Adjustedspearman[i],digits = 3),"\nrho= ",format(aframe$rho[i],digits = 3),"\nweek= ",aframe$time[i])
  myPlot<-ggplot(data1,aes(x=abundanceinSlurry,y=AbundanceinMice))+geom_point(size=0.5)+labs(title=text1,x=paste0("Abundance of SVs in the slurry "),y="Abundance of SVs in the paired mouse fecal pellet")+theme_gray(base_size = 6)
  myList[[index]]<-myPlot
  index<-index+1
}


#Plot
myVector<-vector()
index<-1
x=1
while (x<535){
  myVector[index]<-x
  x<-x+9
  index<-index+1
}
pdf("SlurryMouseTransferPanel2.pdf")
for (i in myVector){
  if(i<532){
    grid.arrange(myList[[i]],myList[[i+1]],myList[[i+2]],myList[[i+3]],myList[[i+4]],myList[[i+5]],myList[[i+6]],myList[[i+7]],myList[[i+8]],ncol=3) 
  } else {
    grid.arrange(myList[[i]],myList[[i+1]],myList[[i+2]],myList[[i+3]],ncol=3,nrow=3)
  }
}
dev.off()

png("RhoSlurryMice.png", units="in", width=6.5, height=6.5,res=300)
i=1
grid.arrange(myList[[i]],myList[[i+1]],myList[[i+2]],myList[[i+3]],myList[[i+4]],myList[[i+5]],myList[[i+6]],myList[[i+7]],myList[[i+8]],ncol=3)
dev.off()

#Frequency of p-values at each time
pdf("histogramSlurryMiceTransferSpearman.pdf")
par(mfrow=c(2,2))
#P values
for (i in 1:4){
  table1<-aframe[aframe$time==i,]
  hist(table1$Adjustedspearman,breaks = 20,xlab="Adjusted p-value",main=paste0("Week ",i),col = "grey",xlim = c(0,1),ylim = c(0,70),cex.lab=1.5,cex.axis=1.5,cex.main=2)
}
# Rho
for (i in 1:4){
  table1<-aframe[aframe$time==i,]
  hist(table1$rho,breaks = 30,xlab="Rho coefficient",main=paste0("Week ",i),col = "grey",xlim = c(-0.6,0.4),ylim = c(0,15),cex.lab=1.5,cex.axis=1.5,cex.main=2)
}
dev.off()

#Time
png("RhoSlurryMiceTime.png", units="in", width=5, height=5,res=300)
ggplot(aframe,aes(y=rho,x=as.factor(time)))+geom_boxplot()+theme_gray(base_size = 12)+
  labs(x="Week",y="Rho")+ylim(c(-0.7,0.3))
dev.off()
result<-aov(lm(rho~factor(time),data = aframe))
TukeyHSD(result)

png("RhoSlurries.png", units="in", width=5, height=5,res=300) #For publication
i=1
table1<-aframe[aframe$time==i,]
table2<-table1[order(table1$rho),]
table2$Slurry<-as.character(table2$Slurry)
table2$Slurry<-factor(table2$Slurry,levels = unique(table2$Slurry))
Group<-as.vector(sapply(as.character(table2$Slurry),function(x){strsplit(x,"[.]")[[1]][1]}))
plot<-ggplot(data=table2,aes(x=factor(Slurry),y=rho))+geom_boxplot()+
  scale_x_discrete(labels=Group)+
  labs(x="Slurries",y="Rho")+theme_classic(base_size = 12)+
  theme(axis.text.x = element_text(angle=90,size = 5))
print(plot)

dev.off()

#How many SVs significantly show negative correlation at week1?
aframe1<-aframe[aframe$time==1,]
sum(aframe1$Adjustedspearman<0.05) #66
length(aframe1$Adjustedspearman) #135 total paired slurry-mouse
mean(aframe1$rho) #-0.20
sd(aframe1$rho)#0.12
sum(aframe1$rho<0)#133
sum(aframe1$rho>0)#2