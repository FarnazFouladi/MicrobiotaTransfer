#Farnaz Fouladi
#08-18-2019

#This R code generates Figure 3 and Supplementary Figure 1. 

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
myT1<-myT[myT$Sample.type=="Human.donor"|myT$Sample.type=="Fecal.slurry",]

pSpearman<-vector()
rho<-vector()
donorName<-vector()
slurryName<-vector()
outerindex<-1

for (donor in Donors){
  
  bugFecal<-myT1[myT1$Sample.type=="Human.donor" & myT1$Donor==donor,]
  donorTable<-myT1[myT1$Donor==donor & myT1$Sample.type=="Fecal.slurry",]
  
  for ( slurry in 1:nrow(donorTable) ){
    
    fecalAbundance<-vector()
    slurryAbundance<-vector()
    index<-1
    bugSlurry<-donorTable[slurry,]
    
    for (i in 1:finishAbundanceIndex){
      
      if (!(bugFecal[,i]==0 & bugSlurry[,i]==0)){
        
        fecalAbundance[index]<-bugFecal[,i]
        slurryAbundance[index]<-bugSlurry[,i]
        index<-index+1
      }
    }
    pSpearman[outerindex]<-cor.test(fecalAbundance,slurryAbundance,method="spearman")$p.value
    rho[outerindex]<-cor.test(fecalAbundance,slurryAbundance,method="spearman")$estimate
    donorName[outerindex]<-donor
    slurryName[outerindex]<-as.character(bugSlurry$Slurry.ID)
    outerindex<-outerindex+1
  }
}
aframe<-data.frame(donorName,slurryName,pSpearman,rho)
aframe$adjustedpSpearman<-p.adjust(aframe$pSpearman,method = "BH")
write.table(aframe, file="FecalSlurryTransfer1.txt", sep="\t",row.names=FALSE)

#Plotting

myList<-list()
index<-1
for (i in 1:nrow(aframe)){
  
  data<-data.frame(fecalAbundance=as.numeric(myT1[myT1$Sample.type=="Human.donor" & myT1$Donor==as.character(aframe$donorName[i]),1:finishAbundanceIndex]),
                   slurryAbundance=as.numeric(myT1[myT1$Sample.type=="Fecal.slurry" & myT1$Donor==as.character(aframe$donorName[i]) & myT1$Slurry.ID==as.character(aframe$slurryName[i]),1:finishAbundanceIndex]))
  data1<-data[!(data$fecalAbundance==0 & data$slurryAbundance==0), ]
  
  text1<-paste0("p-value= ",format(aframe$adjustedpSpearman[i],digits = 3),"\nrho= ",format(aframe$rho[i],digits = 3))
  myPlot<-ggplot(data1,aes(x=fecalAbundance,y=slurryAbundance))+geom_point(size=0.5)+labs(title=text1,x=paste0("Abundance of SVs in the human fecal sample"),y="Abundance of SVs in the paired slurry")+theme_gray(base_size = 6)
  myList[[index]]<-myPlot
  index<-index+1
  
}

pdf("RhoFecalSlurry.pdf")
#For figure in the paper i=1:
#png("RhoFecalSlurry.png", units="in", width=6, height=6,res=300) #i=1
for (i in c(1,10,19,28,37,46,55)){
  grid.arrange(myList[[i]],myList[[i+1]],myList[[i+2]],myList[[i+3]],myList[[i+4]],myList[[i+5]],myList[[i+6]],myList[[i+7]],myList[[i+8]],ncol=3,nrow=3) 
}
dev.off()

png("RhoFecalSlurryB.png", units="in", width=6, height=6,res=300)
aframe1<-aframe[order(aframe$rho),]
aframe1$donorName<-as.character(aframe1$donorName)
aframe1$donorName<-factor(aframe1$donorName,levels = unique(aframe1$donorName))
names<-sapply(as.character(aframe1$donorName),function(x){if (x=="AN703HC"|x=="AN703T1"|x=="AN703T2") {strsplit(x,"AN...")[[1]][2]}else{strsplit(x,"AN..")[[1]][2]}})
ggplot(data=aframe1,aes(x=factor(donorName),y=rho))+geom_boxplot()+scale_x_discrete(labels=names)+
  labs(x="Human donors",y="Rho")+theme_classic(base_size = 14)
dev.off()

#Rho mean
mean(aframe$rho) #0.77
sd(aframe$rho)   #0.15

