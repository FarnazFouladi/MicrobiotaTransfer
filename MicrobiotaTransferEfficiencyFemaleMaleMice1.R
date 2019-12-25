#Farnaz Fouladi
#11-22-19
#Comparing the transfer efficiency between male and female mice. 
#Comparison of Spearman correlation coefficients between male and female mice.

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
myTH<-myT[myT$Sample.type=="Fecal.slurry",]
SlurryNames<-myTH$Slurry.ID1
SlurryNames<-SlurryNames[SlurryNames!="T1.slurry.EAN40"]#For this slurry there is no mouse

#Choose female or male mice
sex="M"
myTMouse_S<-myT[myT$Sex==sex & myT$Sample.type=="Mouse.feces",]
#Removing slurries that were not used for colonization of mice
Slurries_S<-SlurryNames[SlurryNames %in% as.character(myTMouse_S$Slurry.ID1)]

spearman<-vector()
rho<-vector()
time<-vector()
Slurry<-vector()
Mouse<-vector()
outerindex<-1

for (week in 1:4){
  
  myTweek<-myTMouse_S[myTMouse_S$Week==week,]
  
  for (slurry in Slurries_S){
    
    if (slurry %in% myTweek$Slurry.ID1){
    
    bugH<-myTH[myTH$Slurry.ID1==slurry,]
    myTweekSlurry<-myTweek[myTweek$Slurry.ID1==slurry,]
    
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
}

aframe<-data.frame(Slurry,Mouse,time,spearman,rho)
aframe$Adjustedspearman<-p.adjust(aframe$spearman,method ="BH")
write.table(aframe,paste0("SlurryMiceTransferSV_",sex,".txt"),sep = "\t",row.names = FALSE)

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
while (x<length(myList)){
  myVector[index]<-x
  x<-x+9
  index<-index+1
}
pdf(paste0("SlurryMouseTransferPanel2_",sex,".pdf"))
for (i in myVector){
  if(i<length(myList)-1){
    grid.arrange(myList[[i]],myList[[i+1]],myList[[i+2]],myList[[i+3]],myList[[i+4]],myList[[i+5]],myList[[i+6]],myList[[i+7]],myList[[i+8]],ncol=3) 
  } else {
    grid.arrange(myList[[i]],myList[[i+1]],ncol=3,nrow=3)
  }
}
dev.off()


#Box plot for comparing transfer efficiency between male and female mice

aframe_M<-read.table(paste0(output,"SlurryMiceTransferSV_M.txt"),header = TRUE,sep="\t")
aframe_M$sex<-rep("Male",nrow(aframe_M))
aframe_F<-read.table(paste0(output,"SlurryMiceTransferSV_F.txt"),header = TRUE,sep="\t")
aframe_F$sex<-rep("Female",nrow(aframe_F))
aframe_all<-rbind(aframe_M,aframe_F)
aframe_all$time<-factor(aframe_all$time)
#Comparing rho coeficients between female and male mice:
theme_set(theme_classic(base_size = 5))
plot1<-ggplot(aframe_all,aes(x=time,y=rho,color=as.factor(sex)))+geom_boxplot(outlier.shape = NA)+
  labs(y="Rho",color="")+geom_jitter(shape=16, size =0.2,position = position_jitterdodge(0.1),aes(color=sex))+
  annotate(geom="text", x=2, y=0.1, label="adjusted p=0.15",size=2,
           color="black")
  

aframe_all1<-aframe_all[aframe_all$time==1,]
aframe_all2<-aframe_all[aframe_all$time==2,]
aframe_all3<-aframe_all[aframe_all$time==3,]
aframe_all4<-aframe_all[aframe_all$time==4,]

pval<-vector()
pval[1]<-t.test(aframe_all1[aframe_all1$sex=="Female",]$rho,aframe_all1[aframe_all1$sex=="Male",]$rho)$p.value
pval[2]<-t.test(aframe_all2[aframe_all2$sex=="Female",]$rho,aframe_all2[aframe_all2$sex=="Male",]$rho)$p.value
pval[3]<-t.test(aframe_all3[aframe_all3$sex=="Female",]$rho,aframe_all3[aframe_all3$sex=="Male",]$rho)$p.value
pval[4]<-t.test(aframe_all4[aframe_all4$sex=="Female",]$rho,aframe_all4[aframe_all4$sex=="Male",]$rho)$p.value
p.adjust(pval,method = "BH")

pdf("ComparisonOfMaleAndFemaleMice_1.pdf",width = 3,height = 3)
print(plot1)
dev.off()
