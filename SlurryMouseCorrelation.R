#Farnaz Fouladi
#08-18-2019

#This R code generates Supplementary Figures 3C Figure 4 and Table 1 and Supplementary Table3. 

rm(list =ls())

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
Slurry<-(myT[myT$Sample.type=="Fecal.slurry",])$Slurry.ID1
Slurry<-Slurry[Slurry!="T1.slurry.EAN40"]#For this slurry there is no mouse

myT1<-myT[myT$Sample.type=="Mouse.feces"| myT$Sample.type=="Fecal.slurry",]

spearman<-vector()
rho<-vector()
bugnames<-vector()
time<-vector()
outerindex<-1
pdf("SlurryMouseCor.pdf")
par(mfrow=c(3,2))

for ( i in 1:finishAbundanceIndex){
  
  bug1<-myT1[,i]
  if(  sum(bug1!=0) > nrow(myT1)/10) {
    for (week in 1:4){
      myTweek<-myT1[is.na(myT1$Week)|myT1$Week==week,]
      
      human<-vector()
      mice<-vector()
      donorName<-vector()
      index<-1
      donor<-myTweek[,"Donor"]
      bug<-myTweek[,i]
      
      for (slurry in Slurry){
        bugM<-bug[myTweek$Slurry.ID1==slurry & myTweek$Sample.type=="Mouse.feces"]
        bugH<-bug[myTweek$Slurry.ID1==slurry & myTweek$Sample.type=="Fecal.slurry"]
        slurryDonor<-donor[myTweek$Slurry.ID1==slurry & myTweek$Sample.type=="Mouse.feces"]
        
        if (length(bugM)>0){
          for (j in 1:length(bugM)){
            mice[index]<-bugM[j]
            human[index]<-bugH
            donorName[index]<-as.character(slurryDonor[j])
            index<-index+1
          }
        } else {
          print(paste("no mouse for",slurry))
        }
      }
      if (sum(human>0)>0&sum(mice>0)>0){
        
        spearman[outerindex]<-cor.test(human,mice,method="spearman")$p.value
        rho[outerindex]<-cor.test(human,mice,method="spearman")$estimate
        time[outerindex]<-week
        bugnames[outerindex]<-colnames(myTweek)[i]
        
        atext<-paste0(bugnames[outerindex], "\np-value= ", format(spearman[outerindex],digits=3)," rho= ",format(rho[outerindex],digits=3),
                      "\nweek= ",week)
        col=c("red","blue","green","orange","purple","pink","black","lightblue","lightgreen","hotpink","cyan","orchid","tan","grey","gold","firebrick","deeppink")
        col1=col[factor(donorName)]
        plot (jitter(human),jitter(mice),main = atext,cex.main=0.75,col=col1,pch=16,xlab = "Abundance in slurries",ylab = "Abundance in paired mouse fecal pellets")
      } else {
        spearman[outerindex]<-1
        rho[outerindex]<-0
        bugnames[outerindex]<-colnames(myTweek)[i]
        time[outerindex]<-week
        
        plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
        mtext(paste0(bugnames[outerindex]," \nNot detected in slurries or mice\nweek ",week), at=0.4, cex=0.8)
        
      }
      if (sum(human>0)==0&sum(mice>0)==0){
        print("neither in mice or human")
      }
      outerindex<-outerindex+1
    }
    
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend(0.3,1.2, legend =levels(factor(donorName)), pch=16, pt.cex=1, cex=1, bty='n',
           col = col[1:length(levels(factor(donorName)))],xpd = TRUE)
    mtext("Donors", at=0.1, cex=1)
    
    plot.new()
  }
}
dev.off()
aframe<-data.frame(bugnames,time,rho,spearman)
aframe<-aframe[order(aframe$spearman),]
aframe$adjustedSpearman<-p.adjust(aframe$spearman,method = "BH")
write.table(aframe, file=paste("SlurryVsMouseCor.txt",sep=""), sep="\t",row.names=FALSE)

#How many SVs were significantly positively correlated:
aframeT<-aframe[aframe$adjustedSpearman<0.05 & aframe$rho>0,]
summary(aframeT$rho) #0.18-0.99
aframe1<-aframe[aframe$time==1 & aframe$adjustedSpearman<0.05 & aframe$rho>0,] #185
aframe2<-aframe[aframe$time==2 & aframe$adjustedSpearman<0.05 & aframe$rho>0,] #179
aframe3<-aframe[aframe$time==3 & aframe$adjustedSpearman<0.05 & aframe$rho>0,] #173
aframe4<-aframe[aframe$time==4 & aframe$adjustedSpearman<0.05 & aframe$rho>0,] #188
nrow(aframe4)/sum(aframe$time==4) #0.69
summary(aframe4$rho) #0.18-0.99


#Histogram of rho coefficients 
pdf("HistogramRho.pdf")

#png("HistogramRhoTime4.png", units="in", width=10, height=10,res=300)
#par(mfrow=c(2,2))

for (i in 1:4){
  aframeweek<-aframe[aframe$time==i,]
  x<-aframeweek$rho
  h<-hist(x,breaks=100,plot=FALSE)
  
  if (i==1){
    cuts<-cut(h$breaks,c(-Inf,-0.38,-0.24,0.18,0.23,0.29,Inf))
    plot(h, col=ifelse(cuts=="(-0.24,0.18]","red",ifelse(cuts=="(0.18,0.23]","orchid",ifelse((cuts=="(0.23,0.29]"|cuts=="(-0.38,-0.24]"),"blue","turquoise1"))),main =paste0("Week ",i),xlab = "Rho",ylim=c(0,14),cex.lab=1,cex.axis=1,cex.main=1)
  } else if (i==2){
    cuts<-cut(h$breaks,c(-Inf,-0.34,0.18,0.24,0.30,Inf))
    plot(h, col=ifelse(cuts=="(-0.34,0.18]","red",ifelse(cuts=="(0.18,0.24]","orchid",ifelse(cuts=="(0.24,0.3]","blue","turquoise1"))),main =paste0("Week ",i),xlab = "Rho",ylim=c(0,14),cex.lab=1,cex.axis=1,cex.main=1)
  } else if (i==3){
    cuts<-cut(h$breaks,c(-Inf,-0.41,0.18,0.23,0.29,Inf))
    plot(h, col=ifelse(cuts=="(-0.41,0.18]","red",ifelse(cuts=="(0.18,0.23]","orchid",ifelse(cuts=="(0.23,0.29]","blue","turquoise1"))),main =paste0("Week ",i),xlab = "Rho",ylim=c(0,14),cex.lab=1,cex.axis=1,cex.main=1)
  }else {
    cuts<-cut(h$breaks,c(-Inf,-0.40,-0.19,0.18,0.23,0.29,Inf))
    plot(h, col=ifelse(cuts=="(-0.19,0.18]","red",ifelse((cuts=="(0.18,0.23]"|cuts=="(-0.4,-0.19]"),"orchid",ifelse(cuts=="(0.23,0.29]","blue","turquoise1"))),main =paste0("Week ",i),xlab = "Rho",ylim=c(0,14),cex.lab=1,cex.axis=1,cex.main=1)
  }
  par(xpd=TRUE)
  legend(0.6,12,c("adj.p>0.05","0.01<adj.p<0.05","0.001<adj.p<0.01","adj.p<0.001"),pch=15,col=c("red","orchid","blue","turquoise1"),cex = 0.8,
         title = "Significance level")
}
dev.off()

#Table for SVs with high transfer efficiency at week four

taxa<-read.table(paste0(input,"taxForwardReads.txt"),header = TRUE,sep = "\t")
rownames(taxa)<-sapply(seq(1:nrow(taxa)),function(x){paste0("SV",x)})
taxa<-cbind(rownames(taxa),taxa)
aframe4<-aframe[aframe$time==4 & aframe$adjustedSpearman<0.05 & aframe$rho>0,]
mergedDtata<-merge(aframe4,taxa,by.x="bugnames",by.y="rownames(taxa)",sort=FALSE)
write.table(mergedDtata, file=paste("SigSlurryVsMouseCorWeek4.txt",sep=""), sep="\t",row.names=FALSE)

#Boxplot for rho at differet time points
theme_set(theme_gray(base_size = 12))
ggplot(data=aframe,aes(x=factor(time),,y=rho))+geom_boxplot()+labs(x="Week",y="Rho")
