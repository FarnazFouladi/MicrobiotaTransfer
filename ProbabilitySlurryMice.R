#Farnaz Fouladi
#08-18-2019

#This R code generates the probability plots in Figure 4 and Supplementary Tables 1 and 2. 

rm(list =ls())

library(ggplot2) 
library(grid)
library(gridExtra) 

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)
data<-read.table(paste0(output,"svMeta.txt"),sep = "\t",header=TRUE,check.names=FALSE,na.strings = "NA" ,comment.char="")
#Slecect the time point
week=1

finishAbundanceIndex<-which(colnames(data)=="Sample")-1
#Removing low abundant sv
sv<-data[,1:finishAbundanceIndex]
sv<-sv[,colMeans(sv>0)>=0.1]
myT<-cbind(sv,data[,(finishAbundanceIndex+1):ncol(data)])
finishAbundanceIndex<-which(colnames(myT)=="Sample")-1
# Attaching taxanomy
taxa<-read.table("/Users/farnazfouladi/CarrollMouseTransfer/taxForwardReads.txt",header = TRUE,sep = "\t")
rownames(taxa)<-sapply(seq(1:nrow(taxa)),function(x){x=paste0("SV",x)})
taxa1<-taxa[intersect(rownames(taxa),colnames(sv)),]
taxa1$names<-rownames(taxa1)
# Ordering SVs based on their abundance in slurry
myTs<-myT[myT$Sample.type=="Fecal.slurry",]
names<-colnames(myTs[,1:finishAbundanceIndex])
sum<-colSums(myTs[,1:finishAbundanceIndex])
df<-data.frame(names,sum)
df<-df[order(df$sum,decreasing = TRUE),]
df$rank<-rank(df$sum)
# Ordering SVs based on their abundance in mice
myTMouse<-myT[myT$Sample.type=="Mouse.feces" & myT$Week==week,]
names<-colnames(myTMouse[,1:finishAbundanceIndex])
sum<-colSums(myTMouse[,1:finishAbundanceIndex])
dfMouse<-data.frame(names,sum)
dfMouse<-dfMouse[order(dfMouse$sum,decreasing = TRUE),]
dfMouse$rank<-rank(dfMouse$sum)
rankGroup<-sapply(dfMouse$rank,function(x){if (x<100) {x=1} else if (x>=100 & x<200) {x=2} else {x=3} })
dfMouse$rankGroup<-rankGroup
#Ordering the sv columns in the myT table based on the abundance of sv in slurry
meta<-myT[,(finishAbundanceIndex+1):ncol(myT)]
sv<-sv[,match(df$names,colnames(sv))]
myT<-cbind(sv,meta)
# Deteriming the probability of shared/only in mouse/only in slurry/nogroup for each sv at week1
myT1<-myT[myT$Sample.type=="Mouse.feces"| myT$Sample.type=="Fecal.slurry",]
Slurry<-(myT1[myT1$Sample.type=="Fecal.slurry",])$Slurry.ID1
Slurry<-Slurry[Slurry%in% myTMouse$Slurry.ID1]

DF<-data.frame(probability=c("pShared","pMouse","pSlurry","noGroup"))

for (i in 1:finishAbundanceIndex){
  
  shared<-0
  onlyMouse<-0
  onlySlurry<-0
  noGroup<-0
  number<-0
  
  for (slurry in Slurry){
    
    bugS<-myT1[myT1$Sample.type=="Fecal.slurry" & myT1$Slurry.ID1==slurry,i]
    bugM<-myT1[myT1$Sample.type=="Mouse.feces" & myT1$Slurry.ID1==slurry & myT1$Week==week,i]
    number<-number+length(bugM)
    
    for (j in 1:length(bugM))
    {
      bugM_individual<-bugM[j]
      
      if (bugS>0 & bugM_individual>0 ){shared<-shared+1}
      else if (bugS>0 & bugM_individual==0) { onlySlurry<-onlySlurry+1}
      else if (bugS==0 & bugM_individual>0) {onlyMouse<-onlyMouse+1}
      else { noGroup<-noGroup+1}
    }
  }
  
  shared<-shared/number
  onlyMouse<-onlyMouse/number
  onlySlurry<-onlySlurry/number
  noGroup<-noGroup/number
  
  values<-c(shared,onlyMouse,onlySlurry,noGroup)
  DF<-cbind(DF,values)
  colnames(DF)[which(colnames(DF)=="values")]<-colnames(myT)[i]
}

DF<-as.data.frame(t(DF))
colnames(DF)<-c("pShared","pMouse","pSlurry","noGroup")
DF<-DF[-1,]
DF$rank<-df$rank
write.table(DF,paste0("probabilityforEachSeqWeek",week,".txt"),sep="\t")

DF$names<-rownames(DF)
DF1<-merge(DF,dfMouse,by.x ="names",by.y ="names",all.x = TRUE,sort = FALSE)
DF2<-merge(DF1,taxa1,by.x ="names",by.y ="names",all.x = TRUE,sort = FALSE)
row.names(DF2)<-DF2$names
write.table(DF2,paste0("probabilityforEachSeqWeek",week,"withTaxa.txt"),sep="\t")

# Extracting Legend

theme_set(theme_classic(base_size = 8))
SH<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF1$pShared))))+geom_point(aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(title="Shared",x="Rank abundance of SVs in slurries",y="Probability")+scale_size_ordinal(range=c(1,3),labels=c("1-99","100-199","200-279"),name="Rank abundance of SVs in mouse fecal pellets")+
  theme(legend.text = element_text(size = 8))

png("legendProbabilitySlurryMice.png", units="in", width=5, height=5,res=300)
legend <- cowplot::get_legend(SH)
grid.newpage()
grid.draw(legend)
dev.off()


theme_set(theme_grey(base_size = 8))
SH<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF1$pShared))))+geom_point(col="indianred1",aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(x="Rank abundance of SVs in slurries",y="Probability of being shared")+scale_size_ordinal(range=c(1,3))+theme(legend.position = "none")

S<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF1$pSlurry))))+geom_point(col="blue",aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(x="Rank abundance of SVs in slurries",y="Probability of being only in slurries")+scale_size_ordinal(range=c(1,3))+theme(legend.position = "none")


M<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF1$pMouse))))+geom_point(col="orchid",aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(x="Rank abundance of SVs in slurries",y="Probability of being only in mouse fecal pellets")+scale_size_ordinal(range=c(1,3))+theme(legend.position = "none")

png(paste0("probabilitySlurryMiceTransferWeek",week,".png"), units="in", width=6, height=6,res=300)
plot_grid(SH,M,S, align = 'h',ncol=2,nrow=2,label_size = 12,scale = 0.9)
dev.off()


# High abundant slurry SVs that had lower probbaility of jumping into mice (Supplementary Table 1):
DF3<-DF2[as.numeric(as.character(DF2$pShared))<0.25 & DF2$rank.x>200,] #32 bugs
rownames(DF3)<-DF3$names
DF5<-DF2[as.numeric(as.character(DF2$pShared))>=0.25 & DF2$rank.x>200,]
write.table(DF3,paste0("taxawithlowChanceofTransferWeek",week,".txt"),sep="\t",row.names = FALSE)
#Fisher test
DF_Fisher<-data.frame(Group=c("high pShared","low pShared"))
Taxa<-c("Firmicutes","Bacteroidetes","Proteobacteria","Actinobacteria","Verrucomicrobia")
for (taxa in Taxa){
  myMatrix<-matrix(c(sum(DF5$Phylum==taxa),sum(DF3$Phylum==taxa),(length(DF5$Phylum)-sum(DF5$Phylum==taxa)),(length(DF3$Phylum)-sum(DF3$Phylum==taxa))),
                   nrow=2,dimnames = list(transfer=c("high pShared","low pShared"),taxa=c(taxa,paste0("non",taxa))))
  P<-fisher.test(myMatrix)$p.value  
  e<-fisher.test(myMatrix)$estimate
  df<-as.data.frame(myMatrix) 
  DF_Fisher<-cbind(DF_Fisher,df,Names=c("pValue","oddsRatio"),fisherTest=c(P,e))
}
DF_Fisher1<-DF_Fisher[1,which(colnames(DF_Fisher)=="fisherTest")]
AdjustedP<-p.adjust(DF_Fisher1,method = "BH")
DF_Fisher2<-data.frame(taxa=Taxa,AsjustedPvals=AdjustedP)
write.table(DF_Fisher,paste0("lowVersusHighTransferFishertestWeek",week,".txt"),sep = "\t")
write.table(DF_Fisher2,paste0("lowVersusHighTransferFishertestAdjustedPWeek",week,".txt"),sep = "\t")

# low abundant slurry SVs that had high probbaility of finding in only mice pmouse>0.25 (Supplementary Table 2):
DF4<-DF2[as.numeric(as.character(DF2$pMouse))>0.25,]
write.table(DF4,"taxawithhigherChanceofbeingFoundOnlyinMice.txt",sep="\t",row.names = FALSE)
DF6<-DF2[as.numeric(as.character(DF2$pMouse))<=0.25,]

#Fisher test
DF6<-DF6[-which(is.na(DF6$Phylum)),]
DF_Fisher<-data.frame(Group=c("high pmouse","low pmouse"))
Taxa<-c("Firmicutes","Bacteroidetes","Proteobacteria","Actinobacteria","Verrucomicrobia")
for (taxa in Taxa){
  myMatrix<-matrix(c(sum(DF4$Phylum==taxa),sum(DF6$Phylum==taxa),(length(DF4$Phylum)-sum(DF4$Phylum==taxa)),(length(DF6$Phylum)-sum(DF6$Phylum==taxa))),nrow=2,dimnames = list(transfer=c("high pmouse","low pmouse"),taxa=c(taxa,paste0("non",taxa))))
  P<-fisher.test(myMatrix)$p.value  
  e<-fisher.test(myMatrix)$estimate
  df<-as.data.frame(myMatrix) 
  DF_Fisher<-cbind(DF_Fisher,df,Names=c("pValue","oddsRatio"),fisherTest=c(P,e))
}
DF_Fisher1<-DF_Fisher[1,which(colnames(DF_Fisher)=="fisherTest")]
AdjustedP<-p.adjust(DF_Fisher1,method = "BH")
DF_Fisher2<-data.frame(taxa=Taxa,AsjustedPvals=AdjustedP)
write.table(DF_Fisher,paste0("taxawithhigherChanceofbeingFoundOnlyinMiceFishertestWeek",week,".txt"),sep = "\t")
write.table(DF_Fisher2,paste0("taxawithhigherChanceofbeingFoundOnlyinMiceFishertestAdjustedPWeek",week,".txt"),sep = "\t")

# SVs that are correlated between human and mice at week1:
myData<-read.table(paste0(output,"SlurryVsMouseCor.txt"),header = TRUE,sep = "\t")
myData1<-myData[myData$adjustedSpearman<0.05 & myData$time==1,]
CorrelatedSV<-as.character(unique(myData1$bugnames))
rownames(DF2)<-DF2$names
DF7<-DF2[intersect(rownames(DF2),CorrelatedSV),]
write.table(DF7,"CorrelatedSVprobabs.txt",row.names = FALSE,sep="\t")


# SVs that are not correlated between human and mice at week 1:
myData2<-myData[myData$adjustedSpearman>0.05&myData$time==1,]
nonCorrelatedSV<-as.character(myData$bugnames)
DF8<-DF2[intersect(rownames(DF2),nonCorrelatedSV),]
write.table(DF8,"NonCorrelatedSVprobabs.txt",row.names = FALSE,sep="\t")

# How many SVs are never been slurry
DF9<-DF2[as.numeric(as.character(DF2$pShared))==0 & as.numeric(as.character(DF2$pSlurry))==0 & as.numeric(as.character(DF2$pMouse))>0 ,]
write.table(DF5,"onlyinMice.txt",row.names = FALSE,sep="\t")