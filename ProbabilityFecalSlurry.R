#Farnaz Fouladi
#08-18-2019

#This R code generates the probability plots in Figure 2. 

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
# Attaching taxanomy
taxa<-read.table("/Users/farnazfouladi/CarrollMouseTransfer/taxForwardReads.txt",header = TRUE,sep = "\t")
rownames(taxa)<-sapply(seq(1:nrow(taxa)),function(x){x=paste0("SV",x)})
taxa1<-taxa[intersect(rownames(taxa),colnames(sv)),]
taxa1$names<-rownames(taxa1)

# Ordering SVs based on their abundance in fecal samples
myTfecal<-myT[myT$Sample.type=="Human.donor",]
names<-colnames(myTfecal[,1:finishAbundanceIndex])
sum<-colSums(myTfecal[,1:finishAbundanceIndex])
df<-data.frame(names,sum)
df<-df[order(df$sum,decreasing = TRUE),]
df$rank<-rank(df$sum)
# Ordering SVs based on their abundance in Slurry
myTslurry<-myT[myT$Sample.type=="Fecal.slurry",]
names<-colnames(myTslurry[,1:finishAbundanceIndex])
sum<-colSums(myTslurry[,1:finishAbundanceIndex])
dfslurry<-data.frame(names,sum)
dfslurry<-dfslurry[order(dfslurry$sum,decreasing = TRUE),]
dfslurry$rank<-rank(dfslurry$sum)
rankGroup<-sapply(dfslurry$rank,function(x){if (x<100) {x=1} else if (x>=100 & x<200) {x=2} else {x=3} })
dfslurry$rankGroup<-rankGroup
#Ordering the sv columns in the myT table based on the abundance of sv in fecal samples
meta<-myT[,(finishAbundanceIndex+1):ncol(myT)]
sv<-sv[,match(df$names,colnames(sv))]
myT<-cbind(sv,meta)
# Deteriming the probability of shared/only in mouse/only in slurry/nogroup 
myT1<-myT[myT$Sample.type=="Human.donor"| myT$Sample.type=="Fecal.slurry",]
Donor<-as.character(myT[myT$Sample.type=="Human.donor",]$Donor)
Slurry<-as.character(myT[myT$Sample.type=="Fecal.slurry",]$Slurry.ID1)
DF<-data.frame(probability=c("pShared","pSlurry","pFecal","noGroup"))

for (i in 1:finishAbundanceIndex){
  
  shared<-0
  onlySlurry<-0
  onlyFecal<-0
  noGroup<-0
  slurries<-0
  for (donor in Donor){
    
    bugH<-myT1[myT1$Sample.type=="Human.donor" & myT1$Donor==donor,i]
    myTDonor<-myT1[myT1$Sample.type=="Fecal.slurry" &  myT1$Donor==donor,]
    slurries<-slurries+nrow(myTDonor)
    
    for(slurry in 1:nrow(myTDonor)){
      
      bugS<-myTDonor[slurry,i]
      
      if (bugH>0 & bugS>0){ shared<-shared+1} 
      else if (bugH>0 & bugS==0) { onlyFecal<-onlyFecal+1}
      else if (bugH==0 & bugS>0) {onlySlurry<-onlySlurry+1}
      else { noGroup<-noGroup+1}
    }
  }
  
  shared<-shared/slurries
  onlySlurry<-onlySlurry/slurries
  onlyFecal<-onlyFecal/slurries
  noGroup<-noGroup/slurries
  
  values<-c(shared,onlySlurry,onlyFecal,noGroup)
  DF<-cbind(DF,values)
  colnames(DF)[which(colnames(DF)=="values")]<-colnames(myT)[i]
}

DF<-as.data.frame(t(DF))
colnames(DF)<-c("pShared","pSlurry","pFecal","noGroup")
DF<-DF[-1,]
DF$rank<-df$rank
write.table(DF,"probabilityforEachSeqFecaltoSlurry.txt",sep="\t")

DF$names<-rownames(DF)
DF1<-merge(DF,dfslurry,by.x ="names",by.y ="names",all.x = TRUE,sort = FALSE)
DF2<-merge(DF1,taxa1,by.x ="names",by.y ="names",all.x = TRUE,sort = FALSE)

# Extracting Legend

SH<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF$pShared))))+geom_point(aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(title="Shared",x="Rank abundance of SVs in fecal samples",y="Probability")+
  scale_size_ordinal(range=c(1,3),labels=c("6-99","100-199","200-279"),name="Rank abundance of SVs in slurries")+
  theme(legend.text = element_text(size = 8))
png("legendProbabilityFecalSlurry.png", units="in", width=5, height=5,res=300)
legend <- cowplot::get_legend(SH)
grid.newpage()
grid.draw(legend)
dev.off()

#plots
theme_set(theme_grey(base_size = 8))
SH<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF$pShared))))+geom_point(col="indianred1",aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(x="Rank abundance of SVs in human fecal samples",y="Probability of being shared")+scale_size_ordinal(range=c(1,3))+theme(legend.position = "none")

S<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF$pSlurry))))+geom_point(col="orchid",aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(x="Rank abundance of SVs in human fecal samples",y="Probability of being only in slurries")+scale_size_ordinal(range=c(1,3))+theme(legend.position = "none")

F<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF$pFecal))))+geom_point(col="blue",aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(x="Rank abundance of SVs in human fecal samples",y="Probability of being only in human fecal samples")+scale_size_ordinal(range=c(1,3))+theme(legend.position = "none")

png("probabilityHumanSlurryTransfer.png", units="in", width=6, height=6,res=300)
plot_grid(SH,S,F, align = 'h',ncol=2,nrow=2,label_size = 12,scale = 0.9)
dev.off()
