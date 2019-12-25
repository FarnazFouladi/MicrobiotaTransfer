#Farnaz Fouladi
#12-1-2019
#Percent of sequence variants that were shared or detected only in human fecal samples 
#or in mouse fecal pellets.

rm(list =ls())

library(ggplot2) 
library(grid)
library(gridExtra) 

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)
data<-read.table(paste0(output,"svMeta.txt"),sep = "\t",header=TRUE,check.names=FALSE,na.strings = "NA" ,comment.char="")
#Slecect the time point
week=4

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
# Ordering SVs based on their abundance in Fecal samples
myTs<-myT[myT$Sample.type=="Human.donor",]
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
#Ordering the sv columns in the myT table based on the abundance of sv in fecal samples
meta<-myT[,(finishAbundanceIndex+1):ncol(myT)]
sv<-sv[,match(df$names,colnames(sv))]
myT<-cbind(sv,meta)
# Deteriming the probability of shared/only in mouse/only in fecal/nogroup for each sv at week1
myT1<-myT[myT$Sample.type=="Mouse.feces"| myT$Sample.type=="Human.donor",]
Donors<-(myT1[myT1$Sample.type=="Human.donor",])$Donor


DF<-data.frame(probability=c("pShared","pMouse","pFecal","noGroup"))

for (i in 1:finishAbundanceIndex){
  
  shared<-0
  onlyMouse<-0
  onlyFecal<-0
  noGroup<-0
  number<-0
  
  for (donor in Donors){
    
    bugS<-myT1[myT1$Sample.type=="Human.donor" & myT1$Donor==donor,i]
    bugM<-myT1[myT1$Sample.type=="Mouse.feces" & myT1$Donor==donor & myT1$Week==week,i]
    number<-number+length(bugM)
    
    for (j in 1:length(bugM))
    {
      bugM_individual<-bugM[j]
      
      if (bugS>0 & bugM_individual>0 ){shared<-shared+1}
      else if (bugS>0 & bugM_individual==0) { onlyFecal<-onlyFecal+1}
      else if (bugS==0 & bugM_individual>0) {onlyMouse<-onlyMouse+1}
      else { noGroup<-noGroup+1}
    }
  }
  
  shared<-shared/number
  onlyMouse<-onlyMouse/number
  onlyFecal<-onlyFecal/number
  noGroup<-noGroup/number
  
  values<-c(shared,onlyMouse,onlyFecal,noGroup)
  DF<-cbind(DF,values)
  colnames(DF)[which(colnames(DF)=="values")]<-colnames(myT)[i]
}

DF<-as.data.frame(t(DF))
colnames(DF)<-c("pShared","pMouse","pFecal","noGroup")
DF<-DF[-1,]
DF$rank<-df$rank
write.table(DF,paste0("probabilityforEachSeqWeek",week,"FecalMouse.txt"),sep="\t")

DF$names<-rownames(DF)
DF1<-merge(DF,dfMouse,by.x ="names",by.y ="names",all.x = TRUE,sort = FALSE)
DF2<-merge(DF1,taxa1,by.x ="names",by.y ="names",all.x = TRUE,sort = FALSE)
row.names(DF2)<-DF2$names
write.table(DF2,paste0("probabilityforEachSeqWeek",week,"withTaxaFecalMouse.txt"),sep="\t")

# Extracting Legend

theme_set(theme_classic(base_size = 7))
SH<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF1$pShared))))+geom_point(aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(title="Shared",x="Rank abundance of SVs in slurries",y="Probability")+scale_size_ordinal(range=c(1,3),labels=c("1-99","100-199","200-279"),name="Rank abundance of SVs in mouse fecal pellets")+
  theme(legend.text = element_text(size = 8))

pdf("legendProbabilitySlurryMice.pdf", width=2.1, height=1)
legend <- cowplot::get_legend(SH)
grid.newpage()
grid.draw(legend)
dev.off()


theme_set(theme_grey(base_size = 8))
SH<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF1$pShared))))+geom_point(col="indianred1",aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(x="Rank abundance of SVs in human fecal samples",y="Probability of being shared")+scale_size_ordinal(range=c(1,3))+theme(legend.position = "none")

S<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF1$pFecal))))+geom_point(col="blue",aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(x="Rank abundance of SVs in human fecal samples",y="Probability of being only in human fecal samples")+scale_size_ordinal(range=c(1,3))+theme(legend.position = "none")


M<-ggplot(data = DF1,aes(x=rank.x,y=as.numeric(as.character(DF1$pMouse))))+geom_point(col="orchid",aes(size=factor(rankGroup)),shape=1)+
  ylim(y=c(0,1))+labs(x="Rank abundance of SVs in human fecal samples",y="Probability of being only in mouse fecal pellets")+scale_size_ordinal(range=c(1,3))+theme(legend.position = "none")

png(paste0("probabilitySlurryMiceTransferWeek",week,"HumanFecalMouse.png"), units="in", width=6, height=6,res=300)
plot_grid(SH,M,S, align = 'h',ncol=2,nrow=2,label_size = 12,scale = 0.9)
dev.off()

#All time points together (SH1,M1,S1,SH2,M2,S2,SH3,M3,S3,SH4,M4,S4 should be generated)
pdf("probabilityFecalMiceTransferALLWeeks.pdf",width = 11,height = 11)
plot_grid(SH1,M1,S1,SH2,M2,S2,SH3,M3,S3,SH4,M4,S4,ncol=3,nrow=4,label_size = 10,scale = 0.8)
dev.off()


