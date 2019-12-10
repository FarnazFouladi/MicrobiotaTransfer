#Farnaz Fouladi
#12-1-2019
#Comparison of percent of sequence variants that were shared, were detected only human
#fecal samples or only in slurries at different thresholds.

rm(list =ls())

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)
source("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/Rcode/pie2.R")
data<-read.table(paste0(output,"svMeta.txt"),sep = "\t",header=TRUE,check.names=FALSE,na.strings = "NA" ,comment.char="")
finishAbundanceIndex<-which(colnames(data)=="Sample")-1
#Removing low abundant sv
sv<-data[,1:finishAbundanceIndex]
#Select a threshold
threshold<-0.03
sv<-sv[,colMeans(sv>0)>=threshold]
myT<-cbind(sv,data[,(finishAbundanceIndex+1):ncol(data)])
finishAbundanceIndex<-which(colnames(myT)=="Sample")-1
Donors<-as.character(myT$Donor[myT$Sample.type=="Human.donor"])
myT1<-myT[myT$Sample.type=="Human.donor"|myT$Sample.type=="Fecal.slurry",]
Groups<-c("Total bug in Feces","Shared btw feces and slurry","only in feces","only in slurry","None")
df<-data.frame(Groups)
numSlurry<-0
for ( donor in Donors){
  fecal<-myT1[myT1$Sample.type=="Human.donor"& myT1$Donor==donor,]
  donorTable<-myT1[myT1$Sample.type=="Fecal.slurry"&myT1$Donor==donor,]
  slurries<-nrow(donorTable)
  
  for (s in 1:slurries)
    
  {
    numSlurry<-numSlurry+1
    TotalFecalbug<-0
    Shared<-0
    bugFecal<-0
    bugSlurry<-0
    noGroup<-0
    
    slurry<-donorTable[s,]
    
    for (i in 1:finishAbundanceIndex){
      
      fecalBug<-fecal[,i]
      slurryBug<-slurry[,i]
      
      if (fecalBug>0){
        TotalFecalbug<-TotalFecalbug+1
        if (slurryBug>0){
          Shared<-Shared+1
        }else{
          bugFecal<-bugFecal+1
        }
      } else{
        if (slurryBug>0){
          bugSlurry<-bugSlurry+1
        }else{
          noGroup<-noGroup+1
        }
      }
    }
    value<-c(TotalFecalbug,Shared,bugFecal,bugSlurry,noGroup)
    df<-cbind(df,value)
    colnames(df)[which(colnames(df)=="value")]<-paste0(donor,"_",slurry$Slurry.ID)
  }
}
print(paste0("Number of slurries is ",numSlurry)) #63
#How much percentage of total fecal bug is shared with slurry and how much not.
totalBuginFeces<-as.numeric(df[1,2:ncol(df)])
twoGroups<-sweep(df[2:3,2:ncol(df)],2,totalBuginFeces,'/')*100
# Calculating the percentage of bugShared, bugHuman,bugMouse
sum<-colSums(df[2:4,2:ncol(df)])
threeGroups<-sweep(df[2:4,2:ncol(df)],2,sum,'/')*100
df1<-rbind(df[,2:ncol(df)],twoGroups,threeGroups)
rownames(df1)<-c("TotalbugFecal","bugShared","bugFecal","bugSlurry","bugNogroup","sharedBugs","fecalBugs","percentBugShared","percentBugfecal","percentBugSlurry")
mean<-as.vector(rowMeans(df1))
sd<-apply(df1,1,sd)
df2<-cbind(mean,sd,df1)
write.table(df2,paste0("bugFecalSlurry_",threshold,".txt"),sep="\t")

#Pie Chart: paired analysis
meanPercentage<-c(as.numeric(format(df2[8,1],digits = 4)),as.numeric(format(df2[9,1],digits = 3)),as.numeric(format(df2[10,1],digits = 4)))
deviations<-c(as.numeric(format(df2[8,2],digits = 3)),as.numeric(format(df2[9,2],digits = 3)),as.numeric(format(df2[10,2],digits = 3)))
slices<-paste(meanPercentage,"\u00b1",deviations,"%",sep="")
labels<-c("Shared","Only in human fecal sample","Only in slurry")
labels<-paste(labels,"\n",slices)

png(paste0("PieChartHumanFecalSlurry_",threshold,".png"), units="in", width=13, height=13,res=300)
par(mfrow=c(2,1),mar=c(2,2,2,2))
pie2(meanPercentage,labels,col = c("indianred1","blue","orchid"),radius = 0.6,cex=1.2)
title("Paired analysis\n Average percentage of sequence variants\n present in a human fecal sample-slurry pair",line = -3,cex.main=1.1)

#Nonpaired analysis
sharedName<-vector()
slurryOnlyName<-vector()
fecalOnlyName<-vector()
noGroup<-vector()
sharedIndex<-1
slurryIndex<-1
fecalIndex<-1
noGroupIndex<-1

for (i in 1:finishAbundanceIndex){
  bugFecal<-myT[myT$Sample.type=="Human.donor",i]
  bugSlurry<-myT[myT$Sample.type=="Fecal.slurry" & myT$Donor %in% Donors,i]
  
  if (sum(bugFecal>0)>0 & sum(bugSlurry>0)>0){
    sharedName[sharedIndex]<-names(myT)[i] 
    sharedIndex<-sharedIndex+1
  }else if (sum(bugFecal>0)>0 & sum(bugSlurry>0)==0){
    fecalOnlyName[fecalIndex]<-names(myT)[i]
    fecalIndex<-fecalIndex+1
  }else if (sum(bugFecal>0)==0 & sum(bugSlurry>0)>0){
    slurryOnlyName[slurryIndex]<-names(myT)[i] 
    slurryIndex<-slurryIndex+1
  }else {
    noGroup[noGroupIndex]<-names(myT)[i] 
    noGroupIndex<-noGroupIndex+1
  }
}

#Pie chart from non-matched analysis
totalSeq<-length(sharedName)+length(slurryOnlyName)+length(fecalOnlyName)
shared<-length(sharedName)/totalSeq*100
onlyFecal<-length(fecalOnlyName)/totalSeq*100
onlySlurry<-length(slurryOnlyName)/totalSeq*100

slices<-c(as.numeric(format(shared,digits = 4)),as.numeric(format(onlyFecal,digit=2)),as.numeric(format(onlySlurry,digits = 4)))
labels<-c("Shared","Only in human fecal samples","Only in slurries")
labels<-paste(labels,slices)
labels<-paste(labels,"%",sep = "")

pie2(slices,labels,col = c("indianred1","blue","orchid"),radius = 0.6,cex=1.2)
title("Non-paired analysis\nPercentage of sequence variants\n present in human fecal samples and slurries",line = -3,cex.main=1.2)

dev.off()

slurryOnlyName[length(slurryOnlyName)+1:(length(sharedName)-length(slurryOnlyName))]<-NA
fecalOnlyName[length(fecalOnlyName)+1:(length(sharedName)-length(fecalOnlyName))]<-NA
noGroup[length(noGroup)+1:(length(sharedName)-length(noGroup))]<-NA
aframe<-data.frame(sharedName,fecalOnlyName,slurryOnlyName,noGroup)
write.table(aframe,paste0("BugNamesFecalSlurry",threshold,".txt"),sep="\t",row.names = FALSE)

#Comparing different thresholds 
t10<-read.table("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/bugFecalSlurry.txt",sep="\t",header = TRUE)
t7<-read.table("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/bugFecalSlurry_0.07.txt",sep="\t",header = TRUE)
t5<-read.table("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/bugFecalSlurry_0.05.txt",sep="\t",header = TRUE)
t3<-read.table("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/bugFecalSlurry_0.03.txt",sep="\t",header = TRUE)


df_all<-data.frame("10"=c(as.numeric(t10[8,3:ncol(t10)]),as.numeric(t10[9,3:ncol(t10)]),as.numeric(t10[10,3:ncol(t10)])),
                    "7"=c(as.numeric(t7[8,3:ncol(t7)]),as.numeric(t7[9,3:ncol(t7)]),as.numeric(t7[10,3:ncol(t7)])),
                    "5"=c(as.numeric(t5[8,3:ncol(t5)]),as.numeric(t5[9,3:ncol(t5)]),as.numeric(t5[10,3:ncol(t5)])),
                    "3"=c(as.numeric(t3[8,3:ncol(t3)]),as.numeric(t3[9,3:ncol(t3)]),as.numeric(t3[10,3:ncol(t3)])),
                    Group=c(rep("Shared",ncol(t10)-2),rep("OnlyInFecal",ncol(t10)-2),rep("OnlyInSlurry",ncol(t10)-2)))
 


df_all1<-data.frame(percent=c(df_all$X3,df_all$X5,df_all$X7,df_all$X10),
                    Filter=c(rep("3",nrow(df_all)),rep("5",nrow(df_all)),rep("7",nrow(df_all)),rep("10",nrow(df_all))),
                             Group=rep(df_all$Group,4))
                                                                                                
df_all1$Filter<-factor(df_all1$Filter,levels = c("3","5","7","10"))

df_fecalSlurry<-df_all1

theme_set(theme_classic(base_size = 12))
plot1<-ggplot(df_fecalSlurry,aes(x=Group,y=percent,color=Filter))+geom_boxplot(outlier.shape = NA)+
  labs(y="Percent (%)",x="")+geom_jitter(shape=16, size =0.4,position = position_jitterdodge(0.2),aes(color=Filter))+
  ylim(0,100)+
  scale_x_discrete(labels=c("Only in human \nfecal sample","Only in slurry","Shared"))


fecal<-df_all1[df_all1$Group=="OnlyInFecal",]
slurry<-df_all1[df_all1$Group=="OnlyInSlurry",]
shared<-df_all1[df_all1$Group=="Shared",]

TukeyHSD(aov(lm(fecal$percent ~ fecal$Filter)))
TukeyHSD(aov(lm(slurry$percent ~ slurry$Filter)))
TukeyHSD(aov(lm(shared$percent ~ shared$Filter)))



                          
