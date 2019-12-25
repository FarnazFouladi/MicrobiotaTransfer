#Farnaz Fouladi
#12-02-2019
#Comparing shared SVs between ANT1 and ANT2

rm(list =ls())

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)
source("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/Rcode/pie2.R")
data<-read.table(paste0(output,"svMeta.txt"),sep = "\t",header=TRUE,check.names=FALSE,na.strings = "NA" ,comment.char="")
finishAbundanceIndex<-which(colnames(data)=="Sample")-1
#Removing low abundant sv
sv<-data[,1:finishAbundanceIndex]
sv<-sv[,colMeans(sv>0)>=0.1]
myT<-cbind(sv,data[,(finishAbundanceIndex+1):ncol(data)])
finishAbundanceIndex<-which(colnames(myT)=="Sample")-1
Donors<-as.character(myT$Donor[myT$Sample.type=="Human.donor"])
myTH<-myT[myT$Sample.type=="Fecal.slurry",]
myTMouse<-myT[myT$Sample.type=="Mouse.feces",]
#Selecting one patient
patients=c("AN81T1","AN81T2","AN40T1","AN40T2","AN703T1","AN703T2","AN34T1","AN34T2")
for (patient in patients){
myTH1<-myTH[myTH$Donor==patient,]
myTMouse1<-myTMouse[myTMouse$Donor==patient,]

Slurries<-levels(factor(myTH1$Slurry.ID1))
Slurries<-Slurries[Slurries %in% myTMouse1$Slurry.ID1]
bugs<-c("TotalbugSlurry","bugShared","bugSlurry","bugMouse","bugNogroup")
df<-data.frame(bugs)


for (week in 1:4){
  myTM<-myTMouse1[myTMouse1$Week==week,]
  TotalbugSlurry<-0
  bugShared<-0
  bugSlurry<-0
  bugMouse<-0
  bugNogroup<-0
  numMouse<-0
  
  for (slurry in Slurries){
    myTHS<-myTH1[myTH1$Slurry.ID1==slurry,]
    myTMS<-myTM[myTM$Slurry.ID1==slurry,]
    for (j in 1:nrow(myTMS)){
      #print(nrow(myTMS))
      numMouse<-numMouse+1
      bugM<-myTMS[j,]
      for (i in 1:finishAbundanceIndex){
        bugH<-myTHS[,i]
        bugM_individual<-bugM[,i]
        if (bugH>0){
          TotalbugSlurry<-TotalbugSlurry+1 
          if (bugM_individual>0){
            bugShared<-bugShared+1
          }else{
            bugSlurry<-bugSlurry+1
          }
        }else {
          if (bugM_individual>0){
            bugMouse<-bugMouse+1
          }else{
            bugNogroup<-bugNogroup+1
          }
        }
      }
      
      value<-c(TotalbugSlurry,bugShared,bugSlurry,bugMouse,bugNogroup)
      df<-cbind(df,value)
      colnames(df)[which(colnames(df)=="value")]<-paste0(slurry,"week",week,"_",bugM$Sample.ID) #dim(df)5 536
      bugShared<-0
      bugSlurry<-0
      bugMouse<-0
      bugNogroup<-0
      TotalbugSlurry<-0
    }
  }
  print(paste0("Num of mice at week ",week," ",numMouse))
}
#How much percentage of total slurry bug is shared with mouse and how much not.
totalBuginSlurry<-as.numeric(df[1,2:ncol(df)])
twoGroups<-sweep(df[2:3,2:ncol(df)],2,totalBuginSlurry,'/')*100
# Calculating the percentage of bugShared, bugSlurry,bugMouse
sum<-colSums(df[2:4,2:ncol(df)])
threeGroups<-sweep(df[2:4,2:ncol(df)],2,sum,'/')*100
df1<-rbind(df[,2:ncol(df)],twoGroups,threeGroups)
rownames(df1)<-c("TotalbugSlurry","bugShared","bugSlurry","bugMouse","bugNogroup","percentBugShared","percentBugSlurry","percentBugShared3","percentBugSlurry3","percentBugMouse3")

names<-sapply(colnames(df1),function(x){strsplit(x,"week")[[1]][2]})
weeks<-sapply(names,function(x){strsplit(x,"_")[[1]][1]})

meanWeek1<-as.vector(rowMeans(df1[which(weeks=="1")]))
meanWeek2<-as.vector(rowMeans(df1[which(weeks=="2")]))
meanWeek3<-as.vector(rowMeans(df1[which(weeks=="3")]))
meanWeek4<-as.vector(rowMeans(df1[which(weeks=="4")]))
sd1<-apply(df1[which(weeks=="1")],1,sd)
sd2<-apply(df1[which(weeks=="2")],1,sd)
sd3<-apply(df1[which(weeks=="3")],1,sd)
sd4<-apply(df1[which(weeks=="4")],1,sd)
df2<-cbind(meanWeek1,meanWeek2,meanWeek3,meanWeek4,sd1,sd2,sd3,sd4,df1)
write.table(df2,paste0("bugSlurryMouse_",patient,".txt"),sep="\t")
}

#Comparing the shared SVs

patients=c("AN81T1","AN81T2","AN40T1","AN40T2","AN703T1","AN703T2","AN34T1","AN34T2")
df<-data.frame()

for (patient in patients){
  T134<-read.table(paste0(output,"bugSlurryMouse_",patient,".txt"),sep="\t",header=TRUE)
  T134_1<-T134[,9:ncol(T134)]
  T134_1<-as.data.frame(t(T134_1))
  T134_1$Sample<-rownames(T134_1)
  df<-rbind(df,T134_1)
  
}

week<-sapply(df$Sample, function(x)strsplit(x,"week.")[[1]][2])
df$Time<-sapply(week, function(x)strsplit(x,"_")[[1]][1])

df$Treatment<-sapply(df$Sample, function(x)strsplit(x,"[.]")[[1]][1])
patient<-sapply(df$Sample, function(x)strsplit(x,"_")[[1]][2])
df$Patient<-sapply(patient, function(x)strsplit(x,"AN")[[1]][1])
df$Patient<-factor(df$Patient,levels = c("34","40","81","703"))

theme_set(theme_classic(base_size = 10))
plot1<-ggplot(data=df,aes(x=Patient,y=percentBugShared3,color=Treatment))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape=16, size =0.4,position = position_jitterdodge(0.1),aes(color=Treatment))+
  labs(y="Shared (%)",x="Patients")+scale_x_discrete(labels=c("AN_5","AN_6","AN_7","AN_8"))+
  annotate(geom = "text",x=3,y=65,label=c("adjusted p<0.01"))

pval<-vector()
pval[1]<-t.test(df[df$Patient=="34" & df$Treatment=="T1",]$percentBugShared3,
       df[df$Patient=="34" & df$Treatment=="T2",]$percentBugShared3)$p.value

pval[2]<-t.test(df[df$Patient=="40" & df$Treatment=="T1",]$percentBugShared3,
       df[df$Patient=="40" & df$Treatment=="T2",]$percentBugShared3)$p.value

pval[3]<-t.test(df[df$Patient=="81" & df$Treatment=="T1",]$percentBugShared3,
       df[df$Patient=="81" & df$Treatment=="T2",]$percentBugShared3)$p.value

pval[4]<-t.test(df[df$Patient=="703" & df$Treatment=="T1",]$percentBugShared3,
       df[df$Patient=="703" & df$Treatment=="T2",]$percentBugShared3)$p.value

p.adjust(pval,method = "BH")

pdf("ComparisonOfT1andT2.pdf",height = 5,width = 5)
print(plot1)
dev.off()

