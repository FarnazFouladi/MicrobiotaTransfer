#Farnaz Fouladi
#12-1-2019
#Comparison of percent of sequence variants that were shared, were detected only in
#slurries or only in mouse fecal pellets at different thresholds.

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
threshold<-0.07
sv<-sv[,colMeans(sv>0)>=threshold]
myT<-cbind(sv,data[,(finishAbundanceIndex+1):ncol(data)])
finishAbundanceIndex<-which(colnames(myT)=="Sample")-1
Donors<-as.character(myT$Donor[myT$Sample.type=="Human.donor"])
myTH<-myT[myT$Sample.type=="Fecal.slurry",]
myTMouse<-myT[myT$Sample.type=="Mouse.feces",]
Slurries<-levels(factor(myTH$Slurry.ID1))
Slurries<-Slurries[Slurries!="T1.slurry.EAN40"]
bugs<-c("TotalbugSlurry","bugShared","bugSlurry","bugMouse","bugNogroup")
df<-data.frame(bugs)


for (week in 1:4){
  myTM<-myTMouse[myTMouse$Week==week,]
  TotalbugSlurry<-0
  bugShared<-0
  bugSlurry<-0
  bugMouse<-0
  bugNogroup<-0
  numMouse<-0
  
  for (slurry in Slurries){
    myTHS<-myTH[myTH$Slurry.ID1==slurry,]
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
      colnames(df)[which(colnames(df)=="value")]<-paste0(slurry,"week ",week,"_",bugM$Sample.ID) #dim(df)5 536
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
meanWeek1<-as.vector(rowMeans(df1[1:135]))
meanWeek2<-as.vector(rowMeans(df1[136:265]))
meanWeek3<-as.vector(rowMeans(df1[266:400]))
meanWeek4<-as.vector(rowMeans(df1[401:535]))
sd1<-apply(df1[1:135],1,sd)
sd2<-apply(df1[136:265],1,sd)
sd3<-apply(df1[266:400],1,sd)
sd4<-apply(df1[401:535],1,sd)
df2<-cbind(meanWeek1,meanWeek2,meanWeek3,meanWeek4,sd1,sd2,sd3,sd4,df1)
write.table(df2,paste0("bugSlurryMouse_",threshold,".txt"),sep="\t")

#Pie Chart for all weeks: Paired analysis
pdf(paste0("pieSlurrybug",threshold,".pdf"))
for (i in 1:4)
{
  meanPercentage<-c(as.numeric(format(df2[8,i],digit=4)),as.numeric(format(df2[9,i],digit=4)),as.numeric(format(df2[10,i],digit=4)))
  deviations<-c(as.numeric(format(df2[8,i+4],digits = 3)),as.numeric(format(df2[9,i+4],digits = 3)),as.numeric(format(df2[10,i+4],digits = 3)))
  slices<-paste(meanPercentage,"\u00b1",deviations,"%",sep="")
  labels<-c("Shared","Only in Slurry","Only in Mouse")
  labels<-paste(labels,"\n",slices)
  pie(meanPercentage,labels,col = c("indianred1","blue","orchid"),radius = 0.4,cex=0.8)
  title(paste0("Average percentage of sequence variants\n present in either mouse or slurry or both","\nWeek= ",i),line=-3,cex=0.3)
}


#Nonpaired analysis
for (week in 1:4){
  sharedName<-vector()
  slurryOnlyName<-vector()
  mouseOnlyName<-vector()
  noGroup<-vector()
  sharedIndex<-1
  slurryIndex<-1
  mouseIndex<-1
  noGroupIndex<-1
  myTw<-myT[is.na(myT$Week)|myT$Week==week,]
  
  for (i in 1:finishAbundanceIndex){
    
    bugSlurry<-myTw[myTw$Sample.type=="Fecal.slurry" & myTw$Slurry.ID1 %in% Slurries,i]
    bugMouse<-myTw[myTw$Sample.type=="Mouse.feces" & myTw$Slurry.ID1 %in% Slurries,i]
    
    if (sum(bugSlurry>0)>0 & sum(bugMouse>0)>0){
      sharedName[sharedIndex]<-names(myT)[i] 
      sharedIndex<-sharedIndex+1
    }else if (sum(bugSlurry>0)>0 & sum(bugMouse>0)==0){
      slurryOnlyName[slurryIndex]<-names(myT)[i]
      slurryIndex<-slurryIndex+1
    }else if (sum(bugSlurry>0)==0 & sum(bugMouse>0)>0){
      mouseOnlyName[mouseIndex]<-names(myT)[i] 
      mouseIndex<-mouseIndex+1
    }else {
      noGroup[noGroupIndex]<-names(myT)[i] 
      noGroupIndex<-noGroupIndex+1
    }
  }
  
  #Pie chart from non-matched analysis
  totalSeq<-length(sharedName)+length(slurryOnlyName)+length(mouseOnlyName)+length(noGroup)
  shared<-length(sharedName)/totalSeq*100
  onlyMouse<-length(mouseOnlyName)/totalSeq*100
  onlySlurry<-length(slurryOnlyName)/totalSeq*100
  
  slices<-c(as.numeric(format(shared,digits = 4)),as.numeric(format(onlySlurry,digit=2)),as.numeric(format(onlyMouse,digits = 3)))
  labels<-c("Shared","Only in slurries","Only in mouse fecal pellets")
  labels<-paste(labels,slices)
  labels<-paste(labels,"%",sep = "")
  pie(slices,labels,col = c("indianred1","blue","orchid"),radius = 0.4,cex=0.8)
  title(paste0("Non-paired analysis\nPercentage of sequence variants\n present in slurries and mouse fecal pellets \n week= ",week),line = -3,cex.main=1.2)
  
  slurryOnlyName[length(slurryOnlyName)+1:(length(sharedName)-length(slurryOnlyName))]<-NA
  mouseOnlyName[length(mouseOnlyName)+1:(length(sharedName)-length(mouseOnlyName))]<-NA
  noGroup[length(noGroup)+1:(length(sharedName)-length(noGroup))]<-NA
  aframe<-data.frame(sharedName,slurryOnlyName,mouseOnlyName,noGroup)
  write.table(aframe,paste0("BugNamesSlurryMouse_week",week,"_",threshold,".txt"),sep="\t",row.names = FALSE)
}


dev.off()

#Comapring Filters:

#Comparing different thresholds 
t10<-read.table("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/bugSlurryMouse.txt",sep="\t",header = TRUE)
t7<-read.table("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/bugSlurryMouse_0.07.txt",sep="\t",header = TRUE)
t5<-read.table("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/bugSlurryMouse_0.05.txt",sep="\t",header = TRUE)
t3<-read.table("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/bugSlurryMouse_0.03.txt",sep="\t",header = TRUE)


df_all<-data.frame("10"=c(as.numeric(t10[8,3:ncol(t10)]),as.numeric(t10[9,3:ncol(t10)]),as.numeric(t10[10,3:ncol(t10)])),
                   "7"=c(as.numeric(t7[8,3:ncol(t7)]),as.numeric(t7[9,3:ncol(t7)]),as.numeric(t7[10,3:ncol(t7)])),
                   "5"=c(as.numeric(t5[8,3:ncol(t5)]),as.numeric(t5[9,3:ncol(t5)]),as.numeric(t5[10,3:ncol(t5)])),
                   "3"=c(as.numeric(t3[8,3:ncol(t3)]),as.numeric(t3[9,3:ncol(t3)]),as.numeric(t3[10,3:ncol(t3)])),
                   Group=c(rep("Shared",ncol(t10)-2),rep("OnlyInSlurry",ncol(t10)-2),rep("OnlyInMouse",ncol(t10)-2)))



df_all1<-data.frame(percent=c(df_all$X3,df_all$X5,df_all$X7,df_all$X10),
                    Filter=c(rep("3",nrow(df_all)),rep("5",nrow(df_all)),rep("7",nrow(df_all)),rep("10",nrow(df_all))),
                    Group=rep(df_all$Group,4))

df_all1$Filter<-factor(df_all1$Filter,levels = c("3","5","7","10"))

plot2<-ggplot(df_all1,aes(x=Group,y=percent,color=Filter))+geom_boxplot(outlier.shape = NA)+
  labs(y="Percent (%)",x="")+geom_jitter(shape=16, size =0.4,position = position_jitterdodge(0.2),aes(color=Filter))+
  scale_x_discrete(labels=c("Only in mouse \nfecal pellet","Only in slurry","Shared"))+
  ylim(0,100)+
  annotate(geom="text",x=c(1,2,3),y=c(90,90,90),label=c("adjusted p<0.05","adjusted p<0.001","adjusted p<0.001"),
           size=3)


mouse<-df_all1[df_all1$Group=="OnlyInMouse",]
slurry<-df_all1[df_all1$Group=="OnlyInSlurry",]
shared<-df_all1[df_all1$Group=="Shared",]

anova(lm(mouse$percent ~ mouse$Filter))
anova(lm(slurry$percent ~ slurry$Filter))
anova(lm(shared$percent ~ shared$Filter))

TukeyHSD(aov(lm(mouse$percent ~ mouse$Filter)))
TukeyHSD(aov(lm(slurry$percent ~ slurry$Filter)))
TukeyHSD(aov(lm(shared$percent ~ shared$Filter)))

#Plot1 from PieChartHumanFecalSlurry_threshold
pdf("TransferEfficiencyFilters.pdf",height = 6,width = 12)
plot_grid(plot1,plot2,labels = c("A","B"),ncol=2,nrow=1,scale = 0.9)
dev.off()


