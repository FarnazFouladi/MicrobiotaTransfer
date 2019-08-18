#Farnaz Fouladi
#08-18-2019

#This R code generates the pie charts in Figure 4. 

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
write.table(df2,"bugSlurryMouse.txt",sep="\t")

#Pie Chart for all weeks: Paired analysis
pdf("pieSlurrybug.pdf")
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


#For publication only at week1
i=1
meanPercentage<-c(as.numeric(format(df2[8,i],digit=4)),as.numeric(format(df2[9,i],digit=4)),as.numeric(format(df2[10,i],digit=4)))
deviations<-c(as.numeric(format(df2[8,i+4],digits = 3)),as.numeric(format(df2[9,i+4],digits = 3)),as.numeric(format(df2[10,i+4],digits = 3)))
slices<-paste(meanPercentage,"\u00b1",deviations,"%",sep="")
labels<-c("Shared","Only in slurry","Only in mouse fecal pellet")
labels<-paste(labels,"\n",slices)

png("PieChartHumanSlurryMouseFeces.png", units="in", width=13, height=13,res=300)
par(mfrow=c(2,1),mar=c(2,2,2,2))
pie2(meanPercentage,labels,col = c("indianred1","blue","orchid"),radius = 0.6,cex=1.2)
title(paste0("Paired analysis\nAverage percentage of sequence variants\n present in a slurry-mouse fecal pellet pair","\nWeek= ",i),line = -3,cex.main=1.1)

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
  write.table(aframe,paste0("BugNamesSlurryMouse_week",week,".txt"),sep="\t",row.names = FALSE)
}


dev.off()

#For publication
#week=1
#pie2(slices,labels,col = c("indianred1","blue","orchid"),radius = 0.6,cex=1.2)
#title(paste0("Non-paired analysis\nPercentage of sequence variants\n present in either any human slurries or mouse feces or both","\nWeek= ",week),line = -3,cex.main=1.2)
#dev.off()
