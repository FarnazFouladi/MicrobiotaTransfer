#Farnaz Fouladi
#12-1-2019
#Data was first rarefied to the lowest sequence depth and then percent of sequence variants 
#that were shared and were detected only in slurries or only in mouse fecal pellets was determined.

rm(list =ls())

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)
source("/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/Rcode/pie2.R")
data<-read.table(paste0(output,"non-normalizedsvMeta.txt"),sep = "\t",header=TRUE,check.names=FALSE,na.strings = "NA" ,comment.char="")
finishAbundanceIndex<-which(colnames(data)=="Sample")-1
#Rarefied
out<-rrarefy(data[,1:finishAbundanceIndex],
             min(rowSums(data[,1:finishAbundanceIndex])))
out1<-log10(out+1)

#Removing low abundant sv
sv<-out1
sv<-sv[,colMeans(sv>0)>=0.1]
myT<-cbind(sv,data[,(finishAbundanceIndex+1):ncol(data)])
finishAbundanceIndex<-which(colnames(myT)=="Sample")-1
Donors<-as.character(myT$Donor[myT$Sample.type=="Human.donor"])
myTH<-myT[myT$Sample.type=="Fecal.slurry",]
myTMouse<-myT[myT$Sample.type=="Mouse.feces",]
Slurries<-levels(factor(myTH$Slurry.ID1))
Slurries<-Slurries[Slurries!="T1.slurry.EAN40"]
bugs<-c("TotalbugFecal","bugShared","bugSlurry","bugMouse","bugNogroup")
df<-data.frame(bugs)


for (week in 1:4){
  myTM<-myTMouse[myTMouse$Week==week,]
  TotalbugSlurry<-0
  bugShared<-0
  bugSlurry<-0
  bugMouse<-0
  bugNogroup<-0
  numMouse<-0
  
  for (donor in Donors){
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
          TotalbugSlurry<-TotalbugFecal+1 
          if (bugM_individual>0){
            bugShared<-bugSlurry+1
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
      colnames(df)[which(colnames(df)=="value")]<-paste0(donor,"week ",week,"_",bugM$Sample.ID) 
      bugShared<-0
      bugSlurry<-0
      bugMouse<-0
      bugNogroup<-0
      TotalbugSlurry<-0
    }
  }
  print(paste0("Num of mice at week ",week," ",numMouse))
}
#How much percentage of total Slurry bug is shared with mouse and how much not.
totalBuginSlurry<-as.numeric(df[1,2:ncol(df)])
twoGroups<-sweep(df[2:3,2:ncol(df)],2,totalBuginSlurry,'/')*100
# Calculating the percentage of bugShared, bugFecal,bugMouse
sum<-colSums(df[2:4,2:ncol(df)])
threeGroups<-sweep(df[2:4,2:ncol(df)],2,sum,'/')*100
df1<-rbind(df[,2:ncol(df)],twoGroups,threeGroups)
rownames(df1)<-c("TotalbugSlurry","bugShared","bugSlurry","bugMouse","bugNogroup","percentBugShared","percentBugSlurry","percentBugShared3","percentBugSlurry3","percentBugMouse3")
meanWeek1<-as.vector(rowMeans(df1[1:136]))
meanWeek2<-as.vector(rowMeans(df1[137:267]))
meanWeek3<-as.vector(rowMeans(df1[268:403]))
meanWeek4<-as.vector(rowMeans(df1[404:539]))
sd1<-apply(df1[1:136],1,sd)
sd2<-apply(df1[137:267],1,sd)
sd3<-apply(df1[268:403],1,sd)
sd4<-apply(df1[404:539],1,sd)
df2<-cbind(meanWeek1,meanWeek2,meanWeek3,meanWeek4,sd1,sd2,sd3,sd4,df1)
write.table(df2,"bugFecalMouse_rarefied.txt",sep="\t")

#Pie Chart for all weeks: Paired analysis
pdf("pieFecalMousebug_rarefied.pdf")
for (i in 1:4)
{
  meanPercentage<-c(as.numeric(format(df2[8,i],digit=4)),as.numeric(format(df2[9,i],digit=4)),as.numeric(format(df2[10,i],digit=4)))
  deviations<-c(as.numeric(format(df2[8,i+4],digits = 3)),as.numeric(format(df2[9,i+4],digits = 3)),as.numeric(format(df2[10,i+4],digits = 3)))
  slices<-paste(meanPercentage,"\u00b1",deviations,"%",sep="")
  labels<-c("Shared","Only in human fecal samples","Only in mouse fecal pellets")
  labels<-paste(labels,"\n",slices)
  pie(meanPercentage,labels,col = c("indianred1","blue","orchid"),radius = 0.6,cex=1.2)
  title(paste0("Average percentage of sequence variants\n present in either human or mouse fecal samples or both","\nWeek= ",i),line=-3,cex=0.3)
}


#only at week1
#i=1
#meanPercentage<-c(as.numeric(format(df2[8,i],digit=4)),as.numeric(format(df2[9,i],digit=4)),as.numeric(format(df2[10,i],digit=4)))
#deviations<-c(as.numeric(format(df2[8,i+4],digits = 3)),as.numeric(format(df2[9,i+4],digits = 3)),as.numeric(format(df2[10,i+4],digits = 3)))
#slices<-paste(meanPercentage,"\u00b1",deviations,"%",sep="")
#labels<-c("Shared","Only in slurry","Only in mouse fecal pellet")
#labels<-paste(labels,"\n",slices)

#png("PieChartHumanSlurryMouseFeces_rarefied.png", units="in", width=13, height=13,res=300)
#par(mfrow=c(2,1),mar=c(2,2,2,2))
#pie2(meanPercentage,labels,col = c("indianred1","blue","orchid"),radius = 0.6,cex=1.2)
#title(paste0("Paired analysis\nAverage percentage of sequence variants\n present in a slurry-mouse fecal pellet pair","\nWeek= ",i),line = -3,cex.main=1.1)

#Nonpaired analysis
for (week in 1:4){
  sharedName<-vector()
  FecalOnlyName<-vector()
  mouseOnlyName<-vector()
  noGroup<-vector()
  sharedIndex<-1
  FecalIndex<-1
  mouseIndex<-1
  noGroupIndex<-1
  myTw<-myT[is.na(myT$Week)|myT$Week==week,]
  
  for (i in 1:finishAbundanceIndex){
    
    bugFecal<-myTw[myTw$Sample.type=="Human.donor" & myTw$Donor %in% Donors,i]
    bugMouse<-myTw[myTw$Sample.type=="Mouse.feces" & myTw$Donor %in% Donors,i]
    
    if (sum(bugFecal>0)>0 & sum(bugMouse>0)>0){
      sharedName[sharedIndex]<-names(myT)[i] 
      sharedIndex<-sharedIndex+1
    }else if (sum(bugFecal>0)>0 & sum(bugMouse>0)==0){
      FecalOnlyName[FecalIndex]<-names(myT)[i]
      FecalIndex<-FecalIndex+1
    }else if (sum(bugFecal>0)==0 & sum(bugMouse>0)>0){
      mouseOnlyName[mouseIndex]<-names(myT)[i] 
      mouseIndex<-mouseIndex+1
    }else {
      noGroup[noGroupIndex]<-names(myT)[i] 
      noGroupIndex<-noGroupIndex+1
    }
  }
  
  #Pie chart from non-matched analysis
  totalSeq<-length(sharedName)+length(FecalOnlyName)+length(mouseOnlyName)+length(noGroup)
  shared<-length(sharedName)/totalSeq*100
  onlyMouse<-length(mouseOnlyName)/totalSeq*100
  onlyFecal<-length(FecalOnlyName)/totalSeq*100
  
  slices<-c(as.numeric(format(shared,digits = 4)),as.numeric(format(onlyFecal,digit=2)),as.numeric(format(onlyMouse,digits = 3)))
  labels<-c("Shared","Only in human fecal samples","Only in mouse fecal pellets")
  labels<-paste(labels,slices)
  labels<-paste(labels,"%",sep = "")
  pie(slices,labels,col = c("indianred1","blue","orchid"),radius = 0.6,cex=1.2)
  title(paste0("Non-paired analysis\nPercentage of sequence variants\n present in human and mouse fecal samples \n week= ",week),line = -3,cex.main=1.2)
  
}


dev.off()
