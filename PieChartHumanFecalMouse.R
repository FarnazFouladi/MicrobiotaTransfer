#Farnaz Fouladi
#08-18-2019
#This R code generates the pie charts for SVs shared between human and mice fecal pellets. 

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
myTH<-myT[myT$Sample.type=="Human.donor",]
myTMouse<-myT[myT$Sample.type=="Mouse.feces",]
bugs<-c("TotalbugFecal","bugShared","bugFecal","bugMouse","bugNogroup")
df<-data.frame(bugs)

for (week in 1:4){
  myTM<-myTMouse[myTMouse$Week==week,]
  TotalbugFecal<-0
  bugShared<-0
  bugFecal<-0
  bugMouse<-0
  bugNogroup<-0
  numMouse<-0
  
  for (donor in Donors){
    myTHS<-myTH[myTH$Donor==donor,]
    myTMS<-myTM[myTM$Donor==donor,]
    for (j in 1:nrow(myTMS)){
      #print(nrow(myTMS))
      numMouse<-numMouse+1
      bugM<-myTMS[j,]
      for (i in 1:finishAbundanceIndex){
        bugH<-myTHS[,i]
        bugM_individual<-bugM[,i]
        if (bugH>0){
          TotalbugFecal<-TotalbugFecal+1 
          if (bugM_individual>0){
            bugShared<-bugShared+1
          }else{
            bugFecal<-bugFecal+1
          }
        }else {
          if (bugM_individual>0){
            bugMouse<-bugMouse+1
          }else{
            bugNogroup<-bugNogroup+1
          }
        }
      }
      
      value<-c(TotalbugFecal,bugShared,bugFecal,bugMouse,bugNogroup)
      df<-cbind(df,value)
      colnames(df)[which(colnames(df)=="value")]<-paste0(donor,"week ",week,"_",bugM$Sample.ID) 
      bugShared<-0
      bugFecal<-0
      bugMouse<-0
      bugNogroup<-0
      TotalbugFecal<-0
    }
  }
  print(paste0("Num of mice at week ",week," ",numMouse))
}
#How much percentage of total fecal bug is shared with mouse and how much not.
totalBuginFecal<-as.numeric(df[1,2:ncol(df)])
twoGroups<-sweep(df[2:3,2:ncol(df)],2,totalBuginFecal,'/')*100
# Calculating the percentage of bugShared, bugFecal,bugMouse
sum<-colSums(df[2:4,2:ncol(df)])
threeGroups<-sweep(df[2:4,2:ncol(df)],2,sum,'/')*100
df1<-rbind(df[,2:ncol(df)],twoGroups,threeGroups)
rownames(df1)<-c("TotalbugFecal","bugShared","bugFecal","bugMouse","bugNogroup","percentBugShared","percentBugFecal","percentBugShared3","percentBugFecal3","percentBugMouse3")
meanWeek1<-as.vector(rowMeans(df1[1:136]))
meanWeek2<-as.vector(rowMeans(df1[137:267]))
meanWeek3<-as.vector(rowMeans(df1[268:403]))
meanWeek4<-as.vector(rowMeans(df1[404:539]))
sd1<-apply(df1[1:136],1,sd)
sd2<-apply(df1[137:267],1,sd)
sd3<-apply(df1[268:403],1,sd)
sd4<-apply(df1[404:539],1,sd)
df2<-cbind(meanWeek1,meanWeek2,meanWeek3,meanWeek4,sd1,sd2,sd3,sd4,df1)
write.table(df2,"bugFecalMouse.txt",sep="\t")

#Pie Chart for all weeks: Paired analysis
pdf("pieFecalMousebug_suppFigure.pdf",height = 10, width = 10)
par(mfrow=c(2,2))
for (i in 1:4)
{
  meanPercentage<-c(as.numeric(format(df2[8,i],digit=4)),as.numeric(format(df2[9,i],digit=4)),as.numeric(format(df2[10,i],digit=4)))
  deviations<-c(as.numeric(format(df2[8,i+4],digits = 3)),as.numeric(format(df2[9,i+4],digits = 3)),as.numeric(format(df2[10,i+4],digits = 3)))
  slices<-paste(meanPercentage,"\u00b1",deviations,"%",sep="")
  labels<-c("Shared","Only in human fecal sample","Only in mouse fecal pellet")
  labels<-paste(labels,"\n",slices)
  pie2(meanPercentage,labels,col = c("indianred1","blue","orchid"),radius = 0.6,cex=0.9)
  title(paste0("Week= ",i),line=0,cex=0.1)
}

dev.off()
