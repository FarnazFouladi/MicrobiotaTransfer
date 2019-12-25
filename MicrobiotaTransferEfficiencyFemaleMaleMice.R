#Farnaz Fouladi
#11-22-19
#Comparing the transfer efficiency between male and female mice. 
#Comparison of percent SVs that were shared, dectected only in slurries, and only in mice between female and male mice.

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

#Choose female or male mice
sex="F"
myTMouse_S<-myTMouse[myTMouse$Sex==sex,]
Slurries<-levels(factor(myTH$Slurry.ID1))
#Removing slurries that were not used for colonization of mice
Slurries_S<-Slurries[Slurries %in% as.character(myTMouse_S$Slurry.ID1)]
bugs<-c("TotalbugSlurry","bugShared","bugSlurry","bugMouse","bugNogroup")
df<-data.frame(bugs)

for (week in 1:4){
  myTM<-myTMouse_S[myTMouse_S$Week==week,]
  TotalbugSlurry<-0
  bugShared<-0
  bugSlurry<-0
  bugMouse<-0
  bugNogroup<-0
  numMouse<-0
  
  for (slurry in Slurries_S){
    
    if (slurry %in% myTM$Slurry.ID1){
    
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

weeks<-sapply(colnames(df1),function(x){strsplit(x," ")[[1]][2]})
weeksN<-sapply(weeks,function(x){strsplit(x,"_")[[1]][1]})

meanWeek1<-as.vector(rowMeans(df1[1:which(weeksN=="1")[length(which(weeksN=="1"))]]))
meanWeek2<-as.vector(rowMeans(df1[(which(weeksN=="1")[length(which(weeksN=="1"))]+1):which(weeksN=="2")[length(which(weeksN=="2"))]]))
meanWeek3<-as.vector(rowMeans(df1[(which(weeksN=="2")[length(which(weeksN=="2"))]+1):which(weeksN=="3")[length(which(weeksN=="3"))]]))
meanWeek4<-as.vector(rowMeans(df1[(which(weeksN=="3")[length(which(weeksN=="3"))]+1):which(weeksN=="4")[length(which(weeksN=="4"))]]))

sd1<-apply(df1[1:which(weeksN=="1")[length(which(weeksN=="1"))]],1,sd)
sd2<-apply(df1[(which(weeksN=="1")[length(which(weeksN=="1"))]+1):which(weeksN=="2")[length(which(weeksN=="2"))]],1,sd)
sd3<-apply(df1[(which(weeksN=="2")[length(which(weeksN=="2"))]+1):which(weeksN=="3")[length(which(weeksN=="3"))]],1,sd)
sd4<-apply(df1[(which(weeksN=="3")[length(which(weeksN=="3"))]+1):which(weeksN=="4")[length(which(weeksN=="4"))]],1,sd)
df2<-cbind(meanWeek1,meanWeek2,meanWeek3,meanWeek4,sd1,sd2,sd3,sd4,df1)
write.table(df2,paste0("bugSlurryMouse",sex,"Mice.txt"),sep="\t")

#Pie Chart for all weeks: Paired analysis
pdf(paste0("pieSlurrybug",sex,"Mice.pdf"))
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

dev.off()

#Box plot for the shared:

df_Female<-read.table(paste0(output,"bugSlurryMouseFMice.txt"),sep="\t",header = TRUE)
df_F<-t(df_Female)
df_F<-as.data.frame(df_F)
df_F<-df_F[9:nrow(df_F),]
df_F$sex<-rep("Female",nrow(df_F))
splitNames<-sapply(rownames(df_F),function(x){strsplit(x,"week.")[[1]][2]})
week<-sapply(splitNames,function(x){strsplit(x,"_")[[1]][1]})
df_F$week<-week

df_male<-read.table(paste0(output,"bugSlurryMouseMMice.txt"),sep="\t",header = TRUE)
df_M<-t(df_male)
df_M<-as.data.frame(df_M)
df_M<-df_M[9:nrow(df_M),]
df_M$sex<-rep("Male",nrow(df_M))
splitNames<-sapply(rownames(df_M),function(x){strsplit(x,"week.")[[1]][2]})
week<-sapply(splitNames,function(x){strsplit(x,"_")[[1]][1]})
df_M$week<-week

df_all<-rbind(df_F,df_M)
df_all$sex<-as.factor(df_all$sex)

library(cowplot)
library(ggsignif)
theme_set(theme_classic(base_size = 9))
plot1<-ggplot(df_all,aes(x=week,y=percentBugShared3,color=as.factor(sex)))+geom_boxplot(outlier.shape = NA)+
  labs(y="Shared (%)",x="Week",color="")+geom_jitter(shape=16, size =0.4,position = position_jitterdodge(),aes(color=sex))+
  ylim(0,100)+annotate(geom = "text",x=2,y=70,label="adjusted p=0.09",size=2)

plot2<-ggplot(df_all,aes(x=week,y=percentBugSlurry3,color=as.factor(sex)))+geom_boxplot(outlier.shape = NA)+
  labs(y="Only in slurries (%)",x="Week",color="")+geom_jitter(shape=16, size =0.4,position = position_jitterdodge(),aes(color=sex))+
  ylim(0,100)+annotate(geom = "text",x=2,y=60,label="adjusted p=0.09",size=2)

plot3<-ggplot(df_all,aes(x=week,y=percentBugMouse3,color=as.factor(sex)))+geom_boxplot(outlier.shape = NA)+
  labs(y="Only in mouse fecal pellets (%)",x="Week",color="")+geom_jitter(shape=16, size =0.4,position = position_jitterdodge(),aes(color=sex))+
  ylim(0,100)

pdf("ComparisonOfMaleAndFemaleMice.pdf",width = 9,height = 4)
plot_grid(plot1,plot2,plot3,nrow=1,ncol=3,labels = c("A","B","C"),label_size =12)
dev.off()

df_ll1<-df_all[df_all$week==1,]
df_ll2<-df_all[df_all$week==2,]
df_ll3<-df_all[df_all$week==3,]
df_ll4<-df_all[df_all$week==4,]

pval<-vector()
pval[1]<-t.test(df_ll1[df_ll1$sex=="Female",]$percentBugShared3,df_ll1[df_ll1$sex=="Male",]$percentBugShared3)$p.value
pval[2]<-t.test(df_ll2[df_ll2$sex=="Female",]$percentBugShared3,df_ll2[df_ll2$sex=="Male",]$percentBugShared3)$p.value
pval[3]<-t.test(df_ll3[df_ll3$sex=="Female",]$percentBugShared3,df_ll3[df_ll3$sex=="Male",]$percentBugShared3)$p.value
pval[4]<-t.test(df_ll4[df_ll4$sex=="Female",]$percentBugShared3,df_ll4[df_ll4$sex=="Male",]$percentBugShared3)$p.value

pval[5]<-t.test(df_ll1[df_ll1$sex=="Female",]$percentBugSlurry3,df_ll1[df_ll1$sex=="Male",]$percentBugSlurry3)$p.value
pval[6]<-t.test(df_ll2[df_ll2$sex=="Female",]$percentBugSlurry3,df_ll2[df_ll2$sex=="Male",]$percentBugSlurry3)$p.value
pval[7]<-t.test(df_ll3[df_ll3$sex=="Female",]$percentBugSlurry3,df_ll3[df_ll3$sex=="Male",]$percentBugSlurry3)$p.value
pval[8]<-t.test(df_ll4[df_ll4$sex=="Female",]$percentBugSlurry3,df_ll4[df_ll4$sex=="Male",]$percentBugSlurry3)$p.value

pval[9]<-t.test(df_ll1[df_ll1$sex=="Female",]$percentBugMouse3,df_ll1[df_ll1$sex=="Male",]$percentBugMouse3)$p.value
pval[10]<-t.test(df_ll2[df_ll2$sex=="Female",]$percentBugMouse3,df_ll2[df_ll2$sex=="Male",]$percentBugMouse3)$p.value
pval[11]<-t.test(df_ll3[df_ll3$sex=="Female",]$percentBugMouse3,df_ll3[df_ll3$sex=="Male",]$percentBugMouse3)$p.value
pval[12]<-t.test(df_ll4[df_ll4$sex=="Female",]$percentBugMouse3,df_ll4[df_ll4$sex=="Male",]$percentBugMouse3)$p.value

p.adjust(pval,method = "BH")
"
[1] 0.41456468 0.09316835 0.41456468 0.41456468 0.38359246 0.09316835 0.72415714
 [8] 0.57123737 0.72415714 0.72415714 0.72415714 0.90984658"




