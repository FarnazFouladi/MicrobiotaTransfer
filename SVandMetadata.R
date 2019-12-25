#Farnaz Fouladi
#08-18-2019

# This R code prepares the count table and the metadata for downstream analyses
# and produces MDS plots.

rm(list=ls())

library(ggplot2)
library(vegan)
library(ggsignif)
library(cowplot)

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)

data<-read.table(paste0(input,"ForwardReads.txt"),sep="\t",header = TRUE,comment.char = "")

#Sequencing depth
sumCountReads<-data.frame(rownames(data),rowSums(data))
sumCountReads<-sumCountReads[order(sumCountReads$rownames.data.),]

colnames(data)<-sapply(seq(1:ncol(data)), function(i){paste0("SV",i)})
data1=data

#Normalization
Sum<-rowSums(data1)
Average<-mean(Sum)
data1_relab<-sweep(data1,1,rowSums(data1),"/")
data2<-apply(data1_relab,2,function(x){log10(x*Average+1)})

#Uncomment this only if normalization is not required (To calculate Shannon Diversity Index)
#data2=data1

data2<-data2[order(rownames(data2)),]

#Adding metadata
meta<-read.table(paste0(input,"firstfour_metadata.txt"),check.names=FALSE, na.strings="n.a.", comment.char="", header=TRUE, sep="\t",row.names = 1 )
meta1<-meta[intersect(rownames(meta),rownames(data1)),]
meta1<-meta1[order(rownames(meta1)),]
data3<-data.frame(data2,Sample=rownames(meta1),meta1,SumReads=sumCountReads$rowSums.data.)

#Adding "Run" column
Run<-sapply(data3$Sample,function(x){strsplit(as.character(x),"_")[[1]][2]})

#Adding "Donor" group. 
Donor<-paste0(Run,as.character(data3$Mouse.group))
#Adding human donors(Human.donor or Fecal.clurry) to the column "Donor"
for(i in which(data3$Sample.type=="Human.donor")){
  print(paste0("row number ",i, " is ", data3$Sample.ID[i]," with barcode ",data3$Sample[i]))
}

Donor[6]<-"AN40HC"
Donor[8]<-"AN81T1"
Donor[289]<-"AN703HC"
Donor[292]<-"AN703T1"
Donor[295]<-"AN703T2"
Donor[339]<-"AN34T1"
Donor[341]<-"AN34T2"
Donor[343]<-"AN34HC"
Donor[353]<-"AN40T1"
Donor[355]<-"AN81T2"
Donor[399]<-"AN81HC"

for(i in c(which(data3$Sample.type=="Fecal.slurry"))){
  Donor[i]<-paste0(strsplit(as.character(data3$Sample[i]),"_")[[1]][2],strsplit(as.character(data3$Sample.ID[i]),"[.]")[[1]][1])
}

Donor<-sapply(Donor, function(x){if (x=="AN34NA"|x=="AN40NA"|x=="AN81NA"|x=="AN703NA"){x=NA}else{return(x)}})

#Creating the Group column 
Group<-vector()
index<-1
for (i in 1:length(Donor)){
  if (is.na(Donor[i])==FALSE){
    if(Donor[i]=="AN703T2"|Donor[i]=="AN703T1"|Donor[i]=="AN703HC"){
      Group[index]<-strsplit(Donor[i],"AN...")[[1]][2]
    }else {
      Group[index]<-strsplit(Donor[i],"AN..")[[1]][2]
    }
  }else{
    Group[index]<-NA
  }
  index<-index+1
}
data4<-cbind(data3,Run,Donor,Group)
data4$Slurry.ID1<-paste0(data4$Slurry.ID,data4$Run)

#Removing one sample that has only 18 reads, Only one sample has read counts below 1000.
data4<-data4[-which(data4$SumReads<1000),]
summary(data4$SumReads)
"Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
6135   33586   52450   59209   78169  212225 
sd 32557.47"

#Total Number of reads
sum(data4$SumReads)
mean(data4$SumReads)
sd(data4$SumReads)

write.table(data4,"svMeta.txt",sep="\t",row.names = FALSE)
#write.table(data4,"non-normalizedsvMeta.txt",sep="\t",row.names = FALSE)

#Comparing sequencing depth between human feces, slurries and mouse fecal pellets (non-normalized data)
data_sub<-data4[data4$Sample.type=="Human.donor"|data4$Sample.type=="Fecal.slurry"|data4$Sample.type=="Mouse.feces",]
data_sub$Sample_type<-sapply(data_sub$Sample.type, function(x){ifelse(x=="Human.donor",return("Human fecal samples"),
                                                                      (ifelse(x=="Fecal.slurry",return("Slurries"),
                                                                              return("Mouse fecal pellets"))))})
data_sub$Sample_type<-factor(data_sub$Sample_type,levels = c("Human fecal samples","Slurries","Mouse fecal pellets"))
png("SequencingDepth.png", units="in", width=5, height=5,res=300)
pdf("SequencingDepth.pdf")
theme_set(theme_classic(base_size = 12))
ggplot(data = data_sub, aes(x=Sample_type,y=SumReads))+geom_boxplot(fill="grey")+labs(y="Sequence depth")+
  theme(axis.title.x = element_blank())
result<-aov(lm(data_sub$SumReads~data_sub$Sample.type))
TukeyHSD(result)
dev.off()

#Diversity on non-normalized data:
dataN_sub1<-data_sub[is.na(data_sub$Week)|data_sub$Week==1,]
finishAbundanceIndex<-which(colnames(dataN_sub1)=="Sample")-1
shanonDiversity<-diversity(dataN_sub1[,1:finishAbundanceIndex],index = "shannon")
Evenvess<-shanonDiversity/log(specnumber(dataN_sub1[,1:finishAbundanceIndex]))
result<-aov(lm(shanonDiversity~dataN_sub1$Sample_type))
TukeyHSD(result)
"$`dataN_sub1$Sample_type`
diff        lwr         upr     p adj
Fecal.slurry-Human.donor  0.07842659 -0.1777915  0.33464472 0.7506546
Mouse.feces-Human.donor  -0.29126172 -0.5374247 -0.04509873 0.0156554
Mouse.feces-Fecal.slurry -0.36968832 -0.4847107 -0.25466598 0.0000000
"
result<-aov(lm(Evenvess~dataN_sub1$Sample_type))
TukeyHSD(result)
"$`dataN_sub1$Sample_type`
diff         lwr        upr     p adj
Fecal.slurry-Human.donor  0.016106511 -0.02581708 0.05803010 0.6368512
Mouse.feces-Human.donor   0.007996439 -0.03228189 0.04827477 0.8861957
Mouse.feces-Fecal.slurry -0.008110071 -0.02693056 0.01071042 0.5671541
"

theme_set(theme_classic(base_size = 10))
plot1<-ggplot(data=dataN_sub1,aes(x=Sample_type,y=shanonDiversity))+geom_boxplot(fill="grey")+
  geom_boxplot(fill="grey")+labs(x="Sample type",y="Shannon Diversity")+
  geom_signif(y_position = c(4.1,4.3),xmin=c(1,2),xmax=c(3,3),annotations=c("*","***"), tip_length=0,vjust = 0.5,textsize =5)

plot2<-ggplot(data=dataN_sub1,aes(x=Sample_type,y=Evenvess))+geom_boxplot(fill="grey")+
  geom_boxplot(fill="grey")+labs(x="Sample type",y="Evenness")

#Effect of time on diversity
dataN_sub2<-data_sub[data_sub$Sample.type=="Mouse.feces",]
shanonDiversity2<-diversity(dataN_sub2[,1:finishAbundanceIndex],index = "shannon")
plot3<-ggplot(data=dataN_sub2,aes(x=as.factor(Week),y=shanonDiversity2))+geom_boxplot(fill="grey")+labs(x="Time",y="Shannon Diversity")
Evenness2<-shanonDiversity2/log(specnumber(dataN_sub2[,1:finishAbundanceIndex]))
plot4<-ggplot(data=dataN_sub2,aes(x=as.factor(Week),y=Evenness2))+geom_boxplot(fill="grey")+labs(x="Time",y="Evenness")+
  geom_signif(y_position = 0.85,xmin=1,xmax=4,annotations="**", tip_length=0,vjust = 0.5,textsize =5)

png("Diversity.png", units="in", width=7, height=7,res=300)
plot_grid(plot1, plot2,plot3,plot4,labels = "AUTO", align = 'h',ncol=2,nrow=2,label_size = 12,scale = 0.9)
dev.off()

result<-aov(lm(shanonDiversity2~as.factor(dataN_sub2$Week)))
TukeyHSD(result)
"$`as.factor(dataN_sub2$Week)`
diff         lwr        upr     p adj
2-1  0.05190863 -0.05920458 0.16302184 0.6248274
3-1  0.03783841 -0.07216166 0.14783848 0.8120514
4-1  0.07275273 -0.03724734 0.18275280 0.3224978
3-2 -0.01407022 -0.12500476 0.09686432 0.9879467
4-2  0.02084410 -0.09009044 0.13177864 0.9626260
4-3  0.03491432 -0.07490528 0.14473391 0.8454719"
result<-aov(lm(Evenness2~as.factor(dataN_sub2$Week)))
TukeyHSD(result)
"$`as.factor(dataN_sub2$Week)`
diff          lwr        upr     p adj
2-1  0.013960440 -0.001517445 0.02943832 0.0938148
3-1  0.009628970 -0.005693857 0.02495180 0.3686689
4-1  0.019400784  0.004077957 0.03472361 0.0064101
3-2 -0.004331470 -0.019784467 0.01112153 0.8882738
4-2  0.005440344 -0.010012653 0.02089334 0.8011343
4-3  0.009771814 -0.005525873 0.02506950 0.3538129"

#MDS plot on normalised data
theme_set(theme_classic(base_size = 12.5))
mydata<-data4[data4$Sample.type!="Denver.human" & data4$Sample.type!="Positive.control",]
finishAbundanceIndex<-which(colnames(mydata)=="Sample")-1
myMDS<-capscale(mydata[,1:finishAbundanceIndex]~1,distance="bray")
percentVariance<-myMDS$CA$eig/sum(eigenvals(myMDS))*100
df<-data.frame(MDS1=myMDS$CA$u[,1],MDS2=myMDS$CA$u[,2],Sample_type=mydata$Sample.type,Donor=mydata$Donor)
col=c("red","blue","orchid","purple","black","lightblue","hotpink","cyan","pink","tan","darkgrey","gold")
df$Sample_type<-factor(df$Sample_type,levels = c("Human.donor","Fecal.slurry","Mouse.feces"))

png("pcoA_1.png", units="in", width=8, height=5,res=300)
plot1<-ggplot(data=df,aes(x=MDS1,y=MDS2))+geom_point(aes(col=Sample_type))+
  scale_colour_manual(values=col[1:length(levels(factor(df$Sample_type)))],labels=c("Human fecal samples","Slurries","Mouse fecal pellets"),name="Sample type")+
  labs(x=paste0("MDS1 (",format(percentVariance[1],digits = 4),"%)"),y=paste0("MDS2 (",format(percentVariance[2],digits = 4),"%)"))+lims(y=c(min(df$MDS1),max(df$MDS2)),x=c(min(df$MDS1),max(df$MDS2)))
dev.off()


df$Donor_newName<-sapply(as.character(df$Donor),function(x){
  if (substr(x,3,4)=="34" & substr(x,5,6)=="HC" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","1"))
  else if (substr(x,3,4)=="40" & substr(x,5,6)=="HC" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","2"))
  else if (substr(x,3,4)=="81" & substr(x,5,6)=="HC" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","3"))
  else if (substr(x,3,4)=="70" & substr(x,6,7)=="HC" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","4"))
  else if (substr(x,3,4)=="34" & substr(x,5,6)=="T1" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","5"))
  else if (substr(x,3,4)=="40" & substr(x,5,6)=="T1" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","6"))
  else if (substr(x,3,4)=="81" & substr(x,5,6)=="T1" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","7"))
  else if (substr(x,3,4)=="70" & substr(x,6,7)=="T1" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","8"))
  else if (substr(x,3,4)=="34" & substr(x,5,6)=="T2" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","5"))
  else if (substr(x,3,4)=="40" & substr(x,5,6)=="T2" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","6"))
  else if (substr(x,3,4)=="81" & substr(x,5,6)=="T2" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","7"))
  else if (substr(x,3,4)=="70" & substr(x,6,7)=="T2" ) return(paste0(substr(x,nchar(x)-1,nchar(x)),"_","8"))
  
})

png("pcoA_2.png", units="in", width=8, height=5,res=300)
plot2<-ggplot(data=df,aes(x=MDS1,y=MDS2))+geom_point(aes(col=factor(Donor_newName),shape=Sample_type))+
  scale_colour_manual(values=col[1:length(levels(factor(df$Donor_newName)))],name="Donors")+
  scale_shape_manual(values=c(15,16,17),labels=c("Human fecal samples","Slurries","Mouse fecal pellets"),name="Sample type")+
  labs(x="MDS1 (13.11%)",y="MDS2 (11.43%)")+lims(y=c(min(df$MDS1),max(df$MDS2)),x=c(min(df$MDS1),max(df$MDS2)))
dev.off()

adonis(mydata[,1:finishAbundanceIndex]~factor(mydata$Sample.type))

"
Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
mydata$Sample.type   2    16.065  8.0327  43.675 0.11369  0.001 ***
Residuals          681   125.249  0.1839         0.88631           
Total              683   141.314                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
"

#Figure 1 
pdf("pco1AND2.pdf",height = 6,width = 16,useDingbats=FALSE)
plot_grid(plot1,plot2,ncol=2,nrow=1,labels = c("A","B"),scale = 0.9)
dev.off()

#Mouse data and time
mydata<-data4[data4$Sample.type=="Mouse.feces",]
finishAbundanceIndex<-which(colnames(mydata)=="Sample")-1
myMDS<-capscale(mydata[,1:finishAbundanceIndex]~1,distance="bray")
percentVariance<-myMDS$CA$eig/sum(eigenvals(myMDS))*100
df<-data.frame(MDS1=myMDS$CA$u[,1],MDS2=myMDS$CA$u[,2],Donor=mydata$Donor,time=mydata$Week)
col=c("red","blue","orchid","orange")

pdf("pcoA_3.pdf",width = 7,height = 5)
ggplot(data=df,aes(x=MDS1,y=MDS2))+geom_point(aes(col=factor(time)))+
  scale_colour_manual(values=col[1:length(levels(factor(df$time)))],labels=c("Week 1","Week 2","Week 3","Week 4"),name="Time")+
  labs(x=paste0("MDS1 (",format(percentVariance[1],digits = 4),"%)"),y=paste0("MDS2 (",format(percentVariance[2],digits = 4),"%)"))
dev.off()

#How much variation in the mouse microbiome data is explained by fecal samples and slurry samples
adonis(mydata[,1:finishAbundanceIndex] ~ mydata$Donor)
"
Call:
adonis(formula = mydata[, 1:finishAbundanceIndex] ~ mydata$Donor) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
mydata$Donor  11    77.884  7.0804  125.06 0.69878  0.001 ***
Residuals    593    33.573  0.0566         0.30122           
Total        604   111.457                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
"

adonis(mydata[,1:finishAbundanceIndex] ~ factor(mydata$Slurry.ID1))
"
Call:
adonis(formula = mydata[, 1:finishAbundanceIndex] ~ factor(mydata$Slurry.ID1)) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2
factor(mydata$Slurry.ID1)  68    85.072 1.25106  25.415 0.76328
Residuals                 536    26.385 0.04923         0.23672
Total                     604   111.457                 1.00000
Pr(>F)    
factor(mydata$Slurry.ID1)  0.001 ***
Residuals                           
Total                               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

"