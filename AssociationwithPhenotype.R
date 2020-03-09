#Farnaz Fouladi
#08-18-2019

rm(list =ls())

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

#Removing Relative mesenteric fat and mesenteric fat
myT<-myT[,-c(which(colnames(myT)=="Mes.fat.wt"),which(colnames(myT)=="Rel.mes.fat.wt"),
             which(colnames(myT)=="Brown.fat.wt"),which(colnames(myT)=="Rel.brown.fat.wt"),
             which(colnames(myT)=="Lean.mass"),which(colnames(myT)=="Lean.mass.pct.change"))]
myT1<-myT[myT$Sample.type=="Mouse.feces",]

spearman<-vector()
rho<-vector()
bugName<-vector()
variableName<-vector()
time<-vector()
index<-1
pdf("AssociationwithPhenotype.pdf")
par(mfrow=c(3,3))

for (i in 1:finishAbundanceIndex){
  
  bug1<-myT1[,i]
  
  if(  sum(bug1!=0) > nrow(myT1)/10) {
    
    for (week in 1:4){
      
      myTweek<-myT1[myT1$Week==week,]
      bug<-myTweek[,i]
      
      for (j in which(colnames(myT1)=="Body.wt"):which(colnames(myT1)=="Cum.food.consumption"))
      {
        variable<-myTweek[,j]
        spearman[index]<-cor.test(variable,bug,method="spearman")$p.value
        rho[index]<-cor.test(variable,bug,method="spearman")$estimate
        variableName[index]<-colnames(myTweek)[j]
        time[index]<-week
        bugName[index]<-colnames(myTweek)[i]
        
        text<-paste0("p-value= ",format(spearman[index],digits = 3),"\nrho= ",format(rho[index],digits = 3),"\nweek= ",time[index])
        plot(variable,bug,ylab= bugName[index],xlab = variableName[index],main = text,cex.main=0.75)
        index<-index+1
      }
      
      if (week==4){
        for (j in which(colnames(myT)=="Cecum.wt"):which(colnames(myT)=="Rel.gonadal.fat.wt"))
        {
          variable<-myTweek[,j]
          spearman[index]<-cor.test(variable,bug,method="spearman")$p.value
          rho[index]<-cor.test(variable,bug,method="spearman")$estimate
          variableName[index]<-colnames(myTweek)[j]
          time[index]<-week
          bugName[index]<-colnames(myTweek)[i]
          
          text<-paste0("p-value= ",format(spearman[index],digits = 3),"\nrho= ",format(rho[index],digits = 3),"\nweek= ",time[index])
          plot(variable,bug,ylab= bugName[index],xlab = variableName[index],main = text,cex.main=0.75)
          index<-index+1
        }
      }
    }
  }
}
dev.off()
aframe<-data.frame(bugName,time,variableName,spearman,rho)
aframe<-aframe[order(aframe$spearman),]
aframe$adjustedSpearman<-p.adjust(aframe$spearman,method = "BH")
write.table(aframe,"interactionPhenotypeBug.text",sep="\t",row.names=FALSE)
