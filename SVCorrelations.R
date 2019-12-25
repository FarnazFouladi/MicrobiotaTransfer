#Farnaz Fouladi
#08-18-2019

#This R code generates Supplementary Figure 6. 

rm(list = ls())

library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)

output<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/output/"
input<-"/Users/farnazfouladi/Google Drive/MicrobiotaTransfer/input/"
setwd(output)

myT<-read.table(paste0(output,"interactionPhenotypeBug.text"),header = TRUE,sep = "\t")
abundance<-read.table(paste0(output,"svMeta.txt"),sep = "\t",header=TRUE,check.names=FALSE,na.strings = "NA" ,comment.char="")
finishAbundanceIndex<-which(colnames(abundance)=="Sample")-1
# Removing low abundant taxa 
sv<-abundance[,1:finishAbundanceIndex]
sv<-sv[,colMeans(sv>0)>=0.1]
abundance<-cbind(sv,abundance[,(finishAbundanceIndex+1):ncol(abundance)])
finishAbundanceIndex<-which(colnames(abundance)=="Sample")-1
abundance1<-abundance[abundance$Sample.type=="Mouse.feces",]
# Attaching taxanomy
taxa<-read.table(paste0(input,"taxForwardReads.txt"),header = TRUE,sep = "\t")
rownames(taxa)<-sapply(seq(1:nrow(taxa)),function(x){paste0("SV",x)})
taxa<-cbind(rownames(taxa),taxa)
myT1<-merge(myT,taxa,by.x="bugName",by.y="rownames(taxa)",sort=FALSE)

#Change the name of variables
myT1$variableName<-sapply(as.character(myT1$variableName),function(x){
  if (x=="Body.wt") return("Body weight")
  else if (x=="Body.wt.pct.change") return("Body weight percentage change")
  else if (x=="Brown.fat.wt") return("Brown fat weight")
  else if (x=="Cecum.wt") return("Cecum weight")
  else if (x=="Cum.food.consumption") return("Cumulative food intake")
  else if (x=="Fat.mass") return("Fat mass")
  else if (x=="Fat.mass.pct.change") return("Fat mass percentage change")
  else if (x=="Gonadal.fat.wt") return("Gonadal fat weight")
  else if (x=="Rel.brown.fat.wt") return("Relative brown fat weight")
  else if (x=="Rel.cecum.wt") return("Relative cecum weight")
  else if (x=="Rel.gonadal.fat.wt") return("Relative gonadal fat weight")
  else if (x=="Rel.SI.wt") return("Relative small intestinal weight")
  else if (x=="SI.wt") return("Small intestinal weight")
  else if (x=="Wk.food.consumption") return("Weekly food intake")
  else return(print("error"))
})


#Correlated SV with humans
myTCor<-read.table(paste0(output,"SlurryVsMouseCor.txt"),header = TRUE,sep = "\t")
orchid<-myTCor[myTCor$adjustedSpearman<0.05 & myTCor$adjustedSpearman>=0.01,]
blue<-myTCor[myTCor$adjustedSpearman<0.01 & myTCor$adjustedSpearman>=0.001,]
turquoise1<-myTCor[myTCor$adjustedSpearman<0.001,]

Variables<-levels(factor(myT1$variableName))
myList<-list()
pFisher<-vector()
index<-1

theme_set(theme_classic(base_size = 8))

for (week in 1:4){
  myTweek<-myT1[myT1$time==week,]
  orchidWeek<-orchid[orchid$time==week,]
  blueWeek<-blue[blue$time==week,]
  turquoise1Week<-turquoise1[turquoise1$time==week,]
  abundanceWeek<-abundance1[abundance1$Week==week,]
  
  #Ordering based on rank abundance using abundance table at each week
  sum<-colSums(abundanceWeek[,1:finishAbundanceIndex])
  names<-colnames(abundanceWeek[,1:finishAbundanceIndex])
  df<-data.frame(names,sum)
  df<-df[order(df$sum),]
  
  for (variable in Variables){
    
    if (variable %in% as.character(myTweek$variableName)){
      myTweekVar<-myTweek[myTweek$variableName==variable,]
      df1<-df[intersect(df$names,myTweekVar$bugName),]
      df1$rank<-rank(df1$sum)
      myTweekVar<-myTweekVar[match(df1$names,myTweekVar$bugName),]
      myTweekVar$SV<-df1$rank
      
      correlationWithHuman<-vector()
      
      for (i in 1:nrow(myTweekVar)){
        if (myTweekVar$bugName[i] %in% orchidWeek$bugnames )
          correlationWithHuman[i]<-"<0.05"
        else if (myTweekVar$bugName[i] %in% blueWeek$bugnames )  {
          correlationWithHuman[i]<-"<0.01"
        } else if (myTweekVar$bugName[i] %in% turquoise1Week$bugnames ){
          correlationWithHuman[i]<-"<0.001"
        } else {
          correlationWithHuman[i]<-"insignificant"
        }
      }
      
      myTweekVar$correlationWithHumanSV<-as.factor(correlationWithHuman)
      
      #Fisher test
      
      sigPhenotype<-myTweekVar[myTweekVar$adjustedSpearman<0.05,]
      sigPhenoSigCor<-sum(sigPhenotype$correlationWithHumanSV!="insignificant")
      sigPhenoNonSigCor<-sum(sigPhenotype$correlationWithHumanSV=="insignificant")
      
      NonSigPhenotype<-myTweekVar[!myTweekVar$adjustedSpearman<0.05,] 
      NonsigPhenoSigCor<-sum(NonSigPhenotype$correlationWithHumanSV!="insignificant")
      NonsigPhenoNonSigCor<-sum(NonSigPhenotype$correlationWithHumanSV=="insignificant")
      
      myMatrix<-matrix(c(sigPhenoSigCor,NonsigPhenoSigCor,sigPhenoNonSigCor,NonsigPhenoNonSigCor),nrow=2,dimnames = list(Phenotype = c("SigPheno", "nonSigPheno"),Taxa = c("SigTaxa", "NonSigTaxa")))
      
      P<-fisher.test(myMatrix)$p.value  
      e<-fisher.test(myMatrix)$estimate
      
      df_p<-as.data.frame(myMatrix)
      
      fisherTest<-c(P,e)
      Names<-c("pValue","oddsRatio")
      df_p1<-cbind(df_p,Names,fisherTest)
      
      write.table(df_p1,paste0("FishertestWeek",week,variable,".txt"),sep = "\t")
      
      pFisher[index]<-P
      
      vals<-c("insignificant"="red","<0.001"="turquoise1","<0.01"="blue","<0.05"="orchid")
      
      plot<-ggplot(myTweekVar,aes(x=SV,y=log10(adjustedSpearman),col=correlationWithHumanSV))+geom_point(size=1)+geom_hline(yintercept=log10(0.05), linetype="dashed", color="black", size=1)+
        labs(title=paste0(variable,"\nweek= ",week,"\nFisher p-value= ",format(P,digits = 2)),x="Rank abundance of SVs in mouse fecal pellets",y="Adjusted p-value")+
        scale_color_manual(values =vals)+ylim(c(-3,0))+theme(legend.position = "none")
      myList[[index]]<-plot
      index<-index+1
    }
  }
}

p.adjust(sort(pFisher),method = "BH")
#No significant p-value after adjustment

pdf("CorrelationswithPhenotype.pdf",width = 8,height = 8)
for (i in c(1,10,19)){
  grid.arrange(myList[[i]],myList[[i+1]],myList[[i+2]],myList[[i+3]],
            myList[[i+4]],myList[[i+5]],myList[[i+6]],myList[[i+7]],myList[[i+8]],nrow=3,ncol=3)
}
plot_grid(myList[[i+9]],myList[[i+10]],myList[[i+11]],nrow=3,ncol=3,scale = 0.9)
dev.off()

#png("SignificantCorrelationswithPhenotype.png", units="in", width=8, height=8,res=300)
pdf("SignificantCorrelationswithPhenotype.pdf",height = 7,width = 10)
plot_grid(myList[[11]],myList[[17]],myList[[24]],myList[[21]],
          myList[[24]],nrow=2,ncol=3,scale = 0.9)
dev.off()

#Legend Extract
plot<-ggplot(myTweekVar,aes(x=SV,y=log10(adjustedSpearman),col=correlationWithHumanSV))+geom_point(size=1)+geom_hline(yintercept=log10(0.05), linetype="dashed", color="black", size=1)+
  labs(title=paste0(variable,"\nweek= ",week,"\nFisher p-value= ",format(P,digits = 2)),x="Rank abundance of SVs in mice",y="Adjusted p-value")+ theme_classic(base_size = 9)+
  scale_color_manual(values =vals,labels=c("adj.p<0.001","0.001<adj.p<0.01","0.01<adj.p<0.05","adj.p>0.05"),name="Significance of correlations with human donor SVs")+ylim(c(-5,0))

#png("legendCorrelation.png", units="in", width=10, height=10,res=300)
pdf("legendCorrelation.pdf",height = 2, width = 3)
legend <- cowplot::get_legend(plot)
grid.newpage()
grid.draw(legend)
dev.off()
