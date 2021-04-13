require(metagenomeSeq)
require(here)
require(reshape)
library(plyr)
library(ggplot2)
library(ggpubr)



presentThresh=1
fdepthThresh=100
#fdepthThresh

outputFolder <- paste0("outputBoston_Filt",fdepthThresh,"_SampPres", presentThresh)
print(outputFolder)

source(here:::here("scripts","00.1_CAMDA_jc_Boston.r"))
source(here:::here("scripts","00.2_CAMDA_jc_Boston_combine_clean.r"))

rm(list=setdiff(ls(), c( "fdepthThreshs","presentThreshs", "outputFolder")))

fdepthThresh<-gsub("Filt","",strsplit(outputFolder, "_")[[1]][2])
presentThresh<-gsub("SampPres","",strsplit(outputFolder, "_")[[1]][3])
print(outputFolder)
print(presentThresh)
print(fdepthThresh)






########
#### It is required to run the scripts above for the thresholds before proceeding
#### Script to examine feature filtering threshold: # of samples (presG) with at least # (fdepthThreshs) of reads for each taxa
########

fdepthThreshs<-c(1,50,100,500,1000,1500)
presentThreshs<-c(1,2,5)
  
  
outputFolder="outputBoston_Filt1_SampPres1" # no filter
### Kraken2 Bracken
load(here:::here(outputFolder,"brackenBoston.rdata"))

presG <- 1:ncol(MRcounts(metaobj))
fdepthThreshs<-c(1,50,100,500,1000,1500)

mat=matrix(NA, nrow=length(presG), ncol=length(fdepthThreshs))
colnames(mat)=fdepthThreshs
rownames(mat)=presG

for(i in presG){
  for(j in fdepthThreshs){
    mat[i,as.character(j)]= length(which(rowSums(MRcounts(metaobj, norm=F, log=F)>=j)>= i))
  }
}

mat2<-melt(mat)
mat2$X2=factor(mat2$X2, levels=fdepthThreshs)
p0=ggplot(mat2, aes(x=X1, y=value, group=X2)) +
  geom_line(aes(col=X2))+
  geom_point(aes(col=X2))+ggtitle("Feature filtering on Boston-SG-KB")+
  xlab("Required sample counts with feature")+ylab("# of remaining features")



########
#### Script to compare the correlation coefficients between Boston-SG-KB data and Boston-16S data at each taxa level, given the filter settings
####
olfile="Table1_OverlappingCounts.xlsx"


#fdepthThreshs<-c(1,50,100,500,1000,1500)
#presentThreshs<-c(1,2,5)
taxan=c("species","genus", "family", "order","class")
corlist<-lapply(fdepthThreshs, function(fdepthThresh){
  cor1<-lapply(presentThreshs, function(presentThresh){
    
    outputFolder=paste0("outputBoston_Filt",fdepthThresh,"_SampPres",presentThresh)
    
    corComp<-read.xlsx(here:::here(outputFolder, "comb",olfile), 4, header=TRUE)[,2]
    #print(c(fdepthThresh, presentThresh, corComp))
  })
  names(cor1)=presentThreshs
  f1<- ldply(cor1, data.frame)
  f1<- data.frame(taxa=rep(taxan, length(cor1)), f1)
  colnames(f1)=c("taxa","presentThresh","cor")
  f1
})

names(corlist)=fdepthThreshs
df1<- ldply(corlist, data.frame)
colnames(df1)[1]="fdepthThreshs"
df1$taxa=factor(df1$taxa, levels=taxan)
df1$fdepthThreshs=factor(df1$fdepthThreshs, levels=fdepthThreshs)


p1<-ggplot(df1[df1$presentThresh==1,], aes(x=taxa, y=cor, group=fdepthThreshs)) +
  geom_line(aes(col=fdepthThreshs))+
  geom_point(aes(col=fdepthThreshs))+ggtitle("In at least 1 sample")

p2<-ggplot(df1[df1$presentThresh==2,], aes(x=taxa, y=cor, group=fdepthThreshs)) +
  geom_line(aes(col=fdepthThreshs))+
  geom_point(aes(col=fdepthThreshs))+ggtitle("In at least 2 samples")

p5<-ggplot(df1[df1$presentThresh==5,], aes(x=taxa, y=cor, group=fdepthThreshs)) +
  geom_line(aes(col=fdepthThreshs))+
  geom_point(aes(col=fdepthThreshs))+ggtitle("In at least 5 samples")

pdf(here:::here(outputFolder, "FigS2_BostonData-ThresholdEval.pdf"))
ggarrange(p0,p1,p2,p5, common.legend = T, ncol = 2, nrow = 2)
dev.off()
ggarrange(p0,p1,p2,p5, common.legend = T, ncol = 2, nrow = 2)
