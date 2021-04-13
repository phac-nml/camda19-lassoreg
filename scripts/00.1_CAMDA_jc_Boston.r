###########
### processing Boston pilot data
###########

#require(here)

#outputFolder="outputBoston_Filt100_SampPres1"
#outputFolder="outputBoston_Filt1_SampPres1"

if(!dir.exists(here:::here(outputFolder))){dir.create(here:::here(outputFolder))}



fileNameAdd="brackenBoston"


abd<-read.table(here:::here("data/BostonPilot","biom_table_bracken.tsv"), sep="\t", header=T, comment.char = "", skip=1, quote="")
colnames(abd)=gsub("CAMDA2019_","M_",gsub("CAMDA19_MetaSUB_|_annon_profile|_profile|_bracken", "", gsub("CAMDA19_MetaSUB_pilot","P",colnames(abd))))


rownames(abd)=abd[,1]
dim(abd)
#[1]  4919   25



### removing human reads
abdHuman <- abd[which(abd[,1]==9606),]
abd<-abd[-which(abd[,1]==9606),]

OTUannot<-abd[,c(1,ncol(abd))]
OTUannot2<-do.call(rbind,strsplit(as.character(OTUannot$taxonomy), ";"))
colnames(OTUannot2)=c("kingdom","phylum","class","order","family","genus","species")
OTUannot2=gsub("k__|p__|c__|o__|f__|g__|s__| ","",OTUannot2)
OTUannot<-data.frame(OTUannot,OTUannot2)


abd<-data.matrix(abd[,-c(1,ncol(abd))])

library(xlsx)
### special data wrangling 
sampAnnot<-read.xlsx(here:::here("data/BostonPilot","BostonPilot_S1_subset_Julie.xlsx"), 1, header=TRUE)
abd<-abd[ ,as.character(sampAnnot$NA.)]
colnames(abd)=sampAnnot$sampleID
rownames(sampAnnot)=sampAnnot$sampleID


dat<-abd
print(dim(dat))
# 4918   23


###############
#### setup the MRexperiment object
###############
require(metagenomeSeq)
phenotypeData = AnnotatedDataFrame(sampAnnot)
phenotypeData


OTUdata = AnnotatedDataFrame(OTUannot)
OTUdata

### the input data is in the format of  OTU x samplID 
metaobj = newMRexperiment(dat,phenoData=phenotypeData,featureData=OTUdata)

quantile(as.vector(MRcounts(metaobj, norm=T, log=F)), probs=seq(0, 1, 0.05 ))


###############
##### setting parameters before running source_fromMRexperiment.r
###############
#fdepthThresh=100
#presentThresh=1
depthThresh=1 #depreciated
perct=""


keyFeat="sampletype"
feat2="surfacematerial"
toexport=FALSE

source(here:::here("scripts","CAMDA_jc_source_fromMRexperiment.r"))


rm(list=setdiff(ls(), "outputFolder"))



####################################################
####
#### 16S
####
##################################################
require(here)


fileNameAdd="brackenBoston16S"
library(xlsx)


##S2_16SOTU_filtered0.1in1
abd<-read.xlsx(here:::here("data/BostonPilot","BostonPilot_S2_subset_Julie.xlsx"), 1, header=TRUE)
rownames(abd)=abd[,1]
dim(abd)


OTUannot<-abd[,c(1,ncol(abd))]
OTUannot2<-do.call(rbind,strsplit(as.character(OTUannot$taxonomy), ";"))
colnames(OTUannot2)=c("kingdom","phylum","class","order","family","genus","species")
OTUannot2=gsub("k__|p__|c__|o__|f__|g__|s__| ","",OTUannot2)
OTUannot<-data.frame(OTUannot,OTUannot2)

abd<-data.matrix(abd[,-c(1,ncol(abd))])

sampAnnot<-read.xlsx(here:::here("data/BostonPilot","BostonPilot_S1_subset_Julie.xlsx"), 1, header=TRUE)
rownames(sampAnnot)=sampAnnot$sampleID

## should be 0
sum(!colnames(abd)%in%rownames(sampAnnot))
###make sure the sample orders are the same
abd<-abd[ ,rownames(sampAnnot)]


dat<-abd
print(dim(dat))
# 2134   23

############################
#### setup the MRexperiment object
############################
require(metagenomeSeq)
phenotypeData = AnnotatedDataFrame(sampAnnot)
phenotypeData


OTUdata = AnnotatedDataFrame(OTUannot)
OTUdata

### the input data is in the format of  OTU x samplID 
metaobj = newMRexperiment(dat ,phenoData=phenotypeData ,featureData=OTUdata)


fdepthThresh=1
presentThresh=1
depthThresh=1
perct=""

keyFeat="sampletype"
feat2="surfacematerial"
toexport=FALSE

source(here:::here("scripts","CAMDA_jc_source_fromMRexperiment.r"))


#rm(list=setdiff(ls(), "outputFolder"))

###
####################################################
####
#### metaphlan2
####
##################################################
require(here)

fileNameAdd="brackenBostonMetaphlan"
library(xlsx)

##S2_16SOTU_filtered0.1in1
abd<-read.xlsx(here:::here("data/BostonPilot","BostonPilot_S2_subset_Julie.xlsx"), 2, header=TRUE)
rownames(abd)=abd[,1]
dim(abd)
#1340   24

OTUannot<-abd[,1]

OTUannot2<-t(sapply(strsplit(as.character(OTUannot), "\\|"), function(y){
  c(y,rep("",8-length(y)))
}))
colnames(OTUannot2)=c("kingdom","phylum","class","order","family","genus","species","t")
OTUannot2=gsub("k__|p__|c__|o__|f__|g__|s__|t__| ","",OTUannot2)
OTUannot<-data.frame(OTUannot,OTUannot2)
rownames(OTUannot)=OTUannot$OTUannot

OTUannot$species<- sapply(strsplit(OTUannot$species,"_"), function(y){if(length(y)==0){toret=""}else{toret=y[2]}})
##unclassified 115; "" 457


abd<-data.matrix(abd[,-1])

sampAnnot<-read.xlsx(here:::here("data/BostonPilot","BostonPilot_S1_subset_Julie.xlsx"), 1, header=TRUE)
rownames(sampAnnot)=sampAnnot$sampleID

## should be 0
sum(!colnames(abd)%in%rownames(sampAnnot))
###make sure the sample orders are the same
abd<-abd[ ,rownames(sampAnnot)]


dat<-abd
print(dim(dat))
#  1340   23

############################
#### setup the MRexperiment object
############################
require(metagenomeSeq)
phenotypeData = AnnotatedDataFrame(sampAnnot)
phenotypeData


OTUdata = AnnotatedDataFrame(OTUannot)
OTUdata

### the input data is in the format of  OTU x samplID 
metaobj = newMRexperiment(dat ,phenoData=phenotypeData ,featureData=OTUdata)


fdepthThresh=1
presentThresh=1
depthThresh=1
perct=""

keyFeat="sampletype"
feat2="surfacematerial"
toexport=FALSE

source(here:::here("scripts","CAMDA_jc_source_fromMRexperiment.r"))


#rm(list = ls(all.names = TRUE))
#rm(list=setdiff(ls(), "outputFolder"))

print(outputFolder)
