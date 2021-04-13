#source activate mro_env
#srun -c 4 --mem 128G --partition NMLResearch R --vanilla

library(metagenomeSeq)
library(here)


here()


#tlevel="species"
#tlevel="genus"
#tlevel="family"

#fileNameAdd=paste0("bracken_",tlevel)
#outputFolder="output_pairedsingle"
dir.create(outputFolder, showWarnings = FALSE)


abd0<-read.table(here("data/MetaSUB-WGS_paired",paste0("biom_table_bracken_",tlevel,".tsv")), sep="\t", header=T, comment.char = "", skip=1, quote="")
rownames(abd0)=abd0[,1]

abd1<-read.table(here("data/MetaSUB-WGS_single",paste0("biom_table_bracken_",tlevel,".tsv")), sep="\t", header=T, comment.char = "", skip=1, quote="")
rownames(abd1)=abd1[,1]
### removing the singles from paired
abd1<-abd1[,-grep("_2.fastq_bracken|_1.fastq_bracken",colnames(abd1))]


abd<-merge(abd0, abd1, by="taxonomy", all=T, sort=F)
abd[is.na(abd)]=0
abd$X.OTU.ID=abd$X.OTU.ID.x
abd<-abd[,!colnames(abd)%in%c("X.OTU.ID.x","X.OTU.ID.y")]



### removing human reads
if(length(which(abd[,"X.OTU.ID"]==9606))>0){
  abdHuman <- abd[which(abd[,"X.OTU.ID"]==9606),]
  abd<-abd[-which(abd[,"X.OTU.ID"]==9606),]
}

OTUannot<-abd[,c(1,ncol(abd))]
OTUannot2<-strsplit(as.character(OTUannot$taxonomy), ";")
OTUannot2<-do.call(rbind,lapply(OTUannot2, function(y){c(y,rep("",7-length(y)))}))
colnames(OTUannot2)=c("kingdom","phylum","class","order","family","genus","species")
OTUannot2<-gsub("k__|p__|c__|o__|f__|g__|s__| ","",OTUannot2)
OTUannot<-data.frame(OTUannot,OTUannot2)


abd<-data.matrix(abd[,-c(1,ncol(abd))])
colnames(abd)=gsub("CAMDA2019_","M_",gsub("CAMDA19_MetaSUB_|_annon_profile|_profile|_bracken|.fastq", "", gsub("CAMDA19_MetaSUB_pilot","P",colnames(abd))))


dat=t(abd)
print(dim(dat))
#294 6922 species
#294 1736 genus

###
sampAnnot<-rownames(dat)
sampAnnot<-cbind(ID=sampAnnot, t(sapply(strsplit(rownames(dat), "_"), function(y)y[1:2])))
rownames(sampAnnot)<-rownames(dat)
sampAnnot<-data.frame(sampAnnot, t(apply(sampAnnot, 1,function(y){if(length(grep("\\.", y[3]))!=0){mp<-strsplit(y[3],"\\.")}else{mp=c(y[3],NA)}; unlist(mp) })))
colnames(sampAnnot)[2:5]= c("Batch","ID2","City","ID3")

mysSampID<-grep("^M_",rownames(dat))
sampAnnot$City2 <-as.character(sampAnnot$City)
sampAnnot$City2[mysSampID]="_Mys"
sampAnnot$City2=factor(sampAnnot$City2)

sampAnnot$contin<-factor(ifelse(sampAnnot$City2%in%c("AKL","HAM"), "Oceania",
                                ifelse(sampAnnot$City2%in%c("ILR","OFA"), "Africa",
                                       ifelse(sampAnnot$City2%in%c("BER","LON", "MAR", "PXO", "SOF","STO"),"Europe",
                                              ifelse(sampAnnot$City2%in%c("NYC","SAC"), "NorthA",
                                                     ifelse(sampAnnot$City2%in%c("SAO", "BOG"), "SouthA",
                                                            ifelse(sampAnnot$City2%in%c("HGK", "TOK"), "Asia",
                                                                   "Mystery")))))))
continConv<-function(a){
  ifelse(a%in%c("AKL","HAM"), "Oceania",
         ifelse(a%in%c("ILR","OFA"), "Africa",
                ifelse(a%in%c("BER","LON", "MAR", "PXO", "SOF","STO"),"Europe",
                       ifelse(a%in%c("NYC","SAC"), "NorthA",
                              ifelse(a%in%c("SAO", "BOG"), "SouthA",
                                     ifelse(a%in%c("HGK", "TOK"), "Asia",
                                            "Mystery"))))))
}
#write.table(sampAnnot,"sampAnnot_jc.txt",row.names=F,col.names=T, quote=F,sep="\t")

sampAnnotOrig=sampAnnot


###############
#### setup the MRexperiment object
###############
print("Setting up the MRexperiment object")
phenotypeData = AnnotatedDataFrame(sampAnnot)
phenotypeData


OTUdata = AnnotatedDataFrame(OTUannot)
OTUdata

metaobj = newMRexperiment(t(dat),phenoData=phenotypeData,featureData=OTUdata)



#presentThresh=8
#depthThresh=1
perct=""

keyFeat="City2"
feat2="contin"
toexport=FALSE

clplus=factor(continConv(levels(pData(metaobj)[,keyFeat])))
source(here("scripts","CAMDA_jc_source_fromMRexperiment.r"))




