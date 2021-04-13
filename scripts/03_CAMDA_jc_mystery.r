
require(here)

#here()

#tlevel="species"
#tlevel="genus"
#tlevel="family"


#fileNameAdd=paste("Mystery_New_",tlevel); mfold="data/Mystery_paired"; 
#fileNameAdd=paste("Mystery_New_single",tlevel); mfold="data/Mystery_single"


for(fileNameAdd in c(paste0("Mystery_New_",tlevel), paste0("Mystery_New_single_",tlevel))){
  if(fileNameAdd==paste0("Mystery_New_",tlevel)){
    mfold="data/Mystery_paired"
    outputFolder="output_mystery"
    
    
  }else{
    mfold="data/Mystery_single"
    outputFolder="output_mystery_single"
    
  }
  if(!dir.exists(here:::here(outputFolder))){dir.create(here:::here(outputFolder))}
  
  abd<-read.table(here:::here(mfold, paste0("bracken_",tlevel,".txt")), sep="\t", header=T, comment.char = "", skip=1, quote="")
  rownames(abd)=abd[,1]
  
  
  
  ### removing human reads
  if(length(which(abd[,1]==9606))>0){
  	abdHuman <- abd[which(abd[,1]==9606),]
  	abd<-abd[-which(abd[,1]==9606),]
  }
  
  
  OTUannot<-abd[,c(1,ncol(abd))]
  OTUannot2<-strsplit(as.character(OTUannot$taxonomy), ";")
  OTUannot2<-do.call(rbind,lapply(OTUannot2, function(y){c(y,rep("",7-length(y)))}))
  colnames(OTUannot2)=c("kingdom","phylum","class","order","family","genus","species")
  OTUannot2<-gsub("k__|p__|c__|o__|f__|g__|s__| ","",OTUannot2)
  OTUannot<-data.frame(OTUannot,OTUannot2)
  
  
  abd<-data.matrix(abd[,-c(1,ncol(abd))])
  colnames(abd)=gsub("CAMDA2019_","M_",gsub("CAMDA19_MetaSUB_|_annon_profile|_profile|_bracken|_annon.fastq", "", gsub("CAMDA19_MetaSUB_pilot","P",colnames(abd))))
  
  
  dat=t(abd)
  print(dim(dat))
  #40 4796 pair
  #100 5449 single
  
  ###
  sampAnnot<-as.matrix(rownames(dat))
  sampAnnot<- cbind(sampAnnot, t(sapply(strsplit(rownames(dat), "_"), function(y)y[1:3])))
  sampAnnot[is.na(sampAnnot[,4]),4]="single"
  sampAnnot<- data.frame(sampAnnot)
  
  colnames(sampAnnot)=c("ID","Batch","sample","pair")
  rownames(sampAnnot)=sampAnnot[,"ID"]
  
  ###############
  #### setup the MRexperiment object
  ###############
  require(metagenomeSeq)
  phenotypeData = AnnotatedDataFrame(sampAnnot)
  phenotypeData
  
  
  OTUdata = AnnotatedDataFrame(OTUannot)
  OTUdata
  
  metaobj = newMRexperiment(t(dat),phenoData=phenotypeData,featureData=OTUdata)
  
  
  ### no filters for mystery samples, as we want to make a prediction on all samplesm
  fdepthThresh=1
  presentThresh=1
  depthThresh=1
  perct=""
  
  keyFeat="sample"
  feat2="pair"
  toexport=FALSE
  
  source(here:::here("scripts","CAMDA_jc_source_fromMRexperiment.r"))

}