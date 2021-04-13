
###################
### prediction from GLMNET model


##########################################
##### making predictions: currently on training data, use *_gps_scatter_CV_pred_lambda.1se.pdf instead
### Default is the value s="lambda.1se"
##########################################

library(glmnet)
require(metagenomeSeq)
library(here)

#tlevel="species"; #tlevel="genus"
#dataDir="output"; #dataDir="output_pairedsingle"
#modelDir="glmnet_MSE";model2use=paste0(modelDir,"/Multivariate_",tlevel,"_lassolog_TRUE_gps_MSE.rdata")

if(dataDir=="output_pairedsingle"){
  modelDir=paste0(modelDir,"_pairedsingle")
}
modelDir


for (mdataDir in c("output_mystery","output_mystery_single")){
  
  print(mdataDir)
  
  #load mystery data
  if(mdataDir=="output_mystery"){
    load(here:::here(mdataDir, paste0("Mystery_New_",tlevel,".rdata")))
  }else{
    load(here:::here(mdataDir, paste0("Mystery_New_single_",tlevel,".rdata")))
  }
  
  #load model
  print(here:::here(dataDir,model2use))
  load(here:::here(dataDir,model2use))
  
  
  log0=T
  norm0=T
  #norm0=F
  
  print(tlevel)
  dat2test<-MRcounts(metaobj, norm=norm0, log=log0)
  rownames(dat2test)=fData(metaobj)[,tlevel]
  #dat2test<-cbind(MRcounts(wgsList[[tlevel]], norm=F, log=F), MRcounts(wgsList$family, norm=F, log=F))
  dat2test<-dat2test[rownames(dat2test)!="",]  
  rownames(dat2test)[grep("^[0-9]",rownames(dat2test))]=paste0("zzz-",rownames(dat2test)[grep("^[0-9]",rownames(dat2test))])
  
  loadModelledData<-function (tlevel,norm0,log0 ){
  	load(here:::here(dataDir,paste0("bracken_",tlevel,".rdata")))
  	
  	dat2model<-MRcounts(metaobj, norm=norm0, log=log0)
  	rownames(dat2model)=fData(metaobj)[,tlevel]
  	#dat2model<-cbind(MRcounts(wgsList[[tlevel]], norm=F, log=F), MRcounts(wgsList$family, norm=F, log=F))
  	dat2model<-dat2model[rownames(dat2model)!="",]  
  	rownames(dat2model)[grep("^[0-9]",rownames(dat2model))]=paste0("zzz-",rownames(dat2model)[grep("^[0-9]",rownames(dat2model))])
  	
  	dat2model
  }
  
  modelDat <- loadModelledData(tlevel,norm0,log0)
  
  ### fill in 0 for missing taxa
  miss<-rownames(modelDat)[!rownames(modelDat)%in%rownames(dat2test)]
  missMat<-matrix(0, nrow=length(miss), ncol=ncol(dat2test)) ; rownames(missMat)=miss
  
  testDat <- rbind(dat2test, missMat)
  testDat <- testDat[rownames(modelDat),]
  
  
  

  
    ### lambda.1se on full data  #dim 3043    2    1
  if(length(grep("glmnet_Classification",modelDir))>0){
  	predGpsrn<-predict(fit1m3,newx=t(testDat))[,,1]
  	predGps<-predict(fit1m3,newx=t(testDat), type="class")
  	rownames(predGps)=rownames(predGpsrn)
  	
  	### diversity
  	predGpsR<-predict(fit1m3,newx=t(testDat), type="response")
  	require(vegan)
  	mysSimpson<-apply(predGpsR, 1, function(y){diversity(na.omit(y), "simpson")})
  	save(predGps,predGpsR, mysSimpson, file=here:::here(mdataDir, paste0("Prediction_",modelDir,"_",tlevel,norm0,log0,".rdata")))
  
  }else if(length(grep("RFClassification",modelDir))>0){
    ### symbols and duplicated species names (from diff genus) got modified
    tmp<-cbind(colnames(dat4train)[-1], rownames(testDat))
    head(tmp[which(colnames(dat4train)[-1]!=gsub("-|\\(|\\)|:",".",rownames(testDat))),])
    rownames(testDat)=colnames(dat4train)[-1]
    require(randomForest)
    
    # Predicting response variable, using the model (from 10seeds) with highest OOB Accuracy
    seedrf=which.max(sapply(1:10, function(i){ rfResults[[i]]$confusionm$overall["Accuracy"]}))
    predGps <- predict(rfResults[[seedrf]]$model ,t(testDat), type="class")
    predGpsR<-predict(rfResults[[seedrf]]$model ,t(testDat), type = "prob")
    require(vegan)
    mysSimpson<-apply(predGpsR, 1, function(y){diversity(na.omit(y), "simpson")})
    save(predGps,predGpsR, mysSimpson, file=here:::here(mdataDir, paste0("Prediction_",modelDir,"_",tlevel,norm0,log0,".rdata")))
    
  }else if(length(grep("Univariate",model2use)!=0)!=0){
  	
  	predGps<-predict(fit1m3,newx=t(testDat))
  	
  	save(predGps, file=here:::here(mdataDir, paste0("Prediction_",modelDir,"_",tlevel,norm0,log0,".rdata")))
  
  }else{
  
  	predGps<-predict(fit1m3,newx=t(testDat))[,,1]
  	
  	save(predGps, file=here:::here(mdataDir, paste0("Prediction_",modelDir,"_",tlevel,norm0,log0,".rdata")))
  
  }
  print(mdataDir)
  
  write.table(predGps, file=here:::here(mdataDir, paste0("Prediction_",modelDir,"_",tlevel,norm0,log0,".txt")),	col.names=T, row.names=T, quote=F,sep="\t")
  	
  	

}




