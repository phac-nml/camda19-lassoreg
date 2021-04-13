require(here)



#tlevel ="species" 
#outputFolder="output_Mystery"
#models <- c("glmnet_MSE", "glmnet_Classification", "glmnet_Classification_weighted", "RFClassification", "RFClassification_balanced")
models <- c(models, paste0(models,"_pairedsingle"))


########################################	
#### loading true mystery labels and gps
########################################
truth<- read.table(here:::here("data","MysterySamples_Labels_2.txt"), sep="\t", header=F, as.is=T)
colnames(truth)=c("City","LatitudeDMS","LongitudeDMS","Samples")

truth$Latitude=measurements::conv_unit(truth$LatitudeDMS, from = 'deg_min_sec', to ='dec_deg')
truth$Longitude=measurements::conv_unit(truth$LongitudeDMS, from = 'deg_min_sec', to ='dec_deg')


truthAll<-data.frame(do.call(rbind,apply(truth, 1, function(y){
  tt<-unlist(strsplit(y[4],", ")[[1]])
  cbind(tt, sample=gsub("CAMDA2019_","",tt),y[1],y[5],y[6])
})))
colnames(truthAll)[3:5]=c("City","Latitude","Longitude")
truthAll[,"Latitude"]=as.numeric(as.character(truthAll[,"Latitude"]))
truthAll[,"Longitude"]=as.numeric(as.character(truthAll[,"Longitude"]))



##### Any cities had higher SE?

seAdd2<-function(vdat){
  toret<-data.frame(vdat, longSE=(vdat$Longitude.x-vdat$Longitude.y)^2,
                    latSE=(vdat$Latitude.x-vdat$Latitude.y)^2)
  toret<-data.frame(toret,llSE=toret$longSE+toret$latSE)
  toret
}

#### function to get MSE from models 
getMSE <- function(modelDir, tlevel, truthAll){
  finalPred<-read.table(here:::here(outputFolder,paste0("Prediction_",modelDir,"_",tlevel,"_2and1.txt")), header=T, as.is=T)
  	

  mot<-merge(finalPred, truthAll, by.x="ID",by.y="sample") ## just use the pred from paired if available, otherwise, single
 
   
  if(length(grep("glmnet_MSE", modelDir))==1){  
      
     mot<-seAdd2(mot)
    mot$City2=sapply(strsplit(as.character(mot$City)," \\("), function(y)y[1])
    mcities<-read.table(here:::here("data","_Cities_Mystery.txt"), sep="\t", header=T, quote="")
    
    
    mean(mot$longSE);mean(mot$latSE);mean(mot$llSE)
    print(paste(c(modelDir,round(c(mean(mot$longSE),mean(mot$latSE),mean(mot$llSE)),2)), collapse="; "))
    ### glmnet_MSE w SAC: 8468.868, 1435.015; 9903.883
    ### glmnet_MSE; 2285.7; 1103.41; 3389.11
  }
  
  
  ####################################
  ### for classification model
  ####################################
  
  if(length(grep("Classification", modelDir))==1){
    cities<-read.table(here:::here("data","_Cities.txt"), sep="\t", header=T, quote="")
    cities$City3=paste(cities$ContinentShort, cities$Code,sep="_")
    classCol=ifelse(sum(colnames(mot)%in%"s0")==1,"s0", "x") ## s0 from glmnet, x frmo rf
    
    mot<-merge(mot, cities, by.x=classCol,by.y="Code", sort=F)
    #mot<-merge(mot, cities, by.x="s0",by.y="Code", sort=F)
    
    confusionTest <- as.data.frame(table(mot$City3, mot$City.x))
    colnames(confusionTest)[1:2]=c("Prediction","Reference")
    
    mot<-seAdd2(mot)
    mot$City2=sapply(strsplit(as.character(mot$City.x)," \\("), function(y)y[1])
    
    mean(mot$longSE);mean(mot$latSE);mean(mot$llSE)
    print(paste(c(modelDir,round(c(mean(mot$longSE),mean(mot$latSE),mean(mot$llSE)),2)), collapse="; "))
  }
  ### glmnet_Classification w SAC: 6823.992; 1826.642;8650.634
  ### glmnet_Classification; 2554.01; 1091.56; 3645.57
  ### RFClassification: 5610.93; 1634.77; 7245.71
    mot  
}
  


 
allmses <- lapply( models, function(modelDir, tlevel, truthAll){
  getMSE(modelDir, tlevel, truthAll)
}, tlevel=tlevel, truthAll=truthAll)
names(allmses)=models



#[1] "glmnet_MSE; 2629.85; 1037.63; 3667.48"
#[1] "glmnet_Classification; 2340.23; 1199.41; 3539.65"
#[1] "glmnet_Classification_weighted; 2339.42; 1143.58; 3483"
#[1] "RFClassification; 4387.24; 1885.91; 6273.15"
#[1] "RFClassification_balanced; 4174.88; 1555.62; 5730.5"
#[1] "glmnet_MSE_pairedsingle; 20238.83; 1329.48; 21568.31"
#[1] "glmnet_Classification_pairedsingle; 7847.46; 1468.31; 9315.76"
#[1] "glmnet_Classification_weighted_pairedsingle; 7847.46; 1468.31; 9315.76"
#[1] "RFClassification_pairedsingle; 15026.24; 2148.28; 17174.51"
#[1] "RFClassification_balanced_pairedsingle; 20956.34; 2028.23; 22984.57"

mets<- c("latSE", "longSE", "llSE")
meanSE <- sapply(allmses, function(x){ round(apply(x[,mets], 2, mean), 2)})
 # round(c(mean(x$longSE),mean(x$latSE),mean(x$llSE)), 2)})
colnames(meanSE)= models
rownames(meanSE)= mets
print(t(meanSE))


write.table(t(meanSE), file=here:::here(outputFolder, paste0("MSE_comparison_",tlevel,"_FullTable_byModel.txt")), col.names=T, row.names=T, quote=F,sep="\t")






  
  ##### ###   
  ### comparing GLM reg vs classification
  ########
motmse=allmses[["glmnet_MSE"]]
motclass=allmses[["glmnet_Classification"]]
  motboth <- merge(motmse, motclass, by="ID",all=T)
  comp2models<-motboth[,c("ID","City","llSE.x", "llSE.y")]
  comp2models<-comp2models[order(comp2models$City),]
  comp2models
  colnames(comp2models)=c("ID","City","llSE_glmnet_MSE","llSE_glmnet_Class")
  
  write.table(cbind(comp2models[,1:2], round(comp2models[,3:4],2)), file=here:::here(outputFolder, paste0("MSE_comparison_",tlevel,".txt")), col.names=T, row.names=F, quote=F,sep="\t")
  
  write.table(motboth, file=here:::here(outputFolder, paste0("MSE_comparison_mse2class_",tlevel,"_FullTable.txt")), col.names=T, row.names=F, quote=F,sep="\t")
  
  
  
  ### if squared errors differ between reg and class prediction, by city : one-sided Wilcoxon test with BH adj
  gstat<-t(sapply(unique(comp2models$City), function(y){
    tmp<-comp2models[comp2models[,2]==y,c(3,4)]
    if(median(tmp[,1])>median(tmp[,2])){
      pp<-wilcox.test(tmp[,1], tmp[,2], alternative="greater")$p.value
    }else{
      pp<-wilcox.test(tmp[,1], tmp[,2], alternative="less")$p.value
    }
    c(pp, median(tmp[,1]), median(tmp[,2]))
  }))
  rownames(gstat)=unique(comp2models$City)
  colnames(gstat)=c("p","medianLM","medianClass")
  gstatbh<-p.adjust(gstat[,1], method="BH")
  print(data.frame(gstat, gstatbh))
  write.table(data.frame(city=rownames(gstat),gstat, gstatbh), file=here:::here(outputFolder, paste0("MSE_comparison_mse2class_",tlevel,"_FullTable_citylevel.txt")), col.names=T, row.names=F, quote=F,sep="\t")
  
  
  save.image(here:::here(outputFolder, paste0("MSE_comparison_",tlevel,".rdata")))

