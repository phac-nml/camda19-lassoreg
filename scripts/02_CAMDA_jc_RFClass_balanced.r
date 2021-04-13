
require(metagenomeSeq)
require(randomForest)
require(caret)
require(plyr)
require(here)
require(parallel)

#### CAMDA data
##### Classification

#outputFolder="output"
#outputFolder="output_pairedsingle"

load(here:::here(outputFolder,paste0("bracken_",tlevel,".rdata")))
#load(here:::here(outputFolder,"bracken_species.rdata"))
#load(here:::here(outputFolder,"bracken_genus.rdata"))
#load(here:::here(outputFolder,"bracken_family.rdata"))


######################################
### Classification
######################################
cities<-read.table(here:::here("data","_Cities.txt"), sep="\t", header=T, quote="")
tmpc<-merge(data.frame(pData(metaobj)), cities, by.x="City2",by.y="Code",all.x=T, sort=F)
tmpc2<- cbind(City2=as.character(tmpc$City2), City3=paste0(tmpc$ContinentShort,"_", tmpc$City2))
rownames(tmpc2)=tmpc2[,"City2"]


outc <- tmpc$City2

dir2save = "RFClassification_balanced"

if(!dir.exists(here:::here(outputFolder,dir2save))){dir.create(here:::here(outputFolder,dir2save))}
log0=T
norm0=T
#norm0=F

print(tlevel)
dat2model<-MRcounts(metaobj, norm=norm0, log=log0)
rownames(dat2model)=fData(metaobj)[,tlevel]
dat2model<-dat2model[rownames(dat2model)!="",]  
rownames(dat2model)[grep("^[0-9]",rownames(dat2model))]=paste0("zzz-",rownames(dat2model)[grep("^[0-9]",rownames(dat2model))])





########################	 
### model building with the list of training and validation datasets
######################## 
dat4train=data.frame(cate=outc, t(dat2model), stringsAsFactors=F)
seedNums=1:10

  ### for imbalanced data; 2 fold or more counts, sample balanced
  sampCnt<-as.vector(table(dat4train[,1]))
  sampCntMinV<-rep(min(sampCnt), length(sampCnt))

ntree=500 # mtry 64 default
rfResults<-mclapply(seedNums, function(seedNum){
  set.seed(seedNum)
  #Fit Random Forest Model
  rf1 = randomForest(cate ~ .,  
                     ntree = ntree,
                     sampsize = sampCntMinV,
                     data = dat4train)
  print(rf1) 
  
  # Variable Importance
  varImpPlot(rf1,   sort = T, n.var=10, main="Top 10 - Variable Importance")
  var.imp = data.frame(importance(rf1,type=2))
  var.imp$Variables = row.names(var.imp)  
  #print(var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),])
  
  ### this is from OOB result
  print(cmv<-confusionMatrix(data= rf1$predicted,  
                             reference=dat4train$cate))
  
  
  list(model=rf1, var.imp=var.imp, confusionm=cmv)
}, mc.cores=2)
print("randomForestDone")

  

save(ntree, rfResults,dat4train, outc, file=here:::here(outputFolder,paste0(dir2save,"/RFclassification_",tlevel,"log_",log0,".rdata")))




###credit to https://stackoverflow.com/questions/37897252/plot-confusion-matrix-in-r-using-ggplot

pdf(here:::here(outputFolder,paste0(dir2save,"/FigS4_RFclassification_",tlevel,"log_",log0,"_pred_outofbag_contsorted.pdf")), width=7.5, height=6)
for(i in seedNums){
  ### plot in continent cluster
  confusion_matrix <- as.data.frame(rfResults[[i]]$confusionm$table)
  confusion_matrix[,1]=tmpc2[as.character(confusion_matrix[,1]),"City3"]
  confusion_matrix[,2]=tmpc2[as.character(confusion_matrix[,2]),"City3"]
  #oacc <- sum(diag(rfResults[[i]]$model$confusion))/sum(rfResults[[i]]$model$confusion)
  
  plotg<-ggplot(data = confusion_matrix,
                mapping = aes(x = Reference,
                              y = Prediction)) +
    geom_tile(aes(fill = Freq)) +
    geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
    scale_fill_gradient(low = "lavenderblush2",
                        high = "coral2", #lavenderblush3
                        na.value = "lavenderblush",
                        trans = "log")+
    ggtitle(paste("OOB seed",i, ": Overall accuracy =", signif( rfResults[[i]]$confusionm$overall["Accuracy"],4)))+
    theme(axis.text.x = element_text(angle =45, hjust = 1))
  print(plotg)
}
dev.off()
