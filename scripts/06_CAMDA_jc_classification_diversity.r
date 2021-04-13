
library(here)

#outputFolder="output"
#gpsw = "Classification"
#dir2save = paste0("glmnet_",gpsw)
#tlevel="species" ; #tlevel="genus"

load(here:::here(outputFolder,paste0(dir2save,"/Diversity_",tlevel,".rdata")))

wilcox.test(datSimp[datSimp=="L1CO",2], datSimp[datSimp=="CV10",2])$p.value

#############
#### building naive bayes classifier using a nonparametric approach(gaussian kernel) to classify certain (seen) from uncertain(unseen) predictions
#############
#library(naivebayes)
nb_kde <- naive_bayes(group ~ ., datSimp, usekernel = TRUE)
summary(nb_kde)
get_cond_dist(nb_kde)

# Visualize class conditional densities
pdf(here:::here(outputFolder,paste0(dir2save,"/Fig6_Diversity_Fitted_Probabiliy_",gpsw,"_",tlevel,"_combined_class_conditional_densities.pdf")), width=3.5,height=3.5)
plot(nb_kde, "value", arg.num = list(legend.cex = 0.75, xlab="Simpson diversity index", col=c("#F8766D","#00BFC4")), prob = "conditional")
dev.off()

### load the prediction probabiliies of mystery data for uncertainty prediction
load(here:::here("output_mystery","MysterySamples_Labels_3.rdata"))
load(here:::here("output_mystery",paste0("Prediction_",dir2save,"_",tlevel,"TRUETRUE.rdata")))
mysSimpsonFull<-mysSimpsonp<-mysSimpson
load(here:::here("output_mystery_single",paste0("Prediction_",dir2save,"_",tlevel,"TRUETRUE.rdata")))
mysSimpsonFull<-c(mysSimpsonFull,mysSimpson)
mysSimpsonFull<-mysSimpsonFull[-grep("_1$|_2$",names(mysSimpsonFull))]



# Obtain Posterior probabilities
seentf<-nb_kde %prob% data.frame(value=mysSimpsonFull)
rownames(seentf)=names(mysSimpsonFull)
seentf<-data.frame(sample=gsub("M_","",names(mysSimpsonFull)), seentf, seentf=ifelse(seentf[,"CV10"] < seentf[,"L1CO"], "unseen","seen"))

table(seentf$seentf)


### load data annotation
load(here:::here("output_mystery","MysterySamples_Labels_3.rdata"))
seenplus<-merge(seentf, truthAll,by="sample") ## just use the pred from paired if available, otherwise, single
table(seenplus$City,seenplus$seentf)

### sensitivity on the mystery samples (all from new origins)
colSums(table(seenplus$City,seenplus$seentf))[2]/nrow(seenplus)


pdf(here:::here(outputFolder,paste0(dir2save,"/Diversity_NaiveBayes_seenTF_ROC_",gpsw,"_",tlevel,"_mysteryPred.pdf")), width=6,height=5.5) #,"_",reg,"log_",log0
  par(mai=c(1,2.1,2,1))
  newoldCnts<-t(table(seenplus$City,seenplus$seentf))
  barplot(newoldCnts[2:1,], horiz=T, las=1, col=c("#00BFC4","#F8766D"))
  legend("topright",legend=c("New origin","Pre-trained"),fill=c("#00BFC4","#F8766D"), cex=0.8)
dev.off()
#save(seentf,seenplus,mysSimpsonFull,nb_kde,datSimp, file=here:::here(outputFolder,paste0(dir2save,"/diversityPred_Mystery_",tlevel,".rdata")))


#### 
### NB leave1out for TPR and FPR plot
#### 
nbpred<-do.call(rbind,lapply(1:nrow(datSimp), function(i){
	nb_kde0 <- naive_bayes(group ~ ., datSimp[-i,], usekernel = TRUE)
	# Obtain Posterior probabilities
	seentf0<-nb_kde0 %prob% data.frame(value=datSimp[i,"value"])

	seentf0
}))

#library(ROCR)
pred<-prediction(nbpred[,2], datSimp[,1])
perf<-performance(pred,"tpr","fpr")
pdf(here:::here(outputFolder,paste0(dir2save,"/Fig6_Diversity_NaiveBayes_seenTF_ROC_",gpsw,"_",tlevel,".pdf")), width=3.6,height=4) #,"_",reg,"log_",log0
	plot(perf, main="Naive Bayes Diversity Classification")
	abline(a=0, b= 1, col="grey")
	text(0.7,0.2, paste0("AUC = ",signif(performance(pred,"auc")@y.values[[1]],2)))

dev.off()
performance(pred,"auc")@y.values[[1]]
 #AUC = 0.8595096
 
xtab<-table(levels(datSimp[,1])[apply(nbpred, 1, which.max)], datSimp[,1])
#     CV10 L1CO
#CV10  247   59
#L1CO   29  217
library(caret) 
(confus<-confusionMatrix(xtab, positive ="L1CO"))
#Accuracy : 0.8406
#Sensitivity :  0.7862          
#Specificity : 0.8949 

save.image(here:::here(outputFolder,paste0(dir2save,"/image_diversityPred_",tlevel,".rdata")))
