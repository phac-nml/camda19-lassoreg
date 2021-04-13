#seed2use=111

require(metagenomeSeq)
require(glmnet)
require(caret)
require(plyr)
require(naivebayes)
require(parallel)
#require(here)

source(here:::here("scripts","CAMDA_jc_source_glmnet.r"))


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

outc <- tmpc$City2

gpsw = "Classification_weighted"
dir2save = paste0("glmnet_",gpsw)

if(!dir.exists(here:::here(outputFolder,dir2save))){dir.create(here:::here(outputFolder,dir2save))}
log0=T
norm0=T
#norm0=F

print(tlevel)
dat2model<-MRcounts(metaobj, norm=norm0, log=log0)
rownames(dat2model)=fData(metaobj)[,tlevel]
dat2model<-dat2model[rownames(dat2model)!="",]  
rownames(dat2model)[grep("^[0-9]",rownames(dat2model))]=paste0("zzz-",rownames(dat2model)[grep("^[0-9]",rownames(dat2model))])



##########################################
###### GLMNET classification
##########################################

alpha=1;reg="lasso"
#alpha=0.75;reg="elastic0.75"
#alpha=0.5;reg="elastic0.5"
#alpha=0;reg="ridge"

### adding weights
## weight is computed to be 1-(# of samples in the category / total # of samples)
jGetWeights<-function(cate){
  cateW <- 1-(as.matrix(table(cate))/length(cate))
  cateW[as.character(cate),1]
}
w2balance<-jGetWeights(outc)


set.seed(seed2use) 
system.time(cvob2mse<-cv.glmnet(t(dat2model) , y=outc, weights=w2balance, family="multinomial", type.multinomial = "grouped", nfolds =10, alpha=alpha, parallel=T, keep=T) ) 

##########################################
###### plotting the predicted values vs true in cross validation
##########################################
dim(cvob2mse$fit.preval) # 294  16 100
cvPred1<-cvob2mse$fit.preval[,1:length(levels(outc)),which(cvob2mse$lambda==cvob2mse$lambda.1se)]
cvclass<-levels(outc)[apply(cvPred1,1, which.max)]
cvclass<-factor(cvclass, level=levels(outc))
cvconfMat<-table(cvclass,outc)

print(cmv<-confusionMatrix(data= cvclass,  
				reference=outc))
							   
### for multi-class comparisons		
print("Overall accuracies")
overallacc<-cmv$overall["Accuracy"]
print(overallacc)
print(cmv$byClass)
apply(cmv$byClass, 2, median)


require(gplots)
#heatmap.2(table(cvconfMat), scale="none",trace="none", symm=T, Rowv=F, Colv=F, dendrogram="none")#, cellnote=data.matrix(unclass(cvconfMat)))
heatmap.2(cvconfMat, scale="none",trace="none", symm=T, Rowv=F, Colv=F, dendrogram="none")#, cellnote=data.matrix(unclass(cvconfMat)))


#pdf(paste0(dir2save,"/classification_",tlevel,"_",reg,"log_",log0,"_scatter_CV_pred_lambda.1se_MSE_Loss2_",gpsw,".pdf"), height=14)
#	a=gps[,gpsw]; b=cvPred1
#	plot(a,b, xlab=paste0("True ",gpsw),ylab="CV Pred",xlim=range(c(a,b)), ylim=range(c(a,b)), main=paste(gpsw, tlevel,reg, "r2=",signif(cor(a,b)^2,4)))
#	abline(a=0,b=1) 
	
#dev.off()
#cvob2mse$glmnet.fit


fit1m=glmnet(t(dat2model) ,outc, weights=w2balance, family="multinomial", type.multinomial = "grouped", alpha=alpha) 
fit1m2=glmnet(t(dat2model) ,outc, weights=w2balance, family="multinomial", type.multinomial = "grouped", lambda=cvob2mse$lambda.min, alpha=alpha) 
fit1m3=glmnet(t(dat2model) ,outc, weights=w2balance, family="multinomial", type.multinomial = "grouped", lambda=cvob2mse$lambda.1se, alpha=alpha) 

save(cvob2mse, fit1m, fit1m3, tlevel, reg,log0,cmv, cvclass,outc, file=here:::here(outputFolder,paste0(dir2save,"/classification_",tlevel,"_",reg,"log_",log0,"_gps",gpsw,".rdata")))


which1se<-which(round(cvob2mse$glmnet.fit$lambda,6)==round(cvob2mse$lambda.1se,6))
df1se<-cvob2mse$glmnet.fit$df[which1se]
cvm1se<-cvob2mse$cvm[which1se]
print(paste("lambda 1se: ",cvob2mse$lambda.1se, cvob2mse$lambda[which1se]))
print(paste("$cvm MSE value",cvm1se))
print(paste("df",df1se))

				   



## consider using cvob2mse$lambda.1se instead, using too many variables even with lasso

#cvob2mse.f1<-cvob2mse;cvob2mse.f1<-cvob2mse
#!coef(cvob2mse.f1)<-coef(cvob2mse)[[1]]
#!coef(cvob2mse.f2)<-coef(cvob2mse)[[2]]

##########################################
###### plotting the model info and lambda estimation
##########################################
pdf(here:::here(outputFolder,paste0(dir2save,"/classification_",tlevel,"_",reg,"log_",log0,".pdf")))
	plot(fit1m,type.coef="2norm")
	fit1mRange=range(log(fit1m$lambda))
	plot(fit1m, xvar = "lambda", label = TRUE, type.coef = "2norm", xlim=c(fit1mRange[1]-0.7, fit1mRange[2]))
	abline(v=log(cvob2mse$lambda.min), lty=2)
	abline(v=log(cvob2mse$lambda.1se), lty=2)
	#plot(fit1m2, xvar = "lambda", label = TRUE, type.coef = "2norm", main="lambda.min")
	#plot(fit1m3, xvar = "lambda", label = TRUE, type.coef = "2norm", main="lambda.1se")
	
	par(mar=c(4.5,4.5,4,1)) 
	plot(cvob2mse)
	title(paste("classification", tlevel, reg),line=2.5)
	#plot(cvob2mse, xvar='lambda', label=TRUE)
dev.off()

##########################################
##### plot top 30 features
##########################################
#coef1<-coef(cvob2mse)[,1][coef(cvob2mse)[,1]!=0][-1] ## skipping the intercept
#coef1<-coef1[order(abs(coef1),decreasing=T)][1:30]
#coef1<-sort(coef1)


coefList<-lapply(1:length(coef(cvob2mse)), function(i){
	coef1<-coef(cvob2mse)[[i]][,1][coef(cvob2mse)[[i]][,1]!=0][-1] ## skipping the intercept
	coef1<-coef1[order(abs(coef1),decreasing=T)][1:30]
	coef1<-sort(coef1)
	coef1
})
names(coefList)=levels(outc)

pdf(here:::here(outputFolder,paste0(dir2save,"/classification_",tlevel,"_",reg,"log_",log0,"_featImp",gpsw,".pdf")), height=10)
	par(mar=c(4.5,10,4,1)) 
	for(i in 1:length(coefList)){
		barplot(coefList[[i]], horiz=T,las=2, main=paste(names(coefList)[i],"top30"), cex.names=0.9)}
dev.off()



##########################################
### examine abundance of top taxa along lat and long,color
##########################################
 pdf(here:::here(outputFolder,paste0(dir2save,"/classification_",tlevel,"_",reg,"log_",log0,"_featImp_abundance_",gpsw,".pdf")), width=12,height=7)
 par(mfrow=c(2,1), mar=c(4,4,2,2))
for(taxa in unique(unlist(lapply(coefList, function(y)names(y))))){

	### some share the same species names but have different genera, need to edit the script to keep Genera
	for (taxai in which(rownames(dat2model)==taxa)){
		
		tmp<-as.vector(dat2model[taxai,])
		
		dat2gg<-data.frame(outc, tmp)
		colnames(dat2gg)=c("City","abundance")
	
		plot(dat2gg, main=paste(tlevel,taxa), ylab="log2 (Abundance+1)", xlab="City",outline=F, ylim=range(dat2gg[,2]))
		points(jitter(as.numeric(dat2gg[,1])),dat2gg[,2], pch=20, col="blue")
		#library(ggplot2);
		#gg<-ggplot(dat2gg, aes(City, abundance)) + 
		#	geom_point() + 
		#	labs(y = "log2 (Abundance+1)", x = "City");
		#print(gg)
	}
}
dev.off()

		
		
		
####		
####		
#### use jLeaveOutclassification function onwards
####		
####


###############################
#### leave one city out, train the rest and predict samples from the left out city
#############################

#testBer<-jLeaveOut(dat2model, train=tmpc$City2!="BER", leaveOutName="testBER",tmpc, tlevel, alpha, reg, log0,dir2save, seed0=111, pdfout=T)
#pdf(paste0(dir2save,"/classification_",tlevel,"_",reg,"log_",log0,"_gps_scatter_train_pred_LeftOut_all_",gpsw,".pdf"), width=18, height=12)
#	par(mfrow=c(4,6))


l1oResults<-mclapply(levels(cl), function(city){
	print(city)
	jLeaveOutClass(dat2model, train=tmpc$City2!=city, leaveOutName=city, outc, tlevel, alpha, reg, log0,dir2save, seed0=seed2use, toplot=T)
}, mc.cores=2)
#dev.off()
names(l1oResults)=levels(cl)
#save(l1oResults, file=paste0(dir2save,"/Leave1cityout_results_classification_",tlevel,"_",reg,"log_",log0,".rdata"))


#l1oLong<- ldply(sapply(l1oResults, function(y)y$longitude), data.frame)
l1oDF<-ldply(lapply(l1oResults, function(y)y$gps1), data.frame)
l1oDF$.id<-factor(l1oDF$.id)


l1oDF2<-data.frame(l1oDF)
l1oDF2$gpst= factor(l1oDF2$gpst, level=levels(outc))
l1oDF2$pred= factor(l1oDF2$pred, level=levels(outc))

print(cmvl1o<-confusionMatrix(data= l1oDF2$pred, reference=l1oDF2$gpst))

confusion_matrixl1o <- as.data.frame(cmvl1o$table)

pdf(here:::here(outputFolder,paste0(dir2save,"/classification_",tlevel,"_",reg,"log_",log0,"_scatter_train_pred_LeftOutOnly_all_",gpsw,".pdf")), width=7.5, height=6)
  plotg<-ggplot(data = confusion_matrixl1o ,
       mapping = aes(x = Reference,
                     y = Prediction)) +
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "lavenderblush2",
                      high = "coral2", #lavenderblush3
					  na.value = "lavenderblush",
                      trans = "log")+
  ggtitle("Classification: Leave-one-city-out")
  print(plotg)

dev.off()

### plot in continent cluster
tmpc2<- cbind(City2=as.character(tmpc$City2), City3=paste0(tmpc$ContinentShort,"_", tmpc$City2))
rownames(tmpc2)=tmpc2[,"City2"]
confusion_matrixl1o[,1]=tmpc2[as.character(confusion_matrixl1o[,1]),"City3"]
confusion_matrixl1o[,2]=tmpc2[as.character(confusion_matrixl1o[,2]),"City3"]

pdf(here:::here(outputFolder,paste0(dir2save,"/FigureS2_classification_",tlevel,"_",reg,"log_",log0,"_scatter_train_pred_LeftOutOnly_all_",gpsw,"_contSorted.pdf")), width=7.5, height=6)
  plotg<-ggplot(data = confusion_matrixl1o ,
       mapping = aes(x = Reference,
                     y = Prediction)) +
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "lavenderblush2",
                      high = "coral2", #lavenderblush3
					  na.value = "lavenderblush",
                      trans = "log")+
  ggtitle("Classification: Leave-one-city-out")+
			 theme(axis.text.x = element_text(angle =45, hjust = 1))
  print(plotg)

dev.off()


#cityCol<-cbind(City2=levels(cl),col=cols[1:length(levels(cl))])
#rownames(cityCol)=cityCol[,"City2"]


###############################
#### External 10fold CV, train the rest and predict samples from the left out city
#############################
set.seed(seed2use)
cvcnt=10
cvids<-sample(1:cvcnt, ncol(dat2model), replace=T)
#pdf(paste0(dir2save,"/classification_",tlevel,"_",reg,"log_",log0,"_gps_scatter_pred_10CV_all2_",gpsw,".pdf"), width=18, height=12)
#	par(mfrow=c(4,6))
	CV10Results<-mclapply(1:cvcnt, function(cvid, cvids){
		print(cvid)
		
		jLeaveOutClass(dat2model, train=cvids!=cvid, leaveOutName=paste0("CV",cvid), outc, tlevel, alpha, reg, log0,dir2save, seed0=seed2use, toplot=T, weights=w2balance[cvids!=cvid])

	}, cvids=cvids, mc.cores=2)
	names(CV10Results)=paste0("CV",1:cvcnt)

#dev.off()

#save(CV10Results,cvids, cvcnt, file=paste0(dir2save,"/CV10_results_classification_",tlevel,"_",reg,"log_",log0,".rdata"))


#l1oLong<- ldply(sapply(l1oResults, function(y)y$longitude), data.frame)
CV10DF<-ldply(lapply(CV10Results, function(y)y$gps1), data.frame)
CV10DF$.id<-factor(CV10DF$.id)


CV10DF2<-data.frame(CV10DF[,c(1,2)],pred=apply(CV10DF[,-c(1,2)], 1, which.max))
CV10DF2$gpst= factor(levels(outc)[CV10DF2$gpst])
CV10DF2$pred= factor(levels(outc)[CV10DF2$pred])

print(cmvNestedCV<-confusionMatrix(data= CV10DF2$pred, reference=CV10DF2$gpst))
							   



confusion_matrix <- as.data.frame(cmvNestedCV$table)
###credit to https://stackoverflow.com/questions/37897252/plot-confusion-matrix-in-r-using-ggplot
pdf(here:::here(outputFolder,paste0(dir2save,"/classification_",tlevel,"_",reg,"log_",log0,"_pred_10CV_LeftOutOnly_all_",gpsw,".pdf")), width=7.5, height=6)
plotg<-ggplot(data = confusion_matrix,
       mapping = aes(x = Reference,
                     y = Prediction)) +
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "lavenderblush2",
                      high = "coral2", #lavenderblush3
					  na.value = "lavenderblush",
                      trans = "log")+
  ggtitle(paste("Nested 10-fold CV: Overall accuracy =", signif(cmvNestedCV$overall["Accuracy"],4)))
  print(plotg)
dev.off()

### plot in continent cluster
confusion_matrix[,1]=tmpc2[as.character(confusion_matrix[,1]),"City3"]
confusion_matrix[,2]=tmpc2[as.character(confusion_matrix[,2]),"City3"]

pdf(here:::here(outputFolder,paste0(dir2save,"/Figure4_classification_",tlevel,"_",reg,"log_",log0,"_pred_10CV_LeftOutOnly_all_",gpsw,"_contsorted.pdf")), width=7.5, height=6)
plotg<-ggplot(data = confusion_matrix,
       mapping = aes(x = Reference,
                     y = Prediction)) +
  geom_tile(aes(fill = Freq)) +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "lavenderblush2",
                      high = "coral2", #lavenderblush3
					  na.value = "lavenderblush",
                      trans = "log")+
  ggtitle(paste("Nested 10-fold CV: Overall accuracy =", signif(cmvNestedCV$overall["Accuracy"],4)))+
			 theme(axis.text.x = element_text(angle =45, hjust = 1))
  print(plotg)
dev.off()



#################################################################
############# Diversity / Uncertainty ##########################
###simpson index
#############################
#### computing diversity indices with l1o and cv10 prediction probabilities
require(vegan)
### didn't need to do apply to row for diversity actually, it takes matrix, but same result
l1oShannon<-apply(l1oDF[,!colnames(l1oDF)%in%c(".id","gpst","pred")], 1, function(y){diversity(na.omit(y))})
l1oSimpson<-apply(l1oDF[,!colnames(l1oDF)%in%c(".id","gpst","pred")], 1, function(y){diversity(na.omit(y), "simpson")})


CV10Shannon<-apply(CV10DF[,!colnames(CV10DF)%in%c(".id","gpst","pred")], 1, function(y){diversity(na.omit(y))})
CV10Simpson<-apply(CV10DF[,!colnames(CV10DF)%in%c(".id","gpst","pred")], 1, function(y){diversity(na.omit(y), "simpson")})

wilcox.test(CV10Simpson, l1oSimpson, alternative="less")$p.value #p-value 1.471256e-49
wilcox.test(CV10Shannon, l1oShannon, alternative="less")$p.value #p-value 2.895811e-47



pdf(here:::here(outputFolder,paste0(dir2save,"/Diversity_Fitted_Probabiliy_classification_",tlevel,"_",reg,"log_",log0,".pdf")), width=8.5,height=6)
	par(mfrow=c(2,3))
	hranges <- range(c(CV10Shannon, l1oShannon)*1.1)
	sranges <- range(c(CV10Simpson, l1oSimpson)*1.1)

	hist(CV10Shannon, main="CV10_Shannon", xlim=hranges, freq=F)
	hist( CV10Simpson, main="CV10_Simpson", xlim=sranges, freq=F)
	plot(CV10Shannon,CV10Simpson)
	hist(l1oShannon, main="LeftoutCity_Shannon", xlim=hranges, freq=F)
	hist( l1oSimpson, main="LeftoutCity_Simpson", xlim=sranges, freq=F)
	plot(l1oShannon,l1oSimpson)

dev.off()

pdf(here:::here(outputFolder,paste0(dir2save,"/Diversity_Fitted_Probabiliy_classification_",tlevel,"_",reg,"log_",log0,"_Boxplot.pdf")), width=6,height=6)
boxplot(cbind(CV10=CV10Simpson,L1CO=l1oSimpson), main="Simpson diversity", ylim=c(0,1), col=c("#F8766D","#00BFC4"))
boxplot(cbind(CV10=CV10Shannon,L1CO=l1oShannon), main="Shannon diversity", col=c("#F8766D","#00BFC4"))
dev.off()


datSimp<-data.frame( group=c(rep("CV10",length(CV10Simpson)), rep("L1CO",length(l1oSimpson))),value= c(CV10Simpson, l1oSimpson))
save(CV10Simpson, l1oSimpson, l1oShannon,CV10Shannon,datSimp, file=here:::here(outputFolder,paste0(dir2save,"/Diversity_",tlevel,".rdata")))


pdf(here:::here(outputFolder,paste0(dir2save,"/Diversity_Fitted_Probabiliy_classification_",tlevel,"_",reg,"log_",log0,"_combined.pdf")), width=4.5,height=3.5)
	library(ggplot2)
	p<-ggplot(datSimp, aes(x=value, fill=group)) + geom_histogram(alpha=0.5, position="identity") + xlab("Simpson diversity index")+ theme_bw()+
	ggtitle(paste0("One-sided Wilcox p-value = ",signif(wilcox.test(CV10Simpson, l1oSimpson, alternative="less")$p.value,4)))#+ geom_density(alpha=0.7)
	print(p)

	dat=data.frame( group=c(rep("CV10",length(CV10Shannon)), rep("L1CO",length(l1oShannon))),value= c(CV10Shannon, l1oShannon))
	library(ggplot2)
	p<-ggplot(dat, aes(x=value, fill=group)) + geom_histogram(alpha=0.5, position="identity")+ xlab("Shannon diversity index")+  theme_bw()+
	ggtitle(paste0("One-sided Wilcox p-value = ",signif(wilcox.test(CV10Shannon, l1oShannon, alternative="less")$p.value,4)))#+ geom_density(alpha=0.7)
	print(p)
dev.off()

#exp( l1oResults$TOK$gps1[1,-c(1,ncol(l1oResults$TOK$gps1))])/sum(exp( l1oResults$TOK$gps1[1,-c(1,ncol(l1oResults$TOK$gps1))]))


save.image(here:::here(outputFolder,paste0(dir2save,"/image_",tlevel,".rdata")))
