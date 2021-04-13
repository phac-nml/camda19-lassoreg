#seed2use=111

require(glmnet)
require(metagenomeSeq)
require(plyr)
require(here)

source(here:::here("scripts","CAMDA_jc_source_glmnet.r"))


#### CAMDA data
##### "multivarate gaussian - glmnet" tabs


#outputFolder="output"
#outputFolder="output_pairedsingle"


load(here:::here(outputFolder,paste0("bracken_",tlevel,".rdata")))
#load(here:::here(outputFolder,"bracken_species.rdata"))
#load(here:::here(outputFolder,"bracken_genus.rdata"))
#load(here:::here(outputFolder,"bracken_family.rdata"))


######################################
### multivariate gaussian
######################################
cities<-read.table(here:::here("data","_Cities.txt"), sep="\t", header=T, quote="")
tmpc<-merge(data.frame(pData(metaobj)), cities, by.x="City2",by.y="Code",all.x=T, sort=F)

gps<-as.matrix(tmpc[,c("Longitude","Latitude")])

dir2save="glmnet_MSE"

if(!dir.exists(here:::here(outputFolder,dir2save))){dir.create(here:::here(outputFolder,dir2save))}

log0=T
norm0=T


print(tlevel)
dat2model<-MRcounts(metaobj, norm=norm0, log=log0)
rownames(dat2model)=fData(metaobj)[,tlevel]
dat2model<-dat2model[rownames(dat2model)!="",]  
rownames(dat2model)[grep("^[0-9]",rownames(dat2model))]=paste0("zzz-",rownames(dat2model)[grep("^[0-9]",rownames(dat2model))])

if(dir2save=="glmnet_MSE_euro"){

	whichEuro<-tmpc$Continent=="Europe"
	dat2model<-dat2model[,whichEuro]
	tmpc<-tmpc[whichEuro,]
	gps<-gps[whichEuro,]
	
}


##########################################
###### GLMNET MSE
##########################################

alpha=1;reg="lasso"
set.seed(seed2use) 
system.time(cvob2mse<-cv.glmnet(t(dat2model) ,gps, family="mgaussian", type.measure="mse", nfolds =10, alpha=alpha, parallel=T, keep=T) ) 

##########################################
###### plotting the predicted values vs true in cross validation
##########################################
dim(cvob2mse$fit.preval) # 3043    2  100
cvPred1<-cvob2mse$fit.preval[,1,which(cvob2mse$lambda==cvob2mse$lambda.1se)]
cvPred2<-cvob2mse$fit.preval[,2,which(cvob2mse$lambda==cvob2mse$lambda.1se)]
pdf(here:::here(outputFolder,paste0(dir2save,"/Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_scatter_CV_pred_lambda.1se_MSE_Loss2.pdf")), height=14)
	par(mfrow=c(2,1))
	a=gps[,1]; b=cvPred1
	plot(a,b, xlab="True longitude",ylab="CV Pred",xlim=range(c(a,b)), ylim=range(c(a,b)), main=paste("longitude", tlevel,reg, "r2=",signif(cor(a,b)^2,4)))
	abline(lm(b~a)) 
	
	a=gps[,2]; b=cvPred2
	plot(a,b, xlab="True latitude",ylab="CV Pred",xlim=range(c(a,b)), ylim=range(c(a,b)),main=paste("latitude", tlevel,reg,"r2=", signif(cor(a,b)^2,4)))
	abline(lm(b~a)) 
dev.off()
#cvob2mse$glmnet.fit


fit1m=glmnet(t(dat2model) ,gps, family="mgaussian", alpha=alpha) 
fit1m2=glmnet(t(dat2model) ,gps, family="mgaussian", lambda=cvob2mse$lambda.min, alpha=alpha) 
fit1m3=glmnet(t(dat2model) ,gps, family="mgaussian", lambda=cvob2mse$lambda.1se, alpha=alpha) 

save(cvob2mse, fit1m, fit1m3, tlevel, reg,log0, file=here:::here(outputFolder,paste0(dir2save,"/Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_MSE.rdata")))



which1se<-which(round(cvob2mse$glmnet.fit$lambda,6)==round(cvob2mse$lambda.1se,6))
df1se<-cvob2mse$glmnet.fit$df[which1se]
cvm1se<-cvob2mse$cvm[which1se]
print(paste("lambda 1se: ",cvob2mse$lambda.1se, cvob2mse$lambda[which1se]))
print(paste("$cvm MSE value",cvm1se))
print(paste("df",df1se))


##########################################
###### plotting the model info and lambda estimation
##########################################
pdf(here:::here(outputFolder,paste0(dir2save,"/Multivariate_",tlevel,"_",reg,"log_",log0,"_gps.pdf")))
	plot(fit1m,type.coef="2norm")
	fit1mRange=range(log(fit1m$lambda))
	plot(fit1m, xvar = "lambda", label = TRUE, type.coef = "2norm", xlim=c(fit1mRange[1]-0.7, fit1mRange[2]))
	abline(v=log(cvob2mse$lambda.min), lty=2)
	abline(v=log(cvob2mse$lambda.1se), lty=2)
	#plot(fit1m2, xvar = "lambda", label = TRUE, type.coef = "2norm", main="lambda.min")
	#plot(fit1m3, xvar = "lambda", label = TRUE, type.coef = "2norm", main="lambda.1se")
	
	par(mar=c(4.5,4.5,4,1)) 
	plot(cvob2mse)
	title(paste("Multivariate Gaussian", tlevel, reg),line=2.5)
	#plot(cvob2mse, xvar='lambda', label=TRUE)
dev.off()

##########################################
##### plot top 30 features
##########################################
coef1<-coef(cvob2mse)[[1]][,1][coef(cvob2mse)[[1]][,1]!=0][-1] ## skipping the intercept
coef1<-coef1[order(abs(coef1),decreasing=T)][1:30]
coef1<-sort(coef1)

coef2<-coef(cvob2mse)[[2]][,1][coef(cvob2mse)[[2]][,1]!=0][-1] ## skipping the intercept
coef2<-coef2[order(abs(coef2),decreasing=T)][1:30]
coef2<-sort(coef2)

pdf(here:::here(outputFolder,paste0(dir2save,"/Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_featImp.pdf")), height=10)
	par(mar=c(4.5,10,4,1)) 
	barplot(coef1, horiz=T,las=2, main="longitude top30", cex.names=0.9)
	barplot(coef2, horiz=T,las=2, main="latitude top30", cex.names=0.9)
dev.off()


##########################################
### examine abundance of top taxa along lat and long,color
##########################################
 pdf(here:::here(outputFolder,paste0(dir2save,"/Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_featImp_abundance.pdf")), width=12,height=7)
for(taxa in unique(c(rev(names(coef1)),names(coef2)))){

	### some share the same species names but have different genera, need to edit the script to keep Genera
	for (taxai in which(rownames(dat2model)==taxa)){
		
		tmp<-as.vector(dat2model[taxai,])
		
		dat2gg<-cbind(gps, tmp)
		colnames(dat2gg)=c("longitude","latitude","abundance")
		
			#library(ggplot2)
			#jitter <- position_jitter(width = 2, height =2)
			#sp<-ggplot(dat2gg, aes(x=longitude, y=latitude, color=abundance)) + geom_point(shape=1,alpha=0.3,position = jitter) 
			#sp+scale_color_gradient(low="blue", high="red")+ggtitle(taxa)
		
		par(mfrow=c(1,2), mar=c(4,4,2,2))
			#plot(jitter(dat2gg[,1]),jitter(dat2gg[,3], 0.1))
			#plot(jitter(dat2gg[,2]),jitter(dat2gg[,3], 0.1))
		plot(dat2gg[,c(1,3)], main=paste(tlevel,taxa))
		plot(dat2gg[,c(2,3)], main=paste(tlevel,taxa))
	}
}
dev.off()

		
		
		
####		
####		
#### use jLeaveOut function onwards
####		
####
pch2use=c(1: length(levels(cl2)))
pData2<-pData(metaobj)
continents<-pData2[,"contin", drop=F]



###############################
#### leave one city out, train the rest and predict samples from the left out city
#### here tends to break below
#############################

#testBer<-jLeaveOut(dat2model, train=tmpc$City2!="BER", leaveOutName="testBER",tmpc, tlevel, alpha, reg, log0,dir2save,  typemeasure="mse", seed0=111, pdfout=T)
pdf(here:::here(outputFolder,paste0(dir2save,"/Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_scatter_train_pred_LeftOut_all.pdf")), width=18, height=12)
	par(mfrow=c(4,6))
	l1oResults<-lapply(levels(cl), function(city){
		print(city)
		
		jLeaveOut(dat2model, train=tmpc$City2!=city, leaveOutName=city, tmpc, tlevel, alpha, reg, log0,dir2save,  typemeasure="mse", seed0=seed2use, toplot=T)

	})
dev.off()
names(l1oResults)=levels(cl)


l1osamples<-unlist(lapply(l1oResults, function(y){rownames(y[[1]])})) ## getting sample names for coloring
l1ocitycode <- do.call(rbind,strsplit(l1osamples,"_|\\."))[,2]
l1ocontinents<-continents[l1osamples,1]
l1ocontinentShape<-pch2use[l1ocontinents]
#l1oLong<- ldply(sapply(l1oResults, function(y)y$longitude), data.frame)
l1oDF<-ldply(l1oResults, data.frame)
l1oDF$.id<-factor(l1oDF$.id)


cityCol<-cbind(City2=levels(cl),col=cols[1:length(levels(cl))])
rownames(cityCol)=cityCol[,"City2"]

pdf(here:::here(outputFolder,paste0(dir2save,"/FigureS2_Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_scatter_train_pred_LeftOutOnly_all.pdf")), width=7, height=14)
	par(mfrow=c(2,1), cex=1.2)
	a=l1oDF[,"longitude.1"]; b=l1oDF[,"longitude.2"]
	plot(a,b, xlab="True Longitude",ylab="Leave-1-city-out prediction", xlim=range(c(a,b)), ylim=range(c(a,b)),  main=paste("Longitude", tlevel, reg, "r2=",round(cor(a,b)^2,4)), col=cityCol[as.character(l1oDF$.id),"col"], pch=l1ocontinentShape)
	#legend("bottom",levels(l1oDF$.id),col=cityCol[levels(l1oDF$.id),"col"], pch=19,ncol=6) 
	abline(a=0,b=1, lty=2)
	a=l1oDF[,"latitude.1"]; b=l1oDF[,"latitude.2"]
	plot(a,b, xlab="True Latitude",ylab="Leave-1-city-out prediction", xlim=range(c(a,b)), ylim=range(c(a,b)), main=paste("Latitude", tlevel, reg, "r2=",round(cor(a,b)^2,4)), col=cityCol[as.character(l1oDF$.id),"col"], pch=l1ocontinentShape)
	#legend("bottom",levels(l1oDF$.id),col=cityCol[levels(l1oDF$.id),"col"],pch=pch2use[levels(cl2)], ncol=6,cex = 0.8) 
	legend("bottom",levels(cl),col=cityCol[levels(cl),"col"],pch=pch2use[clplus], ncol=6,cex = 0.8)
	abline(a=0,b=1, lty=2)
dev.off()


seAdd<-function(vdat){
	toret<-data.frame(vdat, longSE=(vdat$longitude.1-vdat$longitude.2)^2,
		latSE=(vdat$latitude.1-vdat$latitude.2)^2)
	toret<-data.frame(toret,llSE=toret$longSE+toret$latSE)
	toret
}
l1oDF<-seAdd(l1oDF)
plot(l1oDF$longSE, l1oDF$latSE)
mean(l1oDF$longSE);mean(l1oDF$latSE)

##### Any cities had higher SE?
pdf(here:::here(outputFolder,paste0(dir2save,"/Multivariate_",tlevel,"_",reg,"log_",log0,"_Boxplot_leave1Out_squarederror.pdf")), width=8, height=9)
require(ggpubr)
par(mfrow=c(3,1))
gb1<-ggboxplot(l1oDF, ".id", "longSE", color = ".id",  
		add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  
		ylab = "Longitude Squared Error", xlab="", main="Leave-1-city-out: Longitude Squared Error")+ scale_color_manual(values=cityCol[levels(l1oDF$.id),"col"])+ theme(legend.position="none")
gb2<-ggboxplot(l1oDF, ".id", "latSE", color = ".id",  
		add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  
		ylab = "Latitude Squared Error", xlab="", main="Leave-1-city-out: Latitude Squared Error")+ scale_color_manual(values=cityCol[levels(l1oDF$.id),"col"])+ theme(legend.position="none")
gb3<-ggboxplot(l1oDF, ".id", "llSE", color = ".id",  
		add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  
		ylab = "Squared Error", xlab="", main="Leave-1-city-out: Longitude+Latitude Squared Error")+ scale_color_manual(values=cityCol[levels(l1oDF$.id),"col"])+ theme(legend.position="none")
ggarrange(gb1, gb2, gb3, 
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)+ scale_color_manual(values=cityCol[levels(l1oDF$.id),"col"])+ theme(legend.position="none")
dev.off()



###############################
#### External 10fold CV, train the rest and predict samples from the left out city
#############################

set.seed(seed2use)
cvcnt=10
cvids<-sample(1:cvcnt, ncol(dat2model), replace=T)
pdf(here:::here(outputFolder,paste0(dir2save,"/Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_scatter_pred_10CV_all2.pdf")), width=18, height=12)
	par(mfrow=c(4,6))
	CV10Results<-lapply(1:cvcnt, function(cvid, cvids){
		print(cvid)
		
		jLeaveOut(dat2model, train=cvids!=cvid, leaveOutName=paste0("CV",cvid), tmpc, tlevel, alpha, reg, log0,dir2save,  typemeasure="mse", seed0=seed2use, toplot=T)

	}, cvids=cvids)
	names(CV10Results)=paste0("CV",1:cvcnt)

dev.off()

cvids2<-cbind(order=1:length(cvids), cvids)

CV10samples<-unlist(lapply(CV10Results, function(y){rownames(y[[1]])})) ## getting sample names for coloring
CV10citycode <- do.call(rbind,strsplit(CV10samples,"_|\\."))[,2]
CV10Col<- cityCol[CV10citycode,2] ## getting color of the city from sample names
CV10continents<-continents[CV10samples,1]
CV10continentShape<-pch2use[CV10continents]

#l1oLong<- ldply(sapply(l1oResults, function(y)y$longitude), data.frame)
CV10DF<-ldply(CV10Results, data.frame)
CV10DF$.id<-factor(CV10DF$.id)



pdf(here:::here(outputFolder,paste0(dir2save,"/Figure4_Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_scatter_pred_10CV_LeftOutOnly_all.pdf")), width=7, height=14)
	par(mfrow=c(2,1), cex=1.2)
	a=CV10DF[,"longitude.1"]; b=CV10DF[,"longitude.2"]
	plot(a,b, xlab="True Longitude",ylab="Nested 10-fold CV prediction", xlim=range(c(a,b)), ylim=range(c(a,b)),  main=paste("Longitude", tlevel, reg, "r2=",round(cor(a,b)^2,4)), col=CV10Col, pch=CV10continentShape)
	abline(a=0,b=1, lty=2)
	a=CV10DF[,"latitude.1"]; b=CV10DF[,"latitude.2"]
	plot(a,b, xlab="True Latitude",ylab="Nested 10-fold CV prediction", xlim=range(c(a,b)), ylim=range(c(a,b)), main=paste("Latitude", tlevel, reg, "r2=",round(cor(a,b)^2,4)), col=CV10Col, pch=CV10continentShape)
	abline(a=0,b=1, lty=2)
	legend("bottom",levels(cl),col=cityCol[levels(cl),"col"],pch=pch2use[clplus], ncol=6,cex = 0.8)
dev.off()


CV10DF<-seAdd(CV10DF)
CV10DF<-data.frame(CV10DF, CV10samples,CV10citycode, CV10Col, CV10continents, CV10continentShape, stringsAsFactors = F)

mean(CV10DF$longSE);mean(CV10DF$latSE)
plot(CV10DF$longSE, CV10DF$latSE, col=CV10DF$CV10Col, pch=CV10DF$CV10continentShape)

### don't really need this
pdf(here:::here(outputFolder,paste0(dir2save,"/Multivariate_",tlevel,"_",reg,"log_",log0,"_Boxplot_10CV_squarederror.pdf")), width=8, height=9)
  require(ggpubr)
  par(mfrow=c(3,1))
  gb1<-ggboxplot(CV10DF, ".id", "longSE", #color = "master",  
  		add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  
  		ylab = "Longitude Squared Error", xlab="", main="10CV: Longitude Squared Error")
  gb2<-ggboxplot(CV10DF, ".id", "latSE", #color = "master",  
  		add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  
  		ylab = "Latitude Squared Error", xlab="", main="10CV: Latitude Squared Error")
  gb3<-ggboxplot(CV10DF, ".id", "llSE", #color = "master",  
  		add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  
  		ylab = "Squared Error", xlab="", main="10CV: Longitude+Latitude Squared Error")
  ggarrange(gb1, gb2, gb3, 
            labels = c("A", "B", "C"),
            ncol = 1, nrow = 3)
dev.off()


save(l1oResults,l1oDF,cvids,CV10DF, CV10Results, file=here:::here(outputFolder,paste0(dir2save,"/Test_Validation_l1o_CV10_",tlevel,"_",reg,"log_",log0,".rdata")))
#save(l1oResults,l1oDF,cvids,CV10DF, CV10Results, file=paste0("/Drives/K/chichen/CAMDA/CAMDA_Rproject/CAMDA_DS/output/","Test_Validation_l1o_CV10_",tlevel,"_",reg,"log_",log0,".rdata"))


paste("Leave1cityout MSE",round(mean(l1oDF$llSE),2), ";",
"10-fold nested MSE",round(mean(CV10DF$llSE),2) )
#SecondPass_paired/glmnet_MSE "Leave1cityout MSE 5962.3 ; 10-fold nested MSE 839.9"
#ThirdPass_12/glmnet_MSE "Leave1cityout MSE 7370.12 ; 10-fold nested MSE 888.34"

save.image(here:::here(outputFolder,paste0(dir2save,"/image_",tlevel,".rdata")))