require(glmnet)
require(geosphere)

####### functions


###############################
#### General functions: leave samples out, train the rest and predict leftout samples
#############################
### train: T for train, F for test
jLeaveOut<-function(dat2model00, train, leaveOutName="",locInfo, tlevel, alpha, reg, log0,dir2save,  typemeasure="mse", seed0=111, pdfout=F,toplot=F, saveModel=F, nfolds =10, ...){
	
	gps00<-as.matrix(locInfo[,c("Longitude","Latitude")])
	
	dat2model0 <- dat2model00[,train]
	gps0<-gps00[train,]
	
	dat2test <- dat2model00[,!train]
	gpst<-gps00[!train,]

	
	set.seed(seed0) 
	system.time(cvob2Hav0<-cv.glmnet(t(dat2model0) ,gps0, family="mgaussian", type.measure=typemeasure, nfolds = nfolds, alpha=alpha, parallel=T, keep=T, ...)) 

	##########################################
	###### plotting the predicted values vs true in cross validation
	##########################################
	dim(cvob2Hav0$fit.preval) # 3043    2  100
	cvPred1<-cvob2Hav0$fit.preval[,1,which(cvob2Hav0$lambda==cvob2Hav0$lambda.1se)]
	cvPred2<-cvob2Hav0$fit.preval[,2,which(cvob2Hav0$lambda==cvob2Hav0$lambda.1se)]
		#pdf(paste0(dir2save,"/Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_scatter_CV_pred_lambda.1se_HaversineDistLoss2.pdf"), height=14)
		#	par(mfrow=c(2,1))
		#	a=gps0[,1]; b=cvPred1
		#	plot(a,b, xlab="True longitude",ylab="CV Pred",xlim=range(c(a,b)), ylim=range(c(a,b)), main=paste("longitude", tlevel,reg, "r2=",signif(cor(a,b)^2,4)))
		#	abline(lm(b~a)) 
			
		#	a=gps0[,2]; b=cvPred2
		#	plot(a,b, xlab="True latitude",ylab="CV Pred",xlim=range(c(a,b)), ylim=range(c(a,b)),main=paste("latitude", tlevel,reg,"r2=", signif(cor(a,b)^2,4)))
		#	abline(lm(b~a)) 
		#dev.off()
	#cvob2Hav0$glmnet.fit


	fit1m=glmnet(t(dat2model0) ,gps0, family="mgaussian", alpha=alpha, ...) 
	#fit1m2=glmnet(t(dat2model0) ,gps0, family="mgaussian", lambda=cvob2Hav0$lambda.min, alpha=alpha) 
	fit1m3=glmnet(t(dat2model0) ,gps0, family="mgaussian", lambda=cvob2Hav0$lambda.1se, alpha=alpha, ...) 

	if(saveModel){
		save(cvob2Hav0, fit1m, fit1m3, tlevel, reg, log0, file=here(outputFolder,paste0("Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_Haversine_LeftOut_",leaveOutName,".rdata")))
	}

	##########################################
	##### making predictions: currently on training data, use *_gps_scatter_CV_pred_lambda.1se.pdf instead
	### Default is the value s="lambda.1se"
	### can be used to predict other data later on
	##########################################
	predGpst<-predict(fit1m3,newx=t(dat2test))  ### lambda.1se on full data  #dim 3043    2    1
	if(pdfout) toplot=T
	if(pdfout){
		pdf(here(outputFolder,paste0("Multivariate_",tlevel,"_",reg,"log_",log0,"_gps_scatter_train_pred_LeftOut_",leaveOutName,".pdf")))
	}
	if(toplot){
			a=c(gpst[,1],gps0[,1]); b=c(predGpst[,1,1],cvPred1)
			plot(a,b, xlab="True",ylab="TrainPred",xlim=range(c(a,b)), ylim=range(c(a,b)), col=c(rep("blue",nrow(gpst)), rep("darkgrey",nrow(gps0))), pch=c(rep(19,nrow(gpst)), rep(1,nrow(gps0))), main=paste(leaveOutName,"longitude ",reg,"MSE Haversine" )) #"r2=",signif(cor(a,b)^2,4)
			abline(a=0,b=1) 
			a=c(gpst[,2],gps0[,2]); b=c(predGpst[,2,1],cvPred2)
			plot(a,b, xlab="True",ylab="TrainPred",xlim=range(c(a,b)), ylim=range(c(a,b)), col=c(rep("blue",nrow(gpst)), rep("darkgrey",nrow(gps0))), pch=c(rep(19,nrow(gpst)), rep(1,nrow(gps0))), main=paste(leaveOutName,"latitude",reg,"MSE Haversine"))
			abline(a=0,b=1) 
	}
	if(pdfout){
		dev.off()
	}
		list(longitude=cbind(gpst[,1],predGpst[,1,1]), latitude=cbind(gpst[,2],predGpst[,2,1]))
}



	### toy example
	#multivariate gaussian 
	sampleExp<-function(){
		x=matrix(rnorm(100*20),100,20) 
		y=matrix(rnorm(100*2),100,2) 
		fit1m=glmnet(x,y,family="mgaussian") 
		plot(fit1m,type.coef="2norm")

		predict(fit1m,newx=x[1:10,],s=c(0.01,0.005))
		y[1:10,]

		i=1
		plot(predict(fit1m,newx=x[1:10,],s=c(0.01,0.005))[i,,2], y[i,])

		fit1mcv=cv.glmnet(x,y,family="mgaussian", type.measure="mse")
		fit1mcv3=cv.glmnet(x,y,family="mgaussian",type.measure="hav")
	}

	


###############################
#### General functions: leave samples out, train the rest and predict leftout samples
#### for classification
#############################
### train: T for train, F for test
jLeaveOutClass<-function(dat2model00, train, leaveOutName="",outcome1, tlevel, alpha, reg, log0,dir2save, seed0=111, pdfout=F,toplot=F, saveModel=F, nfolds = 10, ...){
	
	gps00<-factor(outcome1)
	
	dat2model0 <- dat2model00[,train]
	#gps0<-droplevels(gps00[train])
	#gps0<-factor(gps00[train], level=levels(gps00))
	if(sum(table(gps00)==0)){
		gps0<-gps00[train]
	}else{gps0<-droplevels(gps00[train])}
	
	dat2test <- dat2model00[,!train]
	#gpst<-factor(gps00[!train], level=levels(gps00))
	gpst<-gps00[!train]

	
	set.seed(seed0) 
	system.time(cvob2Hav0<-cv.glmnet(t(dat2model0) ,gps0, family="multinomial", type.multinomial = "grouped", nfolds =nfolds, alpha=alpha, parallel=T, keep=T, ...)) 

	##########################################
	###### plotting the predicted values vs true in cross validation
	##########################################
	dim(cvob2Hav0$fit.preval) # 3043    2  100

	fit1m=glmnet(t(dat2model0) ,gps0, family="multinomial", type.multinomial = "grouped", alpha=alpha, ...) 
	#fit1m2=glmnet(t(dat2model0) ,gps0, family="multinomial", type.multinomial = "grouped", lambda=cvob2Hav0$lambda.min, alpha=alpha) 
	fit1m3=glmnet(t(dat2model0) ,gps0, family="multinomial", type.multinomial = "grouped", lambda=cvob2Hav0$lambda.1se, alpha=alpha, ...) 

	if(saveModel){ 
		save(cvob2Hav0, fit1m, fit1m3, tlevel, reg, log0, file=here(outputFolder,paste0(dir2save,"/classification_",tlevel,"_",reg,"log_",log0,"__LeftOut_",leaveOutName,".rdata")))
	}

	##########################################
	##### making predictions: currently on training data, use *_gps_scatter_CV_pred_lambda.1se.pdf instead
	### Default is the value s="lambda.1se"
	### can be used to predict other data later on
	##########################################
	predGpst<-predict(fit1m3,newx=t(dat2test), type="response")  ### lambda.1se on full data  #dim 3043    2    1   ### type edited 19-11-06
	predMax<-factor(levels(gps0)[apply(predGpst[,,1], 1,which.max)], level=levels(gps0))
	if(length(levels(predMax))==length(levels(gpst))){
		print(cmv<-confusionMatrix(data= predMax,  
						reference=gpst))
	}else{
		cmv=NULL
	}
	### for multi-class comparisons		
	#print("Overall accuracies")
	#overallacc<-cmv$overall["Accuracy"]
	#print(overallacc)
	#print(cmv$byClass)
	#apply(cmv$byClass, 2, median)

	
	list(gps1=data.frame(gpst,predGpst[,,1],pred=predMax), cmv=cmv)
}

	
	