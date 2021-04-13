##### Run CAMDA_jc_mysteryPred.r first for both Mystery paired-end and single end data, then run scripts below
##### keeping results in output_mystery folder

library(here)
require(ggpubr)
require(measurements)
load(here:::here("output_mystery", paste0("Mystery_New_",tlevel,".rdata")))

log0=T
norm0=T

#tlevel="species"
#tlevel="genus"

#modelDir="glmnet_MSE"
#modelDir="glmnet_MSE_pairedsingle"
#modelDir="glmnet_Classification"
#modelDir="glmnet_Classification_pairedsingle"
#modelDir="RFClassification"

mfig=""
if(modelDir%in%c("glmnet_MSE_pairedsingle", "glmnet_Classification_weighted_pairedsingle")){mfig="FigS3"}
if(modelDir%in%c("glmnet_MSE", "glmnet_Classification_weighted")) {mfig="Fig5"}

predp<-read.table(here:::here("output_mystery", paste0("Prediction_",modelDir,"_",tlevel,"TRUETRUE.txt")), sep="\t", header=T, quote = "", as.is=T)
preds<-read.table(here:::here("output_mystery_single",paste0("Prediction_",modelDir,"_",tlevel,"TRUETRUE.txt")), sep="\t", header=T, quote = "", as.is=T)

predpInfo<- do.call(rbind,lapply(strsplit(rownames(predp),"_"),function(y)y[1:3]))
colnames(predpInfo)=c("Batch","ID","NULL")

predsInfo<- do.call(rbind,lapply(strsplit(rownames(preds),"_"),function(y)y[1:3]))
colnames(predsInfo)=c("Batch","ID","pair")


#### saving the merged prediction for paired-end and single-end mystery samples
predMerge<-data.frame(rbind(preds,predp),rbind(predsInfo, predpInfo))
#finalMeanPred<-aggregate(predMerge[,1:2],by=list(predMerge[,4]), FUN=mean)

###
finalPred<-predMerge[is.na(predMerge$pair),]  ### take the prediction on pair-end mystery data when available, if not, using prediction on single-end mystery data

write.table(finalPred,
	file=here:::here("output_mystery",paste0("Prediction_",modelDir,"_",tlevel,"_2and1.txt")),
	col.names=T, row.names=F, quote=F,sep="\t")
	
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
save(truthAll, file=here:::here("output_mystery","MysterySamples_Labels_3.rdata"))
#write.table(truthAll,file=paste0("MysterySamples_Labels_3.txt"),col.names=T, row.names=F, quote=F,sep="\t")
	
#mot<-merge(finalMeanPred, truthAll, by.x="Group.1",by.y="sample")  ### using mean of pred on singles and paired
mot<-merge(finalPred, truthAll, by.x="ID",by.y="sample") ## just use the pred from paired if available, otherwise, single

seAdd2<-function(vdat){
	toret<-data.frame(vdat, longSE=(vdat$Longitude.x-vdat$Longitude.y)^2,
		latSE=(vdat$Latitude.x-vdat$Latitude.y)^2)
	toret<-data.frame(toret,llSE=toret$longSE+toret$latSE)
	toret
}


		###Univariate lat: 
		if(length(grep("_Latitude_",modelDir)!=0)!=0){
			pdf(here:::here("output_mystery",paste0("MysterySample_Prediction_",modelDir,"_",tlevel,"_log_",log0,"_norm",norm0,"_result.pdf")), width=7, height=7)
				a=mot$Latitude; b=mot$s0
				plot( a,b, main=modelDir, xlim=range(c(a,b)), ylim=range(c(a,b))) 
				abline(a=0,b=1, lty=2)
			dev.off()	
			seAddUni<-function(vdat){
				toret<-data.frame(vdat,latSE=(vdat$s0-vdat$Latitude)^2)
				toret
			}
			mot<-seAddUni(mot)
			mot$City2=sapply(strsplit(as.character(mot$City)," \\("), function(y)y[1])
			mcities<-read.table(here:::here("data","_Cities_Mystery.txt"), sep="\t", header=T, quote="")

			### order according to levels of the city
			mot<-mot[order(as.numeric( mot$City)),]
			mean(mot$latSE) # 1206.179
			save.image(here:::here("output_mystery",paste0("MysterySample_Prediction_overall_",modelDir,"_",tlevel,norm0,log0,".rdata")))
		}
		if(length(grep("_LatitudeAbs_",modelDir)!=0)!=0){
			pdf(here:::here("output_mystery",paste0("MysterySample_Prediction_",modelDir,"_",tlevel,"_log_",log0,"_norm",norm0,"_result.pdf")), width=7, height=7)
				a=abs(mot$Latitude); b=mot$s0
				plot( a,b, main=modelDir, xlim=range(c(a,b)), ylim=range(c(a,b))) 
				abline(a=0,b=1, lty=2)
			dev.off()	
			seAddUni<-function(vdat){
				toret<-data.frame(vdat,latSE=(vdat$s0-abs(vdat$Latitude))^2)
				toret
			}
			mot<-seAddUni(mot)
			mot$City2=sapply(strsplit(as.character(mot$City)," \\("), function(y)y[1])
			mcities<-read.table(here:::here("data","_Cities_Mystery.txt"), sep="\t", header=T, quote="")
			
			### order according to levels of the city
			mot<-mot[order(as.numeric( mot$City)),]
			mean(mot$latSE) #319
			save.image(here:::here("output_mystery",paste0("MysterySample_Prediction_Overall_",modelDir,"_",tlevel,norm0,log0,".rdata")))
		}
				

#####################################
### moment of truth  - regression
####################################
#### x and y corrected 20190912
if(length(grep("MSE", modelDir))==1){
  pdf(here:::here("output_mystery",paste0(mfig,"ab_MysterySample_Prediction_",modelDir,"_",tlevel,"_log_",log0,"_norm",norm0,"_result.pdf")), width=7, height=14)
  	par(mfrow=c(2,1), cex=1.2)
  		
  	a=mot$Longitude.y; b=mot$Longitude.x
  	plot( a,b, xlab="True Longitude", ylab="Prediction", main=paste("Longitude of Mystery samples", "r2=",round(cor(a,b)^2,4)), xlim=range(c(a,b)), ylim=range(c(a,b)))
  	abline(a=0,b=1, lty=2)
  	#abline(lm(b~a))
  
  	a=mot$Latitude.y; b=mot$Latitude.x
  	plot( a,b, xlab="True Latitude", ylab="Prediction", main=paste("Latitude of Mystery samples", "r2=",round(cor(a,b)^2,4)), xlim=range(c(a,b)), ylim=range(c(a,b)))
  	abline(a=0,b=1, lty=2)
  
  dev.off()
  
  
  ##### Any cities had higher SE?
  
  mot<-seAdd2(mot)
  mot$City2=sapply(strsplit(as.character(mot$City)," \\("), function(y)y[1])
  mcities<-read.table(here:::here("data","_Cities_Mystery.txt"), sep="\t", header=T, quote="")
  
  ### order according to levels of the city
  mot<-mot[order(as.numeric( mot$City)),]
  
  
  pdf(here:::here("output_mystery",paste0(mfig,"de_MysterySample_Prediction_",modelDir,"_",tlevel,"_","log_",log0,"_norm",norm0,"_Boxplot_testDat_squarederror.pdf")), width=6, height=10)
  
  	par(mfrow=c(3,1))
  	gb1<-ggboxplot(mot, "City2", "longSE", #color = "master",  
  			add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  outlier.shape = NA,
  			ylab = "Longitude Squared Error", xlab="", main="Mystery samples: Longitude Squared Error")+
  			 theme(axis.text.x = element_text(angle =45, hjust = 1))
  	gb2<-ggboxplot(mot, "City2", "latSE", #color = "master",  
  			add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  outlier.shape = NA,
  			ylab = "Latitude Squared Error", xlab="", main="Mystery samples: Latitude Squared Error")+
  			 theme(axis.text.x = element_text(angle =45, hjust = 1))
  	gb3<-ggboxplot(mot, "City2", "llSE", #color = "master",  
  			add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  outlier.shape = NA,
  			ylab = "Squared Error", xlab="", main="Mystery samples: Longitude+Latitude Squared Error")+
  			 theme(axis.text.x = element_text(angle =45, hjust = 1))
  	gbg<- ggarrange(gb1, gb2, gb3, labels = c("A", "B", "C"), ncol = 1, nrow = 3)
  	print(gbg)
  dev.off()
  
  mean(mot$latSE)
  mean(mot$longSE)
  mean(mot$llSE)
  save.image(here:::here("output_mystery",paste0("MysterySample_Prediction_overall_",modelDir,"_",tlevel,norm0,log0,".rdata")))
}
	
#####################################
### moment of truth  - classification
####################################

if(length(grep("Classification", modelDir))==1){
    
  cities<-read.table(here:::here("data","_Cities.txt"), sep="\t", header=T, quote="")
  cities$City3=paste(cities$ContinentShort, cities$Code,sep="_")
  classCol=ifelse(sum(colnames(mot)%in%"s0")==1,"s0", "x") ## s0 from glmnet, x frmo rf
  
  mot<-merge(mot, cities, by.x=classCol,by.y="Code", sort=F)
  
  
  
  confusionTest <- as.data.frame(table(mot$City3, mot$City.x))
  colnames(confusionTest)[1:2]=c("Prediction","Reference")
  
  pdf(here:::here("output_mystery",paste0(mfig,"c_MysterySample_Prediction_",modelDir,"_",tlevel,"_log_",log0,"_norm",norm0,"_result.pdf")), width=6, height=5)
    plotg<-ggplot(data = confusionTest ,
         mapping = aes(x = Reference,
                       y = Prediction)) +
    geom_tile(aes(fill = Freq)) +
    geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
    scale_fill_gradient(low = "lavenderblush2",
                        high = "coral2", #lavenderblush3
  					  na.value = "lavenderblush",
                        trans = "log")+
    ggtitle("Classification: Mystery Samples")+
  			 theme(axis.text.x = element_text(angle =45, hjust = 1))
    print(plotg)
  
  dev.off()
  
  mot<-seAdd2(mot)
  mot$City2=sapply(strsplit(as.character(mot$City.x)," \\("), function(y)y[1])
  ### order according to levels of the city
  mot<-mot[order(as.numeric( mot$City.x)),]
  
  
  pdf(here:::here("output_mystery",paste0(mfig,"fg_MysterySample_Prediction_",modelDir,"_",tlevel,"_","log_",log0,"_norm",norm0,"_Boxplot_testDat_squarederror.pdf")), width=6, height=10)

  	par(mfrow=c(3,1))
  	gb1<-ggboxplot(mot, "City2", "longSE", #color = "master",  
  			#add = "jitter",add.params = list(size = 0.2, jitter = 0.25),  outlier.shape = NA,
  			ylab = "Longitude Squared Error", xlab="", main="Mystery samples: Longitude Squared Error")+
  			 theme(axis.text.x = element_text(angle =45, hjust = 1))
  	gb2<-ggboxplot(mot, "City2", "latSE", #color = "master",  
  			#add = "jitter",add.params = list(size = 0.2, jitter = 0.25),   outlier.shape = NA,
  			ylab = "Latitude Squared Error", xlab="", main="Mystery samples: Latitude Squared Error")+
  			 theme(axis.text.x = element_text(angle =45, hjust = 1))
  	gb3<-ggboxplot(mot, "City2", "llSE", #color = "master",  
  			#add = "jitter",add.params = list(size = 0.2, jitter = 0.25),   outlier.shape = NA,
  			ylab = "Squared Error", xlab="", main="Mystery samples: Longitude+Latitude Squared Error")+
  			 theme(axis.text.x = element_text(angle =45, hjust = 1))

  	gbg<- ggarrange(gb1, gb2, gb3, labels = c("A", "B", "C"), ncol = 1, nrow = 3)
  	print(gbg)
  dev.off()
  
  save.image(here:::here("output_mystery",paste0("MysterySample_Prediction_overall_",modelDir,"_",tlevel,norm0,log0,".rdata")))
}


#################
### checking which files were single- versus paired-end
#################
preds2<-data.frame(preds, predsInfo)
motSingle<-merge(truthAll, preds2, by.y="ID",by.x="sample", sort=F)
