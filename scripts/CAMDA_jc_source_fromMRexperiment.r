### post initial MRexperiment object
require(amap)
require(ape)
require(metagenomeSeq)
require(ggplot2)
require(ggpubr)
require(gplots)
require(RColorBrewer)
require(vegan)

#presentThresh=8
#depthThresh=1
#perct=""

#keyFeat="City2"
#feat2="contin"
#toexport=FALSE

print("depthThresh is no longer used.")

metaobj0<-metaobj
### the input data in metaobj: OTU x samplID 




### colors to be refined


darkcols <- c(brewer.pal(8, "Accent"),rev(brewer.pal(8, "Dark2")[-8]), brewer.pal(8,"Paired"))
#darkcols <- c(brewer.pal(8, "Accent")[1:3],"#ffee99",brewer.pal(8, "Accent")[5:8],rev(brewer.pal(8, "Dark2")[-8]), brewer.pal(8,"Paired"))
  
cbbPalette <- c( "#E69F00", "#56B4E9", "#8B2323","#66CD00","#9932CC", "#D55E00", "#00008B", "#CDC0B0")
cols=darkcols
cols2=brewer.pal(10,"Set3")[-2]

#heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
heatmapCols = colorRampPalette(brewer.pal(9, "GnBu"))(30)

cbbPalette=rainbow(30)


### added factor
cl=factor(pData(metaobj)[,keyFeat])
clcol=cols[as.integer(cl)]

cl2=factor(pData(metaobj)[,feat2])
cl2col=cols2[cl2]

if(!exists("clplus")){
	clplus=levels(cl)
}

### histograms for distribution
pdf(here:::here(outputFolder,paste0("Hist_Sums_raw_counts_",fileNameAdd,"_b4filter.pdf")))
	hist(log2(colSums(MRcounts(metaobj, norm = F, log =F))),xlab="log2(Sum of counts)", main="By Samples", breaks=50)
	hist(log2(rowSums(MRcounts(metaobj, norm = F, log =F))),xlab="log2(Sum of counts)", main="By OTUs", breaks=50)
	hist(rowSums(MRcounts(metaobj, norm = F, log =F)!=0),xlab="Sum of samples w nonzero counts", main="By OTUs", breaks=50)
dev.off()


###### boxplot
pdf(here:::here(outputFolder,paste0("Boxplot_",fileNameAdd,"_b4filter.pdf")), width=20)
   boxplot(log2(1+MRcounts(metaobj, norm = F, log =F)), col=clcol, outcol=clcol, ylab="log2(1+Abundance)")
dev.off()




    # old filtering with one read count
    ###############
    ###filter OTU presence or minimum depth 
    ###############
    #present:Features with at least 'present' positive samples.
    #depth:Samples with at least this much depth of coverage
    
    #(dimb4=dim(MRcounts(metaobj0))) ##b4
    #metaobj=filterData(metaobj0, present = presentThresh, depth = depthThresh)
    #(dimaf=dim(MRcounts(metaobj))) ##after

###############
### revision edit on filtering
###filter OTU presence with a chosen read count
###############

### modified filter function with 
### added fdepth for feature count cutoff for presence. A taxa/OTU has to have at least this number of reads to be deemed present within a sample.
filterDataMod <-function (obj, present = 1, fdepth= 1, depth = 1) {
  mat = returnAppropriateObj(obj, norm = FALSE, log = FALSE) > (fdepth-1) 
  cols = which(colSums(MRcounts(obj)) >= depth)
  rows = which(rowSums(mat[, cols]) >= present)
  return(obj[rows, cols])
}

(dimb4=dim(MRcounts(metaobj0))) ##b4
metaobj=filterDataMod(metaobj0, present = presentThresh, fdepth= fdepthThresh,  depth = depthThresh)
(dimaf=dim(MRcounts(metaobj))) 

###############
###Calculating normalization factors
###############
if(perct==""){
	perct = cumNormStat(metaobj) 
	perct
}
metaobj = cumNorm(metaobj, p = signif(perct,4))



###moved below norm on 20190510
### histograms for distribution after filtering
pdf(here:::here(outputFolder,paste0("Hist_Sums_raw_counts_",fileNameAdd,"_affilter.pdf")))
	hist(log2(colSums(MRcounts(metaobj, norm = F, log =F))), xlab="log2(Sum of counts)", main="By Samples", breaks=50)
	hist(log2(rowSums(MRcounts(metaobj, norm = F, log =F))), xlab="log2(Sum of counts)", main="By OTUs", breaks=50)
	hist(rowSums(MRcounts(metaobj, norm = F, log =F)!=0), xlab="Sum of samples w nonzero counts", main="By OTUs", breaks=50)
dev.off()

###### boxplot after filtering
pdf(here:::here(outputFolder,paste0("Boxplot_",fileNameAdd,"_affilter.pdf")), width=20)
 boxplot(log2(1+MRcounts(metaobj, norm = F, log =F)), col=clcol, outcol=clcol, ylab="log2(1+Abundance)")
 dev.off()


pdf(here:::here(outputFolder,paste0("Heatmap_cor_postFiltering_samp_",fileNameAdd,".pdf")), width=20, height=20)
	corv<-cor(MRcounts(metaobj, norm = F, log =F),use="pairwise.complete.obs", method="spearman")
	sum(is.na(corv)) 
	length(corv)
	corv[is.na(corv)]=0
	heatmap.2(corv, trace="none", scale="none", margin=c(10,10), keysize=0.55, main="Spearman", ColSideColors=clcol, col=heatmapCols)
dev.off()

 


###### boxplot after filtering & normalization
pdf(here:::here(outputFolder,paste0("Boxplot_",fileNameAdd,"_affilterNorm_",signif(perct,4),".pdf")), width=20)
 boxplot(log2(1+MRcounts(metaobj, norm = T, log =F)), col=clcol, outcol=clcol, ylab="log2(1+Abundance)")
 dev.off()


###############
### export normalized count matrices and stats
###############
if(toexport){
	nmat = cbind(as.matrix(fData(metaobj)),MRcounts(metaobj, norm = TRUE, log = FALSE))
	exportMat(nmat, file = paste0("normedcnt_",perct,".tsv"))
	rmat =  cbind(as.matrix(fData(metaobj)),MRcounts(metaobj, norm = FALSE, log = FALSE))
	exportMat(rmat, file = paste0("rawcnt_",perct,".tsv"))
	exportStats(metaobj, p=perct,file = paste0("normStats_",perct,".tsv"))
	#head(read.csv(file = paste0(resultDir,"/normStats.tsv"), sep = "\t"))
}

pdf(here:::here(outputFolder,paste0("Hist_Sums_Normed",signif(perct,4),"_counts_afterfilter_",fileNameAdd,"_present",presentThresh,".pdf")))
  hist(log2(colSums(MRcounts(metaobj, norm = T))),xlab="log2(Sum of normed counts)", main="By Samples", breaks=50)
  hist(log2(rowSums(MRcounts(metaobj, norm = T))),xlab="log2(Sum of normed counts)",main="By OTUs", breaks=50)
dev.off()

### printing out the outlying samples in terms of total raw counts after filter
print(head(sort(log2(colSums(MRcounts(metaobj))))))
print(tail(sort(log2(colSums(MRcounts(metaobj))))))

pdf(here:::here(outputFolder,paste0("Boxplot_log2_counts_afterfilter_",fileNameAdd,"_present",presentThresh,".pdf")), width=12)
boxplot( MRcounts(metaobj, norm = F, log = TRUE), col=clcol, outcol=clcol, main="Before Norm")
boxplot( MRcounts(metaobj, norm = TRUE, log = TRUE), col=clcol, outcol=clcol, main="After Norm")
dev.off()

 
### the norm doesnt change Spearman cor really 
####################
###### PCoA plot - spearman
####################
#require(amap)
system.time(corvd<-Dist(t(MRcounts(metaobj, norm = F, log = F)), method="spearman"))
system.time(fitcor <- cmdscale(corvd,eig=TRUE, k=3)) # k is the number of dim


head(fitcor$eig)

save2pdf=F
# plot solution 
x <- fitcor$points[,1]
y <- fitcor$points[,2]

pdf(here:::here(outputFolder,paste0("PCO_SpearmanCor_",fileNameAdd,".pdf")), height=10, width=10)
	pch2use=c(1: length(levels(cl2)))
	plot(x, y, xlab="PCo 1", ylab="PCo 2",  main=paste("PCo", fileNameAdd, "SpearmanCor"), col=clcol, pch=pch2use[cl2])
	legend("bottom",legend=levels(cl),col=cols, pch=pch2use[clplus], ncol=6) #
	legend("left",legend=levels(cl2),pch=pch2use) #,pch=pch2use
	
	plot(x, y, xlab="PCo 1", ylab="PCo 2",  main=paste("PCo", fileNameAdd, "SpearmanCor"), type="n")
	text(x, y, labels = cl, cex=.7, col=clcol ) #col=clcol
	legend("topleft",legend=levels(cl),col=cols,pch=pch2use[clplus])
	
	pairs(fitcor$points,   main=paste("PCo", fileNameAdd, "SpearmanCor"), col=clcol, pch=pch2use[cl2])
	
	plot(x, y, xlab="PCo 1", ylab="PCo 2",  main=paste("PCo", fileNameAdd, "SpearmanCor"), col=clcol, pch=20)
	legend("bottom",legend=levels(cl),col=cols,pch=20, ncol=6)
dev.off()

topEig<-min(50,length(fitcor$eig))
pdf(here:::here(outputFolder,paste0("PCO_SpearmanCor_VarianceExplained_top",topEig,"_",fileNameAdd,"_2.pdf")), height=10, width=10)
barplot(cumsum(fitcor$eig^2/sum(fitcor$eig^2)*100)[1:topEig], main=paste("CumSum of Variance explained, Spearman, top", topEig))
dev.off()




###################### 
#### Bray Curtis
###################### 
##### Bray Curtis before norm
#require(vegan)
df <- vegdist(t(MRcounts(metaobj, norm = F, log = F)), method="bray") # Bray-Curtis
fitf <- cmdscale(df,eig=TRUE, k=3) # k is the number of dim

#### after norm
d <- vegdist(t(MRcounts(metaobj, norm = T, log = F)), method="bray") # Bray-Curtis
fit <- cmdscale(d,eig=TRUE, k=3) # k is the number of dim


save2pdf=F
# plot solution
x <- fitf$points[,1]
y <- fitf$points[,2]


pdf(here:::here(outputFolder,paste0("PCO_Bray-Curtis_beforeNorm_",fileNameAdd,"_2.pdf")), height=10, width=10)
	pch2use=c(1: length(levels(cl2)))
	plot(x, y, xlab="PCo 1", ylab="PCo 2", main=paste("PCo", fileNameAdd), col=clcol, pch=pch2use[cl2])
	legend("bottom",legend = levels(cl),col=cols, pch=pch2use[clplus], ncol=6) #
	legend("left",legend = levels(cl2),pch=pch2use) #,pch=pch2use
	plot(x, y, xlab="PCo 1", ylab="PCo 2", main=paste("PCo", fileNameAdd), type="n")
	text(x, y, labels = cl, cex=.7, col=clcol ) #col=clcol
	legend("topleft",legend = levels(cl),col=cols,pch=pch2use[clplus])
	pairs(fitf$points, main=paste("PCo", fileNameAdd), col=clcol, pch=pch2use[cl2])
	plot(x, y, xlab="PCo 1", ylab="PCo 2", main=paste("PCo", fileNameAdd), col=clcol, pch=20)
	legend("bottom",legend = levels(cl),col=cols,pch=20, ncol=6)

	topEig<-min(50,length(fitf$eig))
	barplot(cumsum(fitf$eig^2/sum(fitf$eig^2)*100)[1:topEig], main=paste("CumSum of Variance explained, top", topEig))
dev.off()

###################### 
#### Bray Curtis
###################### 
#library(gplots)
titleAdd=""
sampBCTT=as.matrix(df)
pdf(here:::here(outputFolder,paste0("Heatmap_sampleDistance_Bray-Curtis_postFilteringBeforeNorm_OTU_",fileNameAdd,"_withContinent.pdf")), width=20, height=20)
  ###Bray–Curtis dissimilarity
  heatmap.2(sampBCTT, trace="none", scale="none", ColSideColors=clcol, RowSideColors=cl2col, main="samp_Bray-Curtis", col=rev(heatmapCols), keysize=0.7, margin=c(10,10))
  legend("left",  legend = levels(cl), col = cols[1:length(levels(cl))], lty= 1, lwd = 5, cex=1, box.col = NA)
  legend("bottomleft",  legend = levels(cl2), col = cols2[1:length(levels(cl2))], lty= 1, lwd = 5, cex=1, box.col = NA)

dev.off()


#library(gplots)
titleAdd=""
sampBCTT=as.matrix(d)


### figure 3A
pdf(here:::here(outputFolder,paste0("Heatmap_sampleDistance_Bray-Curtis_postFilteringNorm",signif(perct,4),"_logCnt_OTU_",fileNameAdd,"_withContinent.pdf")), width=20, height=20)
  ###Bray–Curtis dissimilarity
  heatmap.2(sampBCTT, trace="none", scale="none", ColSideColors=clcol, RowSideColors=cl2col, main="samp_Bray-Curtis", col=rev(heatmapCols), keysize=0.7, margin=c(10,10))
  legend("top",  legend = levels(cl), col = cols[1:length(levels(cl))], lty= 1, lwd = 5, cex=1, box.col = NA, ncol=8, inset=c(0,-.005))
  legend("topleft",  legend = levels(cl2), col = cols2[1:length(levels(cl2))], lty= 1, lwd = 5, cex=1, box.col = NA, inset=c(-.007,0.1))

dev.off()


####################
###### PCoA plot - Bray Curtis
####################

#fit # view results

save2pdf=F
# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]


pdf(here:::here(outputFolder,paste0("PCO_Bray-Curtis_normed",signif(perct,4),"_",fileNameAdd,"_2.pdf")), height=10, width=10)
	pch2use=c(1: length(levels(cl2)))
	plot(x, y, xlab="PCo 1", ylab="PCo 2",  main=paste("PCo", fileNameAdd), col=clcol, pch=pch2use[cl2])
	legend("bottom",legend = levels(cl),col=cols, pch=pch2use[clplus], ncol=6) #
	legend("left",legend = levels(cl2),pch=pch2use) #,pch=pch2use
	
	plot(x, y, xlab="PCo 1", ylab="PCo 2",  main=paste("PCo", fileNameAdd), type="n")
	text(x, y, labels = cl, cex=.7, col=clcol ) #col=clcol
	legend("topleft",legend = levels(cl),col=cols,pch=pch2use[clplus])
	
	pairs(fit$points,   main=paste("PCo", fileNameAdd), col=clcol, pch=pch2use[cl2], upper.panel=NULL, cex.labels = 2, cex.axis=1.5)
	
	plot(x, y, xlab="PCo 1", ylab="PCo 2",  main=paste("PCo", fileNameAdd), col=clcol, pch=20)
	legend("bottom",legend = levels(cl),col=cols,pch=20, ncol=6)

dev.off()

topEig<-min(50,length(fit$eig))
pdf(here:::here(outputFolder,paste0("PCO_Bray-Curtis_normed",signif(perct,4),"_VarianceExplained_top",topEig,"_",fileNameAdd,"_2.pdf")), height=10, width=10)
  barplot(cumsum(fit$eig^2/sum(fit$eig^2)*100)[1:topEig], main=paste("CumSum of Variance explained, top", topEig))
dev.off()

### figure 3B
#require(ape)
pch2use=c(1: length(levels(cl2)))
dimn=5
PCOA <- pcoa(d)

pdf(here:::here(outputFolder,paste0("PCO_Bray-Curtis_normed",signif(perct,4),"_",fileNameAdd,"_APE.pdf")), height=10, width=10)
if(sum(PCOA$values$Relative_eig<0)>0){
  pcoac<-pcoa(d, correction = "cailliez")
  PCOAaxes <- pcoac$vectors[,c(1:dimn)]
  eignPERC<- pcoac$values$Rel_corr_eig[c(1:dimn)]
  pairs(PCOAaxes, main=paste("PCo", fileNameAdd,"- cailliez correction"), col=clcol, pch=pch2use[cl2], 
        cex=1.1, cex.labels = 2, cex.axis=1.5, upper.panel=NULL,
        labels=paste("Dim",1:dimn,"\n",round(eignPERC,3)*100,"%"))  ###B
  legend("topright",legend = levels(cl),col=cols, pch=pch2use[clplus], ncol=6, cex=1, inset = c(0.1, 0.1)) #
  legend("right",legend = levels(cl2),pch=pch2use, cex=1, inset = c(0.1, -0.1)) #,pch=pch2use
  
}else{
  PCOAaxes <- PCOA$vectors[,c(1:dimn)]
  eignPERC<- PCOA$values$Relative_eig[c(1:dimn)]
  pairs(PCOAaxes, main=paste("PCo", fileNameAdd), col=clcol, pch=pch2use[cl2], 
      cex=1.1, cex.labels = 2, cex.axis=1.5, upper.panel=NULL,
        labels=paste("Dim",1:dimn,"\n",round(eignPERC,3)*100,"%"))  ###B
  legend("topright",legend = levels(cl),col=cols, pch=pch2use[clplus], ncol=6, cex=1, inset = c(0.1, 0.1)) #
  legend("right",legend = levels(cl2),pch=pch2use, cex=1, inset = c(0.1, -0.1)) #,pch=pch2use
}
dev.off()


### density for fig 3b
lp<-lapply(1:dimn, function(i){
  contDens<-data.frame(x=PCOAaxes[,i], continent=cl2)
  pden<-ggplot(contDens, aes(x=x, col=continent))+geom_density()+scale_color_manual(values= cols2[1:length(levels(cl2))])+theme_classic() + xlab(paste0("Dim",i))
  pden
})

densplots<-ggarrange(plotlist=lp, ncol=3,nrow=2)#, labels = paste0("Dim",1:dimn)
ggexport(densplots, filename=here:::here(outputFolder,paste0("PCO_Bray-Curtis_normed",signif(perct,4),"_",fileNameAdd,"_APE_Density.pdf")), width=14, height=7)


save.image(here:::here(outputFolder,paste0(fileNameAdd,".rdata")))





