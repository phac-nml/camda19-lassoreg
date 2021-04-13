
###
####################################################
####
#### comparing 16S and shotgun metagenomics data (Boston Pilot)
####
##################################################

require(metagenomeSeq)
require(here)

#outputFolder="outputBoston_Filt100_SampPres1"
#outputFolder="outputBoston"

dir2save="comb"
if(!dir.exists(here:::here(outputFolder,dir2save))){dir.create(here:::here(outputFolder,dir2save))}

tlevels=c("species","genus", "family", "order","class")


setwd(here:::here(outputFolder,dir2save))

### Kraken2 Bracken
load("../brackenBoston.rdata")
wgs<-metaobj
table(fData(wgs)$kingdom)
(cntsBe4_wgs<-apply(fData(wgs), 2, function(y) length(table(y)))[c(tlevels, "phylum")])
### for a 'fair' comparison between shotgun and 16S, removing Archaea, Viruses and Eukaryota
wgs<-wgs[which(fData(wgs)$kingdom=="Bacteria"),]


###metaphlan
load("../brackenBostonMetaphlan.rdata")
wgsM<-metaobj
table(fData(wgsM)$kingdom)
(cntsBe4_wgsM<-apply(fData(wgsM), 2, function(y) length(table(y)))[c(tlevels, "phylum")])
### for a 'fair' comparison between shotgun and 16S, removing Archaea, Viruses and Eukaryota
wgsM<-wgsM[which(fData(wgsM)$kingdom=="Bacteria"),]


### 16S
load("../brackenBoston16S.rdata")
s16<-metaobj
heatmapCols = colorRampPalette(brewer.pal(9, "GnBu"))(50)
cbbPalette <- c( "#E69F00", "#56B4E9", "#8B2323","#66CD00","#9932CC", "#D55E00", "#00008B", "#CDC0B0")



#### returns the list of count sums at least taxa level
extractTaxa2list<-function(x, norm0=T, log0=F, aggfun0 = colSums){
  t0<-c("kingdom", "phylum",  "class" , "order" ,  "family" ,"genus" , "species")
  ret<-lapply(t0, function(y){
    aggTax(x, lvl = y, out = "matrix" ,norm=norm0, log=log0, aggfun = aggfun0)
  })
  names(ret)=t0
  ret
}

### x and y are genus tables to be compared
jFindCommon<-function(x,y){
  xn<-rownames(x)[rownames(x)%in%rownames(y)]
  xn<-setdiff(xn,"")
  cnts=c(nrow(x),nrow(y),length(xn))
  names0=c("x taxa count", "y taxa count", "overlapping count")
  #print(paste(names0, collapse=" "))
  #print(cnts)
  names(cnts)=names0
  list(x[xn,],y[xn,], cnts=cnts)
}




###modified from https://www.r-bloggers.com/example-8-41-scatterplot-with-marginal-histograms/
### need to adjust the matching still
scatterHist = function(x, y, main="", xlab="", ylab="",xlim=NULL,ylim=NULL, cex2use=2.2){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  nbins=10
  xhist = hist(x, breaks=seq(xlim[1],xlim[2],diff(xlim)/nbins), plot=FALSE)
  yhist = hist(y, breaks=seq(ylim[1],ylim[2],diff(ylim)/nbins), plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  
  par(mar=c(4,4,1,1))
  plot(x,y,xlab="",ylab="", xlim=xlim,ylim=ylim, pch=20, cex.axis=cex2use)
  text(xlim[2]*0.03,ylim[2]*0.85, main,cex=cex2use*1.2,adj=c(0,-1.5))
  
  par(mar=c(0,4,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  
  par(mar=c(4,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1.5, outer=TRUE, adj=0, at=0.15,cex=cex2use*0.8)
  mtext(ylab, side=2, line=1.5, outer=TRUE, adj=0, at=0.15,cex=cex2use*0.8)
}



norm0=T
log0=F ### not logged here as doing this later



wgsList<-extractTaxa2list(wgs, norm0=norm0, log0=log0)
s16List<-extractTaxa2list(s16, norm0=norm0, log0=log0)
wgsMList<-extractTaxa2list(wgsM, norm0=norm0, log0=log0)

##############################
### finding common taxa between wgs and 16S: Table 1
####################################
gcomL1<-list()
gcomL3=gcomL2=gcomL1
for (tlevel in c("species","genus", "family", "order","class", "phylum")){
  
  print(tlevel)
  gcomL1[[tlevel]]<-jFindCommon(wgsList[[tlevel]], s16List[[tlevel]])
  gcomL2[[tlevel]]<-jFindCommon(wgsMList[[tlevel]], s16List[[tlevel]])
  gcomL3[[tlevel]]<-jFindCommon(wgsList[[tlevel]],wgsMList[[tlevel]])
}  

gcomL1ol<-sapply(gcomL1, function(y)y[[3]])
gcomL2ol<-sapply(gcomL2, function(y)y[[3]])
gcomL3ol<-sapply(gcomL3, function(y)y[[3]])

olfile="Table1_OverlappingCounts.xlsx"
write.xlsx(gcomL1ol,olfile,sheetName = "SG_KBto16S",col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(gcomL2ol,olfile,sheetName = "SG_MPto16S",col.names = TRUE, row.names = TRUE, append = T)
write.xlsx(gcomL3ol,olfile,sheetName = "SG_KBtoSG_MP",col.names = TRUE, row.names = TRUE, append = T)


##############################
### Figure 1
##############################
corLevels=rep(0,length(tlevels)); names(corLevels)=tlevels
corLevels2=corLevels

pdf(paste0("Fig1u_ScatterPlot_comparisonBySample_16StoShotgun_preNorm",norm0,"_log",log0,"_taxaALL.pdf"), width=8,height=8)
  par(mfrow=c(1,1))
  par0<-par()
  for (tlevel in tlevels){
    
    ##############################
    ### finding common taxa between wgs and 16S
    ####################################
    
    gcom<-jFindCommon(wgsList[[tlevel]], s16List[[tlevel]])
  
    
    techCol=c("purple3","darkgoldenrod3")
    gcomCol=c(rep(techCol[1],ncol(gcom[[1]])), rep(techCol[2],ncol(gcom[[2]])))
    
  
   allin1<-cbind(wgs=log2(as.vector(gcom[[1]])+1),x16s=log2(as.vector(gcom[[2]])+1))
    #print(range(allin1))
    scatterHist(allin1[,1],allin1[,2], main=paste0(tlevel, ": cor=", signif(cor(allin1[,1],allin1[,2]),4)), xlim=range(allin1),ylim=range(allin1), ylab="log2(16S abundance +1)", xlab="log2(SG-KB abundance +1)")
    corLevels[tlevel]=cor(allin1[,1],allin1[,2])
    
    ## conf. interval for Pearson correlation coefficient
    print(paste(tlevel, paste(signif(cor.test(allin1[,1],allin1[,2], conf.level = 0.95)$conf.int, 2), collapse=" ")))
   par(par0)
  }  
dev.off()

print(corLevels)


limUsed=range(allin1)

### SG-MP and 16S Figure 1 lower

pdf(paste0("Fig1l_ScatterPlot_comparisonBySample_16StoShotgun_preNorm",norm0,"_log",log0,"_taxaALL_SG-MP.pdf"), width=8,height=8)
par(mfrow=c(1,1))
par0<-par()
for (tlevel in tlevels){
  
  gcom<-jFindCommon(wgsMList[[tlevel]], s16List[[tlevel]])
  
  techCol=c("purple3","darkgoldenrod3")
  gcomCol=c(rep(techCol[1],ncol(gcom[[1]])), rep(techCol[2],ncol(gcom[[2]])))
  
  allin2<-cbind(wgs=log2(as.vector(gcom[[1]])+1),x16s=log2(as.vector(gcom[[2]])+1))
  #print(range(allin2))
  scatterHist(allin2[,1],allin2[,2], main=paste0(tlevel, ": cor=", signif(cor(allin2[,1],allin2[,2]),4)), xlim=limUsed,ylim=limUsed, ylab="log2(16S abundance +1)", xlab="log2(SG-MP abundance +1)")
  corLevels2[tlevel]=cor(allin2[,1],allin2[,2])
  
  ## conf. interval for Pearson correlation coefficient
  print(paste(tlevel, paste(signif(cor.test(allin2[,1],allin2[,2], conf.level = 0.95)$conf.int, 2), collapse=" ")))
  par(par0)
}  
dev.off()

write.xlsx(data.frame(kb2s16=corLevels,mg216s=corLevels2),olfile,sheetName = "SG_KBto16S_cor",col.names = TRUE, row.names = TRUE, append = T)

##############################
### Figure 2
##############################
tlevel="genus"
gcom<-jFindCommon(wgsList[[tlevel]], s16List[[tlevel]])

techCol=c("purple3","darkgoldenrod3")
gcomCol=c(rep(techCol[1],ncol(gcom[[1]])), rep(techCol[2],ncol(gcom[[2]])))


pdf(paste0("Fig2.Heatmaps_comparison16StoWGS_preNorm",norm0,"_log",log0,"_taxa",tlevel,".pdf"), width=10, height=10)
  
  a=cbind(gcom[[1]],gcom[[2]])
  colnames(a)=paste(colnames(a),rep(c("wgs","16s"), c(23,23)), sep="_")
  require(gplots)
  ##Bray Curtis
  library(vegan)
  baSamp=vegdist(t(a),method="bray")
  heatmap.2(as.matrix(baSamp), trace="none", scale="none", main="samp_Bray-Curtis", keysize=0.7, margin=c(10,10), col=rev(heatmapCols), ColSideColors=gcomCol)
  legend("topright", title = "Tech",legend=c("SG-KB","16S"), fill=techCol, cex=0.8, box.lty=0)
  heatmap.2(as.matrix(baSamp), trace="none", scale="none", main="samp_Bray-Curtis", labRow="", labCol="", keysize=0.7, margin=c(5,3), col=rev(heatmapCols), ColSideColors=gcomCol, RowSideColors=c("lightsalmon1","wheat4")[c(cl,cl)])
  legend("bottom", title = "Tech & Source",legend=c("SG-KB","16S",  levels(cl)), fill=c(techCol,"lightsalmon1","wheat4"), cex=0.8, box.lty=0, horiz=T)
  heatmap.2(as.matrix(baSamp), trace="none", scale="none", main="samp_Bray-Curtis", labRow="", labCol="", keysize=0.7, margin=c(5,3), col=rev(heatmapCols), ColSideColors=gcomCol, RowSideColors=c("springgreen3", "orangered2", "dodgerblue1", "darkgrey")[c(cl2,cl2)])
  legend("bottom", title = "Tech & Source",legend=c("SG-KB","16S",  levels(cl),  levels(cl2)), fill=c(techCol,"lightsalmon1","wheat4","springgreen3", "orangered2", "dodgerblue1", "darkgrey"), cex=0.8, box.lty=0, ncol=4)
  
  clsource=gcomCol #rep(1:2,c(23,23))
  pch2use=c(1: length(levels(cl2)))
  
  
  require(ape)
  dimn=5
  PCOA <- pcoa(baSamp)
  if(sum(PCOA$values$Relative_eig<0)>0){
    pcoac<-pcoa(baSamp, correction = "cailliez")
    PCOAaxes <- pcoac$vectors[,c(1:dimn)]
    eignPERC<- pcoac$values$Rel_corr_eig[c(1:dimn)]
    pairs(PCOAaxes,   main=paste("PCo_Bray-Curtis_16S+Shotgun: Surface - cailliez correction"), 
          col=clsource, pch=pch2use[cl2], cex=1.5, cex.labels = 2, cex.axis=1.5, upper.panel=NULL,
          labels=paste("Dim",1:dimn,"\n",round(eignPERC,3)*100))  ###B
    pairs(PCOAaxes,   main=paste("PCo_Bray-Curtis_16S+Shotgun: Source"), col=clsource, pch=pch2use[cl])
    plot(PCOAaxes[,1], PCOAaxes[,2], xlab="PCo 1", ylab="PCo 2",  main=paste("PCo_Bray-Curtis_16S+Shotgun"), col=clsource, pch=pch2use[cl2])
    legend("bottom",c("Shotgun","16S"),col=techCol,pch=19, ncol=6) #
    legend("bottom",levels(cl2),pch=pch2use) #,pch=pch2use
    
  }else{
    PCOAaxes <- PCOA$vectors[,c(1:dimn)]
    eignPERC<- PCOA$values$Relative_eig[c(1:dimn)]
    pairs(PCOAaxes,   main=paste("PCo_Bray-Curtis_16S+Shotgun: Surface"), 
          col=clsource, pch=pch2use[cl2], cex=1.5, cex.labels = 2, cex.axis=1.5, upper.panel=NULL,
          labels=paste("Dim",1:dimn,"\n",round(eignPERC,3)*100))  ###B
    pairs(PCOAaxes,   main=paste("PCo_Bray-Curtis_16S+Shotgun: Source"), col=clsource, pch=pch2use[cl])
    plot(PCOAaxes[,1], PCOAaxes[,2], xlab="PCo 1", ylab="PCo 2",  main=paste("PCo_Bray-Curtis_16S+Shotgun"), col=clsource, pch=pch2use[cl2])
    legend("bottom",c("Shotgun","16S"),col=techCol,pch=19, ncol=6) #
    legend("bottom",levels(cl2),pch=pch2use) #,pch=pch2use
  }
  
  ### old
  fit <- cmdscale(baSamp,eig=TRUE, k=5) # k is the number of dim
  # plot solution 
  x <- fit$points[,1]
  y <- fit$points[,2]
  pairs(fit$points,   main=paste("PCo_Bray-Curtis_16S+Shotgun: Surface"), col=clsource, pch=pch2use[cl2], cex=1.5, cex.labels = 2, cex.axis=1.5, upper.panel=NULL)  ###B
  #legend("bottom",levels(cl2),pch=pch2use, cex=0.5) #didn't work grab from indiv plot
  pairs(fit$points,   main=paste("PCo_Bray-Curtis_16S+Shotgun: Source"), col=clsource, pch=pch2use[cl])
  plot(x, y, xlab="PCo 1", ylab="PCo 2",  main=paste("PCo_Bray-Curtis_16S+Shotgun"), col=clsource, pch=pch2use[cl2])
  legend("bottom",c("Shotgun","16S"),col=techCol,pch=19, ncol=6) #
  legend("bottom",levels(cl2),pch=pch2use) #,pch=pch2use

dev.off()
#rm(list = ls(all.names = TRUE))
