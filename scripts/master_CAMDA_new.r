#############
## This is the master script on the CAMDA project to call on scripts for specific components
#############
#cd /home/CSCScience.ca/chichen/CAMDA/CAMDA_Rproject/camda19-lassoreg
#source activate R3.6.0
#srun -c 4 --mem 128G --partition NMLResearch R --vanilla
#srun -c 2 --mem 64G --partition NMLResearch R --vanilla

library(here)
here()


###############################################################
###
### main metasub data
###
###############################################################

###****************************************###
###
### 01 data preprocessing
###
###****************************************###

############### paired only                    ###############

### taxa level
#tlevel="species"; #tlevel="genus" #tlevel="family"
tlevels=c("species","genus","family")
for(tlevel in tlevels){
  ### filtering threshold for metagenomeSeq
  fdepthThresh=100
  presentThresh=8
  depthThresh=1
  ### folder to save results
  fileNameAdd=paste0("bracken_",tlevel)
  outputFolder="output"
  ### Run script to process data
  source(here:::here("scripts","01_CAMDA_jc_paired.R"))
  rm(list=setdiff(ls(), c("tlevel","tlevels")))
}

############### single and paired (include SAC) ###############
### taxa level
#tlevel="species"; #tlevel="genus" #tlevel="family"
tlevels=c("species","genus","family")
for(tlevel in tlevels){
    ### filtering threshold for metagenomeSeq
  fdepthThresh=100
  presentThresh=8
  depthThresh=1
  ### folder to save results
  fileNameAdd=paste0("bracken_",tlevel)
  outputFolder="output_pairedsingle"
  ### Run script to process data
  source(here:::here("scripts","01_CAMDA_jc_paired_single.R"))
  rm(list=setdiff(ls(), c("tlevel","tlevels")))
}


###****************************************###
###
### 02 Vis & Run ML to build model
###
###****************************************###

############### paired only                    ###############

### glmnet lasso-regularized classification
### taxa level
#tlevel="species"; #tlevel="genus" #tlevel="family"
tlevels=c("species","genus","family")
for(tlevel in tlevels){
  ### data source
  outputFolder="output"
  #weighted
  seed2use=111
  source(here:::here("scripts","02_CAMDA_jc_glmnet_class_weighted.r"))
  rm(list=setdiff(ls(), c("outputFolder","tlevel","tlevels")))
}

### glmnet lasso-regularized regression
### taxa level
#tlevel="species"; #tlevel="genus" #tlevel="family"
tlevels=c("species","genus","family")
for(tlevel in tlevels){
  ### data source
  outputFolder="output"
  seed2use=111
  source(here:::here("scripts","02_CAMDA_jc_glmnet_mse.r"))
  rm(list=setdiff(ls(), c("outputFolder","tlevel","tlevels")))
}

### RF classification w balanced sampling 
### taxa level
#tlevel="species"; #tlevel="genus" #tlevel="family"
tlevels=c("species","genus","family")
for(tlevel in tlevels){
  ### data source
  outputFolder="output"
  source(here:::here("scripts","02_CAMDA_jc_RFClass_balanced.r"))
  rm(list=setdiff(ls(), c("outputFolder","tlevel","tlevels")))
}


############### single and paired (include SAC) ###############
### glmnet lasso-regularized classification
### taxa level
tlevels=c("species","genus","family")
for(tlevel in tlevels){
  ### data source
  outputFolder="output_pairedsingle"
   #weighted
  seed2use=111
  source(here:::here("scripts","02_CAMDA_jc_glmnet_class_weighted.r"))
  rm(list=setdiff(ls(), c("outputFolder","tlevel","tlevels")))
}


### glmnet lasso-regularized regression
### taxa level
tlevels=c("species","genus","family")
for(tlevel in tlevels){
  ### data source
  outputFolder="output_pairedsingle"
  seed2use=111
  source(here:::here("scripts","02_CAMDA_jc_glmnet_mse.r"))
  rm(list=setdiff(ls(), c("outputFolder","tlevel","tlevels")))
}

    
### RF classification balanced sampling 
### taxa level
tlevels=c("species","genus","family")
for(tlevel in tlevels){
  ### data source
  outputFolder="output_pairedsingle"
  source(here:::here("scripts","02_CAMDA_jc_RFClass_balanced.r"))
  rm(list=setdiff(ls(), c("outputFolder","tlevel","tlevels")))
}

########
#### MSE comparisons at modelling level: Table 2
########
outputFolder="output"
for(ll in c("species","genus","family")){
  load(here:::here(outputFolder,paste0("glmnet_MSE/Multivariate_",ll,"_lassolog_TRUE_gps_MSE.rdata")))
  which1se<-which(round(cvob2mse$glmnet.fit$lambda,6)==round(cvob2mse$lambda.1se,6))
  df1se<-cvob2mse$glmnet.fit$df[which1se]
  cvm1se<-cvob2mse$cvm[which1se]
  print(ll)
  print(paste("lambda 1se: ",cvob2mse$lambda.1se, cvob2mse$lambda[which1se]))
  print(paste("$cvm MSE value",cvm1se))
  print(paste("df",df1se))
  
  load(here:::here(outputFolder,paste0("glmnet_MSE/Test_Validation_l1o_CV10_",ll,"_lassolog_TRUE.rdata")))
  print(paste("Leave1cityout MSE",round(mean(l1oDF$llSE),2), ";",
        "10-fold nested MSE",round(mean(CV10DF$llSE),2) ))
}
rm(list=setdiff(ls(), c("outputFolder","tlevel","tlevels")))

###****************************************###
###
### 03 Process mystery samples
###
###****************************************###
### taxa level
#tlevel="species"; #tlevel="genus" #tlevel="family"
tlevels=c("species","genus","family")
for(tlevel in tlevels){
  source(here:::here("scripts","03_CAMDA_jc_mystery.r"))
  rm(list=setdiff(ls(), c("tlevel","tlevels")))
}

###****************************************###
###
### 04 Predition on  mystery samples
###
###****************************************###

### taxa level
tlevel="species" #tlevel="genus"

for(dataDir in c("output", "output_pairedsingle")){
  modelDir="glmnet_MSE";model2use=paste0(modelDir,"/Multivariate_",tlevel,"_lassolog_TRUE_gps_MSE.rdata")
  source(here:::here("scripts","04_CAMDA_jc_mysteryPred.r"))
  rm(list=setdiff(ls(), c("dataDir","tlevel")))
 
  modelDir="glmnet_Classification_weighted";model2use=paste0(modelDir,"/classification_",tlevel,"_lassolog_TRUE_gpsClassification_weighted.rdata")
  source(here:::here("scripts","04_CAMDA_jc_mysteryPred.r"))
  rm(list=setdiff(ls(), c("dataDir","tlevel")))
 
  modelDir="RFClassification_balanced";model2use=paste0(modelDir,"/RFclassification_",tlevel,"log_TRUE.rdata")
  source(here:::here("scripts","04_CAMDA_jc_mysteryPred.r"))
  rm(list=setdiff(ls(), c("dataDir","tlevel")))
}

###****************************************###
###
### 05  Merging PairedSingle and Paired Model Predition on  mystery samples/approach
###
###****************************************###

tlevel="species" ;#tlevel="genus"


models=c("glmnet_MSE", "glmnet_Classification_weighted", "RFClassification_balanced")

for( modelDir in c(models, paste0(models,"_pairedsingle"))){
  source(here:::here("scripts","05_CAMDA_jc_mysteryPred_merge.r"))
  rm(list=setdiff(ls(), "tlevel"))
}


###****************************************###
###
### 06  Merging Preditions from different approaches on  mystery samples
###
###****************************************###
require(naivebayes)
require(ROCR)

tlevel="species" ;#tlevel="genus"
outputFolder="output"
gpsw = "Classification_weighted"
dir2save = paste0("glmnet_",gpsw)

source(here:::here("scripts","06_CAMDA_jc_classification_diversity.r"))
rm(list=setdiff(ls(), c("tlevel", "outputFolder")))

###****************************************###
###
### 07  Report MSEs
###
###****************************************###
tlevel ="species" 
outputFolder="output_mystery"
models <- c("glmnet_MSE",  "glmnet_Classification_weighted", "RFClassification_balanced")
source(here:::here("scripts","07_CAMDA_jc_Final_MSE_comparison.r"))

quit(save = "no")

