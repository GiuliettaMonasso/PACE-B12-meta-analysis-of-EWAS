############################################################## LOAD LIBRARIES #######################################################

library(foreign)
library(data.table)# to process results
library(MASS) # rlm function for robust linear regression
library(sandwich) #Huberís estimation of the standard error
library(lmtest) # to use coeftest
library(parallel) # to use multicore approach - part of base R
library(plyr) # to process results
library(matrixStats) #for rowSums, rowIQRs etc.
library(tableone)
library(limma)
library(stats)



############################################################## STUDY PARAMETERS #######################################################

## Please only change cohort and analysis.date below ##
cohort <- "GENR" #change to name of your cohort/study
analysis.date <- "20201030" # change to the date at which you perform the analyses


## You do not have to chagne anything in this part of the script, even if you are not able to run one or multiple models ##
maternal <- "maternal" # refers to mothers in who vitamin B12 concentrations were measured in pregnancy, if available in your study
newborn <- "newborn" # refers to newborns in who vitamin B12 concentrations were measured, if available in your study
childhood.model <- "childhood.model"



############################################################## MODELS ###################################################################

## If you are not able to run one or multiple models, please delete the R code for this/these model(s) further below ##

############################################################## SET MODEL PARAMETERS #######################################################

## Please make sure these names are identical to the names in your .dat phenotype file! ###
## Here, you only have to change anything in the script if you want to add selection factors (and include these in your phenotype file as well)


sample.id <-c("sample.id") ### Case identifier
traits <-c("mat.b12", "newborn.b12") ### The exposure variables of interest, remove the exposure that is not available in your cohort
cell.names <- c("bcell", "mono", "cd4t", "cd8t", "gran", "nk") ### Houseman reference set for cell type correction
batch <-c("batch") ### Adjust for batch effects by including the most important covariate(s), such as plate, in the models
sampling <-c("gest.age.sampling") ## only for maternal model, remove if your exposure is cord blood b12 !!
main.covariates <- c("mat.age", "mat.ses", "mat.bmi", "mat.smoking", "parity", "sex", "child.age") ### New covariate: child age at blood sampling.


############################################################## PHENOTYPE FILE #######################################################

## Load phenotype data- you can adjust the original phenotype file you have used to run the analyses. Please adjust as follows:
## -	Estimated cell types: Use the Houseman reference set for cell type correction (including cd8t, cd4t, nk, bcell, mono, gran), instead of the Salas cell types (please remove these from phenotype file).
## -	Please add child age at blood sampling (child.age, continuous in years) to the models.

dataB12 <- read.table("mat.edu.dichotomous.childhood.b12.20200915.dat", header = T, stringsAsFactors=FALSE) ### between quotes, specify the name of your phenotype file
names(dataB12)

## dataB12 is a .dat file (tab delimited) of 24 columns containing information on: sample IDs (1st column), 4 exposures (vitamin b12 concentrations of mothers and newborns), cell counts (Houseman reference, 6 columns), gestational age at sampling, maternal covariates, parity, child sex, child age, folate concentrations, homocysteine concentration (2x2 columns).
## There should be no extra columns, except if you have added selection factors. If this is the case, for questions about how to change the R script, plese contact g.monasso@erasmusmc.nl

## This part is to make sure that categorical covariates are treated as factors and continuous covariates as numeric.
dataB12$mat.b12 = as.numeric(dataB12$mat.b12)
dataB12$newborn.b12 = as.numeric(dataB12$newborn.b12)
dataB12$batch = as.factor(dataB12$batch)
dataB12$bcell = as.numeric(dataB12$bcell)
dataB12$mono = as.numeric(dataB12$mono)
dataB12$cd4t = as.numeric(dataB12$cd4t)
dataB12$cd8t = as.numeric(dataB12$cd8t)
dataB12$gran = as.numeric(dataB12$gran)
dataB12$nk = as.numeric(dataB12$nk)
dataB12$gest.age.sampling = as.numeric(dataB12$gest.age.sampling)
dataB12$mat.age = as.numeric(dataB12$mat.age)
dataB12$mat.ses = as.factor(dataB12$mat.ses)
dataB12$mat.bmi = as.numeric(dataB12$mat.bmi)
dataB12$mat.smoking = as.factor(dataB12$mat.smoking)
dataB12$parity = as.factor(dataB12$parity)
dataB12$sex = as.factor(dataB12$sex)
dataB12$child.age = as.numeric(dataB12$child.age)


## Summarizes your data before excluding any NA's
summary(dataB12)


## Check if all variables are present in phenotype file
for(i in 1:length(c(sample.id, traits, batch, cell.names, sampling, main.covariates))) {
  print(ifelse(c(sample.id, traits, batch, cell.names, sampling, main.covariates)[i] %in% colnames(dataB12)==FALSE,
               paste("CAUTION: the variable called",c(sample.id, traits, batch, cell.names, sampling, main.covariates)[i],"is missing from pheno file"),
               paste("variable called",c(sample.id, traits, batch, cell.names, sampling, main.covariates)[i],"is present in pheno file")))
}



############################################################## METHYLATION FILE #######################################################
load("/home/579005/GENR3/Methylation/GENR_450KMETH_Norm_Release3/GENR_450KMETH_Release3_Betas_3IQR_5y_20190813.RData")
# load("...") ## load methylation file, eg. load("3IQR-trimmed-methylationfile.Rdata") 
## This is a file with the methylation betas as you have them, with CpGs as rows and participants as columns. 
## We expect this file to be 3IQR-trimmed. If you have any questions about running the 3IQR trim, please contact Giulietta: g.monasso@erasmusmc.nl


## Rename methylation file as betaFINAL and make sure it is a matrix
betaFINAL<-x 
betaFINAL <- as.matrix(betaFINAL)


## Check N cases with methylation data available (before excluding NA's) and order the files according to ID
index<-which(colnames(betaFINAL) %in% dataB12$sample.id)
length(index)
beta<-betaFINAL[,index]
beta<-beta[,order(colnames(beta))]
dataB12<-dataB12[order(dataB12$sample.id),] 
ncol(beta)
nrow(beta)
betaFINAL<-beta


################################################ ANALYSIS FOR SELECTION OF CPGS - HITS FROM MAIN MODEL #################################

# Extract only CpGs of interest = 34 hits main model
df<-read.table("B12_hits.txt", header = T, strings=F) 
dim(df)


#Transpose Methylation data
tx <- t(betaFINAL)


### Select these 34 CpGs from the methylation file
select<-tx[,colnames(tx) %in% df$MarkerName]
dim(select) ### subset of CpGs available in your cohort. 

### transpose again
betaFINAL<-t(select)



################### FROM HERE DELETE THOSE PARTS OF THE SCRIPT THAT ARE MEANT TO RUN MODELS YOU CANNOT RUN ###################



############################################################## MATERNAL MAIN MODEL - CHILDHOOD/ADOLESCENCE FOLLOW-UP ANALYSIS #######################################################

### IF YOUR COHORT DOES NOT HAVE MATERNAL BLOOD B12 CONCENTRATIONS IN PREGNANCY AVAILABLE, YOU CAN DELETE THIS PART OF THE SCRIPT 


## Maternal main model: remove cases with NA for any of mat.b12 and covariates
intersecting.samples <- intersect(dataB12$sample.id,colnames(betaFINAL))
dataB12_main_maternal <- na.omit(dataB12[which(dataB12$sample.id %in% intersecting.samples),unique(c(sample.id, "mat.b12", batch, cell.names, sampling, main.covariates))])
k_main_maternal = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_main_maternal] = rowMedians(betaFINAL, na.rm=TRUE)[k_main_maternal[,1]]
  
summary(dataB12_main_maternal) ## Summarizes data for maternal main model
dim(dataB12_main_maternal) 


## Maternal main model: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataB12_main_maternal$sample.id)
length(index)
beta_main_maternal<-betaFINAL[,index]
beta_main_maternal<-beta_main_maternal[,order(colnames(beta_main_maternal))]
dataB12_main_maternal<-dataB12_main_maternal[order(dataB12_main_maternal$sample.id),]
ncol(beta_main_maternal)
nrow(beta_main_maternal)
all.equal(as.character(dataB12_main_maternal$sample.id),as.character(colnames(beta_main_maternal))) #Return must be TRUE, not FALSE!


betaFINAL_main_maternal <- beta_main_maternal

## Transpose the trimmed methylation data (rows are newborns and columns are CpGs)
beta_matrix_main_maternal<-t(betaFINAL_main_maternal)
dim(beta_matrix_main_maternal)
all.equal(as.character(dataB12_main_maternal$sample.id),rownames(beta_matrix_main_maternal))  #Return must be TRUE, not FALSE!


## Checks if phenotype file and in methylation file into the same order
dataB12_main_maternal <- dataB12_main_maternal[dataB12_main_maternal$sample.id %in% intersect(dataB12_main_maternal$sample.id, rownames(beta_matrix_main_maternal)),]
dataB12_main_maternal <- dataB12_main_maternal[match(rownames(beta_matrix_main_maternal), dataB12_main_maternal$sample.id),]
sum(dataB12_main_maternal$sample.id == rownames(beta_matrix_main_maternal))/nrow(beta_matrix_main_maternal) #check that the IDs are in the same order, should return 1
nprobes_main_maternal = ncol(beta_matrix_main_maternal)

### run EWA for model 1: MATERNAL MAIN MODEL 
results_main_maternal = data.frame(probeID_main_maternal=colnames(beta_matrix_main_maternal),
                 	beta=rep(0, times=nprobes_main_maternal),
                 	se=rep(0, times=nprobes_main_maternal),
                 	p_val=rep(0, times=nprobes_main_maternal),
                 	n=rep(0, times=nprobes_main_maternal))

for(i in 1:nprobes_main_maternal){
  tryCatch({
	CpG_main_maternal = beta_matrix_main_maternal[,i]
	rlm.fit = rlm(CpG_main_maternal ~ mat.b12 +
                	batch +
                	bcell +
                	mono + 
                  cd4t +
                	cd8t + 
                  gran + 
                  nk + 
                  gest.age.sampling +
                  mat.age + 
                  mat.ses +
                  mat.bmi +
                  mat.smoking +
                  parity +
                  sex +
                  child.age,
              	data=dataB12_main_maternal,
              	maxit=200) 
	test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
	beta_main_maternal = test[2,"Estimate"]
	SE_main_maternal = test[2,"Std. Error"]
	PVAL_main_maternal = test[2,"Pr(>|z|)"]
	N_main_maternal = length(rlm.fit$residual)
    
	set(results_main_maternal, i, 2L, beta_main_maternal)
	set(results_main_maternal, i, 3L, SE_main_maternal)
	set(results_main_maternal, i, 4L, PVAL_main_maternal)
	set(results_main_maternal, i, 5L, N_main_maternal)
  }, error = function(err) {
	message("Error in ", colnames(beta_matrix_main_maternal)[i])
	set(results_main_maternal, i, 2L, NA)
	set(results_main_maternal, i, 3L, NA)
	set(results_main_maternal, i, 4L, NA)
	set(results_main_maternal, i, 5L, NA)
  })
 
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_main_maternal, "%\n")}
}

cat("Sorting and saving result file for model 1 (maternal childhood model)...\n\n")
results_main_maternal = results_main_maternal[order(results_main_maternal$p_val),]
cat("Saving the EWAS results for model 1 (maternal childhood model)...\n")

## saves results in .csv file
write.csv(results_main_maternal, file=paste0(cohort,".",maternal,".",childhood.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data:\n\n")
str(results_main_maternal)
cat("\n")
summary(results_main_maternal) # Some summary data to check everything looks good

## Summarise pheno data and save summary as .csv file
main_maternal_B12.tableone <- as.data.frame(print(CreateTableOne(data=dataB12_main_maternal[,-1],factorVars=c("mat.ses", "mat.smoking", "parity", "sex", "batch")),stringsAsFactors=FALSE))
write.csv(main_maternal_B12.tableone,file=paste0(cohort,".",maternal,".",childhood.model,".summary.",analysis.date,".csv"))

print("Maternal main total B12 analysis for childhood/adolescence follow-up analysis completed") 

############################################################## NEWBORN MAIN MODEL - CHILDHOOD/ADOLESCENCE FOLLOW-UP ANALYSIS #######################################################

### IF YOUR COHORT DOES NOT HAVE CORD BLOOD B12 CONCENTRATIONS AVAILABLE, YOU CAN DELETE THIS PART OF THE SCRIPT 


## Newborn main model: remove cases with NA for any of newborn.b12 and main covariates
intersecting.samples <- intersect(dataB12$sample.id,colnames(betaFINAL))
dataB12_main_newborn <- na.omit(dataB12[which(dataB12$sample.id %in% intersecting.samples),unique(c(sample.id, "newborn.b12", batch, cell.names, main.covariates))])
k_main_newborn = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_main_newborn] = rowMedians(betaFINAL, na.rm=TRUE)[k_main_newborn[,1]]
  
summary(dataB12_main_newborn) ### Summarizes data for newborn main model 
dim(dataB12_main_newborn) ###


## newborn main model: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataB12_main_newborn$sample.id)
length(index)
beta_main_newborn<-betaFINAL[,index]
beta_main_newborn<-beta_main_newborn[,order(colnames(beta_main_newborn))]
dataB12_main_newborn<-dataB12_main_newborn[order(dataB12_main_newborn$sample.id),]
ncol(beta_main_newborn)
nrow(beta_main_newborn)
all.equal(as.character(dataB12_main_newborn$sample.id),as.character(colnames(beta_main_newborn))) #Return must be TRUE, not FALSE!


betaFINAL_main_newborn <- beta_main_newborn

## Transpose the trimmed methylation data (rows are newborns and columns are CpGs)
beta_matrix_main_newborn<-t(betaFINAL_main_newborn)
dim(beta_matrix_main_newborn)
all.equal(as.character(dataB12_main_newborn$sample.id),rownames(beta_matrix_main_newborn))  #Return must be TRUE, not FALSE!


## Checks if phenotype file and in methylation file into the same order
dataB12_main_newborn <- dataB12_main_newborn[dataB12_main_newborn$sample.id %in% intersect(dataB12_main_newborn$sample.id, rownames(beta_matrix_main_newborn)),]
dataB12_main_newborn <- dataB12_main_newborn[match(rownames(beta_matrix_main_newborn), dataB12_main_newborn$sample.id),]
sum(dataB12_main_newborn$sample.id == rownames(beta_matrix_main_newborn))/nrow(beta_matrix_main_newborn) #check that the IDs are in the same order, should return 1
nprobes_main_newborn = ncol(beta_matrix_main_newborn)

### run EWA for model 5: NEWBORN MAIN MODEL 
results_main_newborn = data.frame(probeID_main_newborn=colnames(beta_matrix_main_newborn),
                 	beta=rep(0, times=nprobes_main_newborn),
                 	se=rep(0, times=nprobes_main_newborn),
                 	p_val=rep(0, times=nprobes_main_newborn),
                 	n=rep(0, times=nprobes_main_newborn))

for(i in 1:nprobes_main_newborn){
  tryCatch({
	CpG_main_newborn = beta_matrix_main_newborn[,i]
	rlm.fit = rlm(CpG_main_newborn ~ newborn.b12 +
                	batch +
                	bcell +
                	mono + 
                  cd4t +
                	cd8t + 
                  gran + 
                  nk + 
                  mat.age + 
                  mat.ses +
                  mat.bmi +
                  mat.smoking +
                  parity +
                  sex + 
                  child.age,
              	data=dataB12_main_newborn,
              	maxit=200) 
	test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
	beta_main_newborn = test[2,"Estimate"]
	SE_main_newborn = test[2,"Std. Error"]
	PVAL_main_newborn = test[2,"Pr(>|z|)"]
	N_main_newborn = length(rlm.fit$residual)
    
	set(results_main_newborn, i, 2L, beta_main_newborn)
	set(results_main_newborn, i, 3L, SE_main_newborn)
	set(results_main_newborn, i, 4L, PVAL_main_newborn)
	set(results_main_newborn, i, 5L, N_main_newborn)
  }, error = function(err) {
	message("Error in ", colnames(beta_matrix_main_newborn)[i])
	set(results_main_newborn, i, 2L, NA)
	set(results_main_newborn, i, 3L, NA)
	set(results_main_newborn, i, 4L, NA)
	set(results_main_newborn, i, 5L, NA)
  })
 
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_main_newborn, "%\n")}
}

cat("Sorting and saving result file for model 5 (newborn childhood model)...\n\n")
results_main_newborn = results_main_newborn[order(results_main_newborn$p_val),]
cat("Saving the EWAS results for model 5 (newborn childhood model)...\n")

## saves results in .csv file
write.csv(results_main_newborn, file=paste0(cohort,".",newborn,".",childhood.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data:\n\n")
str(results_main_newborn)
cat("\n")
summary(results_main_newborn) # Some summary data to check everything looks good


## Summarise pheno data and save summaries as .csv files
main_newborn_B12.tableone <- as.data.frame(print(CreateTableOne(data=dataB12_main_newborn[,-1],factorVars=c("mat.ses", "mat.smoking", "parity", "sex", "batch")),stringsAsFactors=FALSE))
write.csv(main_newborn_B12.tableone,file=paste0(cohort,".",newborn,".",childhood.model,".summary.",analysis.date,".csv"))

print("Newborn main b12 analysis for childhood/adolescence follow-up analysis completed") ## This completes the analysis of model 5. If you are able to run model 6, please do not delete any code from here.


####################################################################### COMPLETE ANALYSES #############################################################

print("All analyses completed")

