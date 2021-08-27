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
analysis.date <- "20200224" # change to the date at which you perform the analyses


## You do not have to chagne anything in this part of the script, even if you are not able to run one or multiple models ##
maternal <- "maternal" # refers to mothers in who vitamin B12 concentrations were measured in pregnancy, if available in your study
newborn <- "newborn" # refers to newborns in who vitamin B12 concentrations were measured, if available in your study
main.model <- "main.model"
folate.model <- "folate.model"
homocysteine.model <- "homocysteine.model"
active.model <- "active.model"


############################################################## MODELS ###################################################################

## If you are not able to run one or multiple models, please delete the R code for this/these model(s) further below ##

############################################################## SET MODEL PARAMETERS #######################################################

## Please make sure these names are identical to the names in your .dat phenotype file! ###
## Here, you only have to change anything in the script if you want to add selection factors (and include these in your phenotype file as well)


sample.id <-c("sample.id") ### Case identifier
traits <-c("mat.b12", "mat.active.b12", "newborn.b12", "newborn.active.b12") ### The exposure variables of interest
cell.names <- c("bcell", "mono", "cd4t", "cd8t", "gran", "nk", "nRBC") ### Salas reference set for cell type correction
batch <-c("batch") ### Adjust for batch effects by including the most important covariate(s), such as plate, in the models
sampling <-c("gest.age.sampling")
main.covariates <- c("mat.age", "mat.ses", "mat.bmi", "mat.smoking", "parity", "sex") ### Will be included in each model
folate <- c("mat.folate", "newborn.folate") ### For sensitivity analyses: Mat.folate will be included in the maternal folate model; newborn.folate will be included in the newborn folate model.
homocysteine <- c("mat.homocysteine", "newborn.homocysteine") ### For sensitivity analyses: Mat.homocysteine will be included in the maternal folate model; newborn.homocysteine will be included in the newborn folate model.



############################################################## PHENOTYPE FILE #######################################################

## Load phenotype data
dataB12 <- read.table("EWA_PACE_B12.dat", header = T, stringsAsFactors=FALSE) ### between quotes, specify the name of your phenotype file
names(dataB12)

## dataB12 is a .dat file (tab delimited) of 24 columns containing information on: sample IDs (1st column), 4 exposures (vitamin b12 concentrations of mothers and newborns), cell counts (Salas reference, 7 columns), gestational age at sampling, maternal covariates, parity, child sex, folate concentrations, mocysteine concentration (2x2 columns).
## There should be no extra columns, except if you have added selection factors. If this is the case, for questions about how to change the R script, plese contact g.monasso@erasmusmc.nl

## This part is to make sure that categorical covariates are treated as factors and continuous covariates as numeric.
dataB12$mat.b12 = as.numeric(dataB12$mat.b12)
dataB12$mat.active.b12 = as.numeric(dataB12$mat.active.b12)
dataB12$newborn.b12 = as.numeric(dataB12$newborn.b12)
dataB12$newborn.active.b12 = as.numeric(dataB12$newborn.active.b12)
dataB12$batch = as.factor(dataB12$batch)
dataB12$bcell = as.numeric(dataB12$bcell)
dataB12$mono = as.numeric(dataB12$mono)
dataB12$cd4t = as.numeric(dataB12$cd4t)
dataB12$cd8t = as.numeric(dataB12$cd8t)
dataB12$gran = as.numeric(dataB12$gran)
dataB12$nk = as.numeric(dataB12$nk)
dataB12$nRBC = as.numeric(dataB12$nRBC)
dataB12$gest.age.sampling = as.numeric(dataB12$gest.age.sampling)
dataB12$mat.age = as.numeric(dataB12$mat.age)
dataB12$mat.ses = as.factor(dataB12$mat.ses)
dataB12$mat.bmi = as.numeric(dataB12$mat.bmi)
dataB12$mat.smoking = as.factor(dataB12$mat.smoking)
dataB12$parity = as.factor(dataB12$parity)
dataB12$sex = as.factor(dataB12$sex)
dataB12$mat.folate = as.numeric(dataB12$mat.folate)
dataB12$newborn.folate = as.numeric(dataB12$newborn.folate)
dataB12$mat.homocysteine = as.numeric(dataB12$mat.homocysteine)
dataB12$newborn.homocysteine = as.numeric(dataB12$newborn.homocysteine)


## Summarizes your data before excluding any NA's
summary(dataB12)


## Check if all variables are present in phenotype file
for(i in 1:length(c(sample.id, traits, batch, cell.names, sampling, main.covariates, folate, homocysteine))) {
  print(ifelse(c(sample.id, traits, batch, cell.names, sampling, main.covariates, folate, homocysteine)[i] %in% colnames(dataB12)==FALSE,
               paste("CAUTION: the variable called",c(sample.id, traits, batch, cell.names, sampling, main.covariates, folate, homocysteine)[i],"is missing from pheno file"),
               paste("variable called",c(sample.id, traits, batch, cell.names, sampling, main.covariates, folate, homocysteine)[i],"is present in pheno file")))
}



############################################################## METHYLATION FILE #######################################################

load("...") ## load methylation file, eg. load("3IQR-trimmed-methylationfile.Rdata") 
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



################### FROM HERE DELETE THOSE PARTS OF THE SCRIPT THAT ARE MEANT TO RUN MODELS YOU CANNOT RUN ###################



############################################################## MODEL 1: MATERNAL MAIN MODEL #######################################################

## IF NOT AVAILABLE IN YOUR COHORT, YOU CAN DELETE THIS PART OF THE SCRIPT (UNTIL ##### MODEL 4: MATERNAL ACTIVE B12 MODEL ####)

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
                  nRBC +
                  gest.age.sampling +
                  mat.age + 
                  mat.ses +
                  mat.bmi +
                  mat.smoking +
                  parity +
                  sex,
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

cat("Sorting and saving result file for model 1 (maternal main model)...\n\n")
results_main_maternal = results_main_maternal[order(results_main_maternal$p_val),]
cat("Saving the EWAS results for model 1 (maternal main model)...\n")

## saves results in .csv file
write.csv(results_main_maternal, file=paste0(cohort,".",maternal,".",main.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data:\n\n")
str(results_main_maternal)
cat("\n")
summary(results_main_maternal) # Some summary data to check everything looks good

## Calculate the epigenetic inflation factor lambda and save it in a file:
lambda_main_maternal = median(qchisq(results_main_maternal$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for model 1 (main maternal model) is: ", lambda_main_maternal, "\n")
write.table(lambda_main_maternal, file=paste0(cohort,".",maternal,".",main.model,".lambda.",analysis.date,".txt"), quote=FALSE, row.names=F)

## Summarise pheno data and save summary as .csv file
main_maternal_B12.tableone <- as.data.frame(print(CreateTableOne(data=dataB12_main_maternal[,-1],factorVars=c("mat.ses", "mat.smoking", "parity", "sex", "batch")),stringsAsFactors=FALSE))
write.csv(main_maternal_B12.tableone,file=paste0(cohort,".",maternal,".",main.model,".summary.",analysis.date,".csv"))

print("Model 1: Maternal main total b12 analysis completed") ## This completes the analysis of model 1. If you are able to run model 2, please do not delete any code from here



############################################################## MODEL 2: MATERNAL FOLATE MODEL #######################################################

### IF NOT AVAILABLE IN YOUR COHORT, YOU CAN DELETE THIS PART OF THE SCRIPT (UNTIL ##### MODEL 3: MATERNAL HOMOCYSTEINE MODEL ####)


### maternal folate model: remove cases with NA for any of b12 and covariates including mat.folate
intersecting.samples <- intersect(dataB12$sample.id,colnames(betaFINAL))
dataB12_folate_maternal <- na.omit(dataB12[which(dataB12$sample.id %in% intersecting.samples),unique(c(sample.id, "mat.b12", batch, cell.names, sampling, main.covariates, "mat.folate"))])
k_folate_maternal = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_folate_maternal] = rowMedians(betaFINAL, na.rm=TRUE)[k_folate_maternal[,1]]


summary(dataB12_folate_maternal)## Summarizes data for maternal folate model
dim(dataB12_folate_maternal) 


## Maternal folate model: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataB12_folate_maternal$sample.id)
length(index)
beta_folate_maternal<-betaFINAL[,index]
beta_folate_maternal<-beta_folate_maternal[,order(colnames(beta_folate_maternal))]
dataB12_folate_maternal<-dataB12_folate_maternal[order(dataB12_folate_maternal$sample.id),]
ncol(beta_folate_maternal)
nrow(beta_folate_maternal)
all.equal(as.character(dataB12_folate_maternal$sample.id),as.character(colnames(beta_folate_maternal))) #Return must be TRUE, not FALSE!


betaFINAL_folate_maternal<-beta_folate_maternal


## Transpose the trimmed methylation data (rows are newborns and columns are CpGs)
beta_matrix_folate_maternal<-t(betaFINAL_folate_maternal)
dim(beta_matrix_folate_maternal)
all.equal(as.character(dataB12_folate_maternal$sample.id),rownames(beta_matrix_folate_maternal))  #Return must be TRUE, not FALSE!


#Put the IDs in phenotype_file and in methylation_file into the same order
dataB12_folate_maternal <- dataB12_folate_maternal[dataB12_folate_maternal$sample.id %in% intersect(dataB12_folate_maternal$sample.id, rownames(beta_matrix_folate_maternal)),]
dataB12_folate_maternal <- dataB12_folate_maternal[match(rownames(beta_matrix_folate_maternal), dataB12_folate_maternal$sample.id),]
sum(dataB12_folate_maternal$sample.id == rownames(beta_matrix_folate_maternal))/nrow(beta_matrix_folate_maternal) #check that the IDs are in the same order, should return 1
nprobes_folate_maternal = ncol(beta_matrix_folate_maternal)


### run EWA for MODEL 2: MATERNAL FOLATE MODEL 
results_folate_maternal = data.frame(probeID_folate_maternal=colnames(beta_matrix_folate_maternal),
                 	beta=rep(0, times=nprobes_folate_maternal),
                 	se=rep(0, times=nprobes_folate_maternal),
                 	p_val=rep(0, times=nprobes_folate_maternal),
                 	n=rep(0, times=nprobes_folate_maternal))

for(i in 1:nprobes_folate_maternal){
  tryCatch({
	CpG_folate_maternal = beta_matrix_folate_maternal[,i]
	rlm.fit = rlm(CpG_folate_maternal ~ mat.b12 +
                  batch +
                	bcell +
                	mono + 
                  cd4t +
                	cd8t + 
                  gran + 
                  nk + 
                  nRBC +
                  gest.age.sampling +
                  mat.age + 
                  mat.ses +
                  mat.bmi +
                  mat.smoking +
                  parity +
                  sex +
                  mat.folate,
              	data=dataB12_folate_maternal,
              	maxit=200) 
	test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
	BETA_folate_maternal = test[2,"Estimate"]
	SE_folate_maternal = test[2,"Std. Error"]
	PVAL_folate_maternal = test[2,"Pr(>|z|)"]
	N_folate_maternal = length(rlm.fit$residual)
    
	set(results_folate_maternal, i, 2L, BETA_folate_maternal)
	set(results_folate_maternal, i, 3L, SE_folate_maternal)
	set(results_folate_maternal, i, 4L, PVAL_folate_maternal)
	set(results_folate_maternal, i, 5L, N_folate_maternal)
  }, error = function(err) {
	message("Error in ", colnames(beta_matrix_folate_maternal)[i])
	set(results_folate_maternal, i, 2L, NA)
	set(results_folate_maternal, i, 3L, NA)
	set(results_folate_maternal, i, 4L, NA)
	set(results_folate_maternal, i, 5L, NA)
  })
 
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_folate_maternal, "%\n")}
}


cat("Sorting and saving result file for maternal folate model 2...\n\n")
results_folate_maternal = results_folate_maternal[order(results_folate_maternal$p_val),]
cat("Saving the EWAS results for maternal folate model 2...\n")

## saves results in .csv file
write.csv(results_folate_maternal, file=paste0(cohort,".",maternal,".",folate.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)


cat("Done\n\n")
cat("Checking and summary data for maternal folate model 2:\n\n")
str(results_folate_maternal)
cat("\n")
summary(results_folate_maternal) #some summary data to check everything looks good


#Calculate the epigenetic inflation factor lambda for maternal folate model 2:
lambda_folate_maternal = median(qchisq(results_folate_maternal$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for the maternal folate model 2 is: ", lambda_folate_maternal, "\n")
write.table(lambda_folate_maternal, file=paste0(cohort,".",maternal,".",folate.model,".lambda.",analysis.date,".txt"), col.names=T, row.names=F, quote=F)


# Summarise pheno data and save summary as .csv file
folate_maternal_B12.tableone <- as.data.frame(print(CreateTableOne(data=dataB12_folate_maternal[,-1],factorVars=c("mat.ses", "mat.smoking", "parity", "sex", "batch")),stringsAsFactors=FALSE))
write.csv(folate_maternal_B12.tableone,file=paste0(cohort,".",maternal,".",folate.model,".summary.",analysis.date,".csv"))

print("Model 2: Maternal folate analysis completed") ## This completes the analysis of model 2. If you are able to run model 2, please do not delete any code from here.




############################################################## MODEL 3: MATERNAL HOMOCYSTEINE MODEL #######################################################

### IF NOT AVAILABLE IN YOUR COHORT, YOU CAN DELETE THIS PART OF THE SCRIPT (UNTIL ##### MODEL 4: MATERNAL ACTIVE MODEL ####)


### maternal homocysteine model: remove cases with NA for any of b12 and covariates including mat.homocysteine
intersecting.samples <- intersect(dataB12$sample.id,colnames(betaFINAL))
dataB12_homocysteine_maternal <- na.omit(dataB12[which(dataB12$sample.id %in% intersecting.samples),unique(c(sample.id, "mat.b12", batch, cell.names, sampling, main.covariates, "mat.homocysteine"))])
k_homocysteine_maternal = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_homocysteine_maternal] = rowMedians(betaFINAL, na.rm=TRUE)[k_homocysteine_maternal[,1]]


summary(dataB12_homocysteine_maternal)## Summarizes data for maternal homocysteine model
dim(dataB12_homocysteine_maternal) 


## Maternal homocysteine model: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataB12_homocysteine_maternal$sample.id)
length(index)
beta_homocysteine_maternal<-betaFINAL[,index]
beta_homocysteine_maternal<-beta_homocysteine_maternal[,order(colnames(beta_homocysteine_maternal))]
dataB12_homocysteine_maternal<-dataB12_homocysteine_maternal[order(dataB12_homocysteine_maternal$sample.id),]
ncol(beta_homocysteine_maternal)
nrow(beta_homocysteine_maternal)
all.equal(as.character(dataB12_homocysteine_maternal$sample.id),as.character(colnames(beta_homocysteine_maternal))) #Return must be TRUE, not FALSE!


betaFINAL_homocysteine_maternal<-beta_homocysteine_maternal


## Transpose the trimmed methylation data (rows are newborns and columns are CpGs)
beta_matrix_homocysteine_maternal<-t(betaFINAL_homocysteine_maternal)
dim(beta_matrix_homocysteine_maternal)
all.equal(as.character(dataB12_homocysteine_maternal$sample.id),rownames(beta_matrix_homocysteine_maternal))  #Return must be TRUE, not FALSE!


#Put the IDs in phenotype_file and in methylation_file into the same order
dataB12_homocysteine_maternal <- dataB12_homocysteine_maternal[dataB12_homocysteine_maternal$sample.id %in% intersect(dataB12_homocysteine_maternal$sample.id, rownames(beta_matrix_homocysteine_maternal)),]
dataB12_homocysteine_maternal <- dataB12_homocysteine_maternal[match(rownames(beta_matrix_homocysteine_maternal), dataB12_homocysteine_maternal$sample.id),]
sum(dataB12_homocysteine_maternal$sample.id == rownames(beta_matrix_homocysteine_maternal))/nrow(beta_matrix_homocysteine_maternal) #check that the IDs are in the same order, should return 1
nprobes_homocysteine_maternal = ncol(beta_matrix_homocysteine_maternal)


### run EWA for MODEL 3: MATERNAL HOMOCYSTEINE MODEL 
results_homocysteine_maternal = data.frame(probeID_homocysteine_maternal=colnames(beta_matrix_homocysteine_maternal),
                 	beta=rep(0, times=nprobes_homocysteine_maternal),
                 	se=rep(0, times=nprobes_homocysteine_maternal),
                 	p_val=rep(0, times=nprobes_homocysteine_maternal),
                 	n=rep(0, times=nprobes_homocysteine_maternal))

for(i in 1:nprobes_homocysteine_maternal){
  tryCatch({
	CpG_homocysteine_maternal = beta_matrix_homocysteine_maternal[,i]
	rlm.fit = rlm(CpG_homocysteine_maternal ~ mat.b12 +
                  batch +
                	bcell +
                	mono + 
                  cd4t +
                	cd8t + 
                  gran + 
                  nk + 
                  nRBC +
                  gest.age.sampling +
                  mat.age + 
                  mat.ses +
                  mat.bmi +
                  mat.smoking +
                  parity +
                  sex +
                  mat.homocysteine,
              	data=dataB12_homocysteine_maternal,
              	maxit=200) 
	test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
	BETA_homocysteine_maternal = test[2,"Estimate"]
	SE_homocysteine_maternal = test[2,"Std. Error"]
	PVAL_homocysteine_maternal = test[2,"Pr(>|z|)"]
	N_homocysteine_maternal = length(rlm.fit$residual)
    
	set(results_homocysteine_maternal, i, 2L, BETA_homocysteine_maternal)
	set(results_homocysteine_maternal, i, 3L, SE_homocysteine_maternal)
	set(results_homocysteine_maternal, i, 4L, PVAL_homocysteine_maternal)
	set(results_homocysteine_maternal, i, 5L, N_homocysteine_maternal)
  }, error = function(err) {
	message("Error in ", colnames(beta_matrix_homocysteine_maternal)[i])
	set(results_homocysteine_maternal, i, 2L, NA)
	set(results_homocysteine_maternal, i, 3L, NA)
	set(results_homocysteine_maternal, i, 4L, NA)
	set(results_homocysteine_maternal, i, 5L, NA)
  })
 
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_homocysteine_maternal, "%\n")}
}


cat("Sorting and saving result file for maternal homocysteine model 3...\n\n")
results_homocysteine_maternal = results_homocysteine_maternal[order(results_homocysteine_maternal$p_val),]
cat("Saving the EWAS results for maternal homocysteine model 3...\n")

## saves results in .csv file
write.csv(results_homocysteine_maternal, file=paste0(cohort,".",maternal,".",homocysteine.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)


cat("Done\n\n")
cat("Checking and summary data for maternal homocysteine model 3:\n\n")
str(results_homocysteine_maternal)
cat("\n")
summary(results_homocysteine_maternal) #some summary data to check everything looks good


#Calculate the epigenetic inflation factor lambda for maternal homocysteine model 3:
lambda_homocysteine_maternal = median(qchisq(results_homocysteine_maternal$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for the maternal homocysteine model 3 is: ", lambda_homocysteine_maternal, "\n")
write.table(lambda_homocysteine_maternal, file=paste0(cohort,".",maternal,".",homocysteine.model,".lambda.",analysis.date,".txt"), col.names=T, row.names=F, quote=F)


# Summarise pheno data and save summary as .csv file
homocysteine_maternal_B12.tableone <- as.data.frame(print(CreateTableOne(data=dataB12_homocysteine_maternal[,-1],factorVars=c("mat.ses", "mat.smoking", "parity", "sex", "batch")),stringsAsFactors=FALSE))
write.csv(homocysteine_maternal_B12.tableone,file=paste0(cohort,".",maternal,".",homocysteine.model,".summary.",analysis.date,".csv"))

print("Model 3: Maternal homocysteine analysis completed") ## This completes the analysis of model 3. If you are able to run model 4, please do not delete any code from here.



############################################################## MODEL 4: MATERNAL ACTIVE MODEL #######################################################

### IF NOT AVAILABLE IN YOUR COHORT, YOU CAN DELETE THIS PART OF THE SCRIPT (UNTIL ##### MODEL 5: NEWBORN MAIN MODEL ####)


### maternal active B12 model: remove cases with NA for any of maternal active B12 and main covariates
intersecting.samples <- intersect(dataB12$sample.id,colnames(betaFINAL))
dataB12_active_maternal <- na.omit(dataB12[which(dataB12$sample.id %in% intersecting.samples),unique(c(sample.id, "mat.active.b12", batch, cell.names, sampling, main.covariates))])
k_active_maternal = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_active_maternal] = rowMedians(betaFINAL, na.rm=TRUE)[k_active_maternal[,1]]

summary(dataB12_active_maternal)## Summarizes data for maternal active b12 model
dim(dataB12_active_maternal) 


## Maternal active b12 model 4: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataB12_active_maternal$sample.id)
length(index)
beta_active_maternal<-betaFINAL[,index]
beta_active_maternal<-beta_active_maternal[,order(colnames(beta_active_maternal))]
dataB12_active_maternal<-dataB12_active_maternal[order(dataB12_active_maternal$sample.id),]
ncol(beta_active_maternal)
nrow(beta_active_maternal)
all.equal(as.character(dataB12_active_maternal$sample.id),as.character(colnames(beta_active_maternal))) #Return must be TRUE, not FALSE!

betaFINAL_active_maternal<-beta_active_maternal


## Transpose the trimmed methylation data (rows are newborns and columns are CpGs)
beta_matrix_active_maternal<-t(betaFINAL_active_maternal)
dim(beta_matrix_active_maternal)
all.equal(as.character(dataB12_active_maternal$sample.id),rownames(beta_matrix_active_maternal))  #Return must be TRUE, not FALSE!


#Put the IDs in phenotype file and in methylation file into the same order
dataB12_active_maternal <- dataB12_active_maternal[dataB12_active_maternal$sample.id %in% intersect(dataB12_active_maternal$sample.id, rownames(beta_matrix_active_maternal)),]
dataB12_active_maternal <- dataB12_active_maternal[match(rownames(beta_matrix_active_maternal), dataB12_active_maternal$sample.id),]
sum(dataB12_active_maternal$sample.id == rownames(beta_matrix_active_maternal))/nrow(beta_matrix_active_maternal) #check that the IDs are in the same order, should return 1
nprobes_active_maternal = ncol(beta_matrix_active_maternal)


### run EWA for MODEL 4: MATERNAL ACTIVE MODEL 
results_active_maternal = data.frame(probeID_active_maternal=colnames(beta_matrix_active_maternal),
                 	beta=rep(0, times=nprobes_active_maternal),
                 	se=rep(0, times=nprobes_active_maternal),
                 	p_val=rep(0, times=nprobes_active_maternal),
                 	n=rep(0, times=nprobes_active_maternal))

for(i in 1:nprobes_active_maternal){
  tryCatch({
	CpG_active_maternal = beta_matrix_active_maternal[,i]
	rlm.fit = rlm(CpG_active_maternal ~ mat.active.b12 +
                	batch +
                	bcell +
                	mono + 
                  cd4t +
                	cd8t + 
                  gran + 
                  nk + 
                  nRBC +
                  gest.age.sampling +
                  mat.age + 
                  mat.ses +
                  mat.bmi +
                  mat.smoking +
                  parity +
                  sex,
              	data=dataB12_active_maternal,
              	maxit=200) 
	test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
	BETA_active_maternal = test[2,"Estimate"]
	SE_active_maternal = test[2,"Std. Error"]
	PVAL_active_maternal = test[2,"Pr(>|z|)"]
	N_active_maternal = length(rlm.fit$residual)
    
	set(results_active_maternal, i, 2L, BETA_active_maternal)
	set(results_active_maternal, i, 3L, SE_active_maternal)
	set(results_active_maternal, i, 4L, PVAL_active_maternal)
	set(results_active_maternal, i, 5L, N_active_maternal)
  }, error = function(err) {
	message("Error in ", colnames(beta_matrix_active_maternal)[i])
	set(results_active_maternal, i, 2L, NA)
	set(results_active_maternal, i, 3L, NA)
	set(results_active_maternal, i, 4L, NA)
	set(results_active_maternal, i, 5L, NA)
  })
 
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_active_maternal, "%\n")}
}


cat("Sorting and saving result file for maternal active b12 model 4 ...\n\n")
results_active_maternal = results_active_maternal[order(results_active_maternal$p_val),] 
cat("Saving the EWAS results for maternal active b12 model 4...\n")

## saves results in .csv file
write.csv(results_active_maternal, file=paste0(cohort,".",maternal,".",active.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data for maternal active b12 model 4:\n\n")
str(results_active_maternal)
cat("\n")
summary(results_active_maternal) #some summary data to check everything looks good

#Calculate the epigenetic inflation factor lambda for maternal active b12 model 4:
lambda_active_maternal = median(qchisq(results_active_maternal$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for maternal active b12 model 4 is: ", lambda_active_maternal, "\n")
write.table(lambda_active_maternal, file=paste0(cohort,".",maternal,".",active.model,".lambda.",analysis.date,".txt"), col.names=T, row.names=F, quote=F)


# Summarise pheno data and save summary as .csv file
active_maternal_B12.tableone <- as.data.frame(print(CreateTableOne(data=dataB12_active_maternal[,-1],factorVars=c("mat.ses", "mat.smoking", "parity", "sex", "batch")),stringsAsFactors=FALSE))
write.csv(active_maternal_B12.tableone,file=paste0(cohort,".",maternal,".",active.model,".summary.",analysis.date,".csv"))

print("Model 4: maternal active B12 analysis completed") ## This completes the analysis of model 4 If you are able to run model 5, please do not delete any code from here.



############################################################## MODEL 5: NEWBORN MAIN MODEL #######################################################

### IF NOT AVAILABLE IN YOUR COHORT, YOU CAN DELETE THIS PART OF THE SCRIPT (UNTIL ##### MODEL 8: NEWBORN ACTIVE MODEL ####)


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
                  nRBC +
                  mat.age + 
                  mat.ses +
                  mat.bmi +
                  mat.smoking +
                  parity +
                  sex,
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

cat("Sorting and saving result file for model 5 (newborn main model)...\n\n")
results_main_newborn = results_main_newborn[order(results_main_newborn$p_val),]
cat("Saving the EWAS results for model 5 (newborn main model)...\n")

## saves results in .csv file
write.csv(results_main_newborn, file=paste0(cohort,".",newborn,".",main.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data:\n\n")
str(results_main_newborn)
cat("\n")
summary(results_main_newborn) # Some summary data to check everything looks good

## Calculate the epigenetic inflation factor lambda and save it in a txt file:
lambda_main_newborn = median(qchisq(results_main_newborn$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for model 5 (main newborn model) is: ", lambda_main_newborn, "\n")
write.table(lambda_main_newborn, file=paste0(cohort,".",newborn,".",main.model,".lambda.",analysis.date,".txt"), quote=FALSE, row.names=F)

## Summarise pheno data and save summaries as .csv files
main_newborn_B12.tableone <- as.data.frame(print(CreateTableOne(data=dataB12_main_newborn[,-1],factorVars=c("mat.ses", "mat.smoking", "parity", "sex", "batch")),stringsAsFactors=FALSE))
write.csv(main_newborn_B12.tableone,file=paste0(cohort,".",newborn,".",main.model,".summary.",analysis.date,".csv"))

print("Model 5: newborn main b12 analysis completed") ## This completes the analysis of model 5. If you are able to run model 6, please do not delete any code from here.



############################################################## MODEL 6: NEWBORN FOLATE MODEL #######################################################

### IF NOT AVAILABLE IN YOUR COHORT, YOU CAN DELETE THIS PART OF THE SCRIPT (UNTIL ##### MODEL 7: NEWBORN HOMOCYSTEINE MODEL ####)


## Newborn folate model: remove cases with NA for any of newborn.b12 and covariates including folate
intersecting.samples <- intersect(dataB12$sample.id,colnames(betaFINAL))
dataB12_folate_newborn <- na.omit(dataB12[which(dataB12$sample.id %in% intersecting.samples),unique(c(sample.id, "newborn.b12", batch, cell.names, main.covariates, "newborn.folate"))])
k_folate_newborn = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_folate_newborn] = rowMedians(betaFINAL, na.rm=TRUE)[k_folate_newborn[,1]]
  
summary(dataB12_folate_newborn) ### Summarizes data for newborn folate model
dim(dataB12_folate_newborn) 


## newborn folate model: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataB12_folate_newborn$sample.id)
length(index)
beta_folate_newborn<-betaFINAL[,index]
beta_folate_newborn<-beta_folate_newborn[,order(colnames(beta_folate_newborn))]
dataB12_folate_newborn<-dataB12_folate_newborn[order(dataB12_folate_newborn$sample.id),]
ncol(beta_folate_newborn)
nrow(beta_folate_newborn)
all.equal(as.character(dataB12_folate_newborn$sample.id),as.character(colnames(beta_folate_newborn))) #Return must be TRUE, not FALSE!


betaFINAL_folate_newborn <- beta_folate_newborn

## Transpose the trimmed methylation data (rows are newborns and columns are CpGs)
beta_matrix_folate_newborn<-t(betaFINAL_folate_newborn)
dim(beta_matrix_folate_newborn)
all.equal(as.character(dataB12_folate_newborn$sample.id),rownames(beta_matrix_folate_newborn))  #Return must be TRUE, not FALSE!


## Checks if phenotype file and in methylation file into the same order
dataB12_folate_newborn <- dataB12_folate_newborn[dataB12_folate_newborn$sample.id %in% intersect(dataB12_folate_newborn$sample.id, rownames(beta_matrix_folate_newborn)),]
dataB12_folate_newborn <- dataB12_folate_newborn[match(rownames(beta_matrix_folate_newborn), dataB12_folate_newborn$sample.id),]
sum(dataB12_folate_newborn$sample.id == rownames(beta_matrix_folate_newborn))/nrow(beta_matrix_folate_newborn) #check that the IDs are in the same order, should return 1
nprobes_folate_newborn = ncol(beta_matrix_folate_newborn)

### run EWA for model 6: NEWBORN FOLATE MODEL
results_folate_newborn = data.frame(probeID_folate_newborn=colnames(beta_matrix_folate_newborn),
                 	beta=rep(0, times=nprobes_folate_newborn),
                 	se=rep(0, times=nprobes_folate_newborn),
                 	p_val=rep(0, times=nprobes_folate_newborn),
                 	n=rep(0, times=nprobes_folate_newborn))

for(i in 1:nprobes_folate_newborn){
  tryCatch({
	CpG_folate_newborn = beta_matrix_folate_newborn[,i]
	rlm.fit = rlm(CpG_folate_newborn ~ newborn.b12 +
                	batch +
                	bcell +
                	mono + 
                  cd4t +
                	cd8t + 
                  gran + 
                  nk + 
                  nRBC +
                  mat.age + 
                  mat.ses +
                  mat.bmi +
                  mat.smoking +
                  parity +
                  sex +
                  newborn.folate,
              	data=dataB12_folate_newborn,
              	maxit=200) 
	test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
	beta_folate_newborn = test[2,"Estimate"]
	SE_folate_newborn = test[2,"Std. Error"]
	PVAL_folate_newborn = test[2,"Pr(>|z|)"]
	N_folate_newborn = length(rlm.fit$residual)
    
	set(results_folate_newborn, i, 2L, beta_folate_newborn)
	set(results_folate_newborn, i, 3L, SE_folate_newborn)
	set(results_folate_newborn, i, 4L, PVAL_folate_newborn)
	set(results_folate_newborn, i, 5L, N_folate_newborn)
  }, error = function(err) {
	message("Error in ", colnames(beta_matrix_folate_newborn)[i])
	set(results_folate_newborn, i, 2L, NA)
	set(results_folate_newborn, i, 3L, NA)
	set(results_folate_newborn, i, 4L, NA)
	set(results_folate_newborn, i, 5L, NA)
  })
 
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_folate_newborn, "%\n")}
}

cat("Sorting and saving result file for model 6 (newborn folate model)...\n\n")
results_folate_newborn = results_folate_newborn[order(results_folate_newborn$p_val),]
cat("Saving the EWAS results for model 6 (newborn folate model)...\n")

## saves results in .csv file
write.csv(results_folate_newborn, file=paste0(cohort,".",newborn,".",folate.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data:\n\n")
str(results_folate_newborn)
cat("\n")
summary(results_folate_newborn) # Some summary data to check everything looks good

## Calculate the epigenetic inflation factor lambda and save it in a txt file:
lambda_folate_newborn = median(qchisq(results_folate_newborn$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for model 6 (folate newborn model) is: ", lambda_folate_newborn, "\n")
write.table(lambda_folate_newborn, file=paste0(cohort,".",newborn,".",folate.model,".lambda.",analysis.date,".txt"), quote=FALSE, row.names=F)

## Summarise pheno data and save summary as .csv file
folate_newborn_B12.tableone <- as.data.frame(print(CreateTableOne(data=dataB12_folate_newborn[,-1],factorVars=c("mat.ses", "mat.smoking", "parity", "sex", "batch")),stringsAsFactors=FALSE))
write.csv(folate_newborn_B12.tableone,file=paste0(cohort,".",newborn,".",folate.model,".summary.",analysis.date,".csv"))

print("Model 6: newborn folate b12 analysis completed") ## This completes the analysis of model 6. If you are able to run model 7, please do not delete any code from here.



############################################################## MODEL 7: NEWBORN HOMOCYSTEINE MODEL #######################################################

### IF NOT AVAILABLE IN YOUR COHORT, YOU CAN DELETE THIS PART OF THE SCRIPT (UNTIL ##### MODEL 8: NEWBORN ACTIVE MODEL ####)


## Newborn homocysteine model: remove cases with NA for any of newborn.b12 and covariates including homocysteine
intersecting.samples <- intersect(dataB12$sample.id,colnames(betaFINAL))
dataB12_homocysteine_newborn <- na.omit(dataB12[which(dataB12$sample.id %in% intersecting.samples),unique(c(sample.id, "newborn.b12", batch, cell.names, main.covariates, "newborn.homocysteine"))])
k_homocysteine_newborn = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_homocysteine_newborn] = rowMedians(betaFINAL, na.rm=TRUE)[k_homocysteine_newborn[,1]]
  
summary(dataB12_homocysteine_newborn) ### Summarizes data for newborn homocysteine model
dim(dataB12_homocysteine_newborn) 


## newborn homocysteine model: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataB12_homocysteine_newborn$sample.id)
length(index)
beta_homocysteine_newborn<-betaFINAL[,index]
beta_homocysteine_newborn<-beta_homocysteine_newborn[,order(colnames(beta_homocysteine_newborn))]
dataB12_homocysteine_newborn<-dataB12_homocysteine_newborn[order(dataB12_homocysteine_newborn$sample.id),]
ncol(beta_homocysteine_newborn)
nrow(beta_homocysteine_newborn)
all.equal(as.character(dataB12_homocysteine_newborn$sample.id),as.character(colnames(beta_homocysteine_newborn))) #Return must be TRUE, not FALSE!


betaFINAL_homocysteine_newborn <- beta_homocysteine_newborn

## Transpose the trimmed methylation data (rows are newborns and columns are CpGs)
beta_matrix_homocysteine_newborn<-t(betaFINAL_homocysteine_newborn)
dim(beta_matrix_homocysteine_newborn)
all.equal(as.character(dataB12_homocysteine_newborn$sample.id),rownames(beta_matrix_homocysteine_newborn))  #Return must be TRUE, not FALSE!


## Checks if phenotype file and in methylation file into the same order
dataB12_homocysteine_newborn <- dataB12_homocysteine_newborn[dataB12_homocysteine_newborn$sample.id %in% intersect(dataB12_homocysteine_newborn$sample.id, rownames(beta_matrix_homocysteine_newborn)),]
dataB12_homocysteine_newborn <- dataB12_homocysteine_newborn[match(rownames(beta_matrix_homocysteine_newborn), dataB12_homocysteine_newborn$sample.id),]
sum(dataB12_homocysteine_newborn$sample.id == rownames(beta_matrix_homocysteine_newborn))/nrow(beta_matrix_homocysteine_newborn) #check that the IDs are in the same order, should return 1
nprobes_homocysteine_newborn = ncol(beta_matrix_homocysteine_newborn)

### run EWA for model 7: NEWBORN HOMOCYSTEINE MODEL
results_homocysteine_newborn = data.frame(probeID_homocysteine_newborn=colnames(beta_matrix_homocysteine_newborn),
                 	beta=rep(0, times=nprobes_homocysteine_newborn),
                 	se=rep(0, times=nprobes_homocysteine_newborn),
                 	p_val=rep(0, times=nprobes_homocysteine_newborn),
                 	n=rep(0, times=nprobes_homocysteine_newborn))

for(i in 1:nprobes_homocysteine_newborn){
  tryCatch({
	CpG_homocysteine_newborn = beta_matrix_homocysteine_newborn[,i]
	rlm.fit = rlm(CpG_homocysteine_newborn ~ newborn.b12 +
                	batch +
                	bcell +
                	mono + 
                  cd4t +
                	cd8t + 
                  gran + 
                  nk + 
                  nRBC +
                  mat.age + 
                  mat.ses +
                  mat.bmi +
                  mat.smoking +
                  parity +
                  sex +
                  newborn.homocysteine,
              	data=dataB12_homocysteine_newborn,
              	maxit=200) 
	test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
	beta_homocysteine_newborn = test[2,"Estimate"]
	SE_homocysteine_newborn = test[2,"Std. Error"]
	PVAL_homocysteine_newborn = test[2,"Pr(>|z|)"]
	N_homocysteine_newborn = length(rlm.fit$residual)
    
	set(results_homocysteine_newborn, i, 2L, beta_homocysteine_newborn)
	set(results_homocysteine_newborn, i, 3L, SE_homocysteine_newborn)
	set(results_homocysteine_newborn, i, 4L, PVAL_homocysteine_newborn)
	set(results_homocysteine_newborn, i, 5L, N_homocysteine_newborn)
  }, error = function(err) {
	message("Error in ", colnames(beta_matrix_homocysteine_newborn)[i])
	set(results_homocysteine_newborn, i, 2L, NA)
	set(results_homocysteine_newborn, i, 3L, NA)
	set(results_homocysteine_newborn, i, 4L, NA)
	set(results_homocysteine_newborn, i, 5L, NA)
  })
 
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_homocysteine_newborn, "%\n")}
}

cat("Sorting and saving result file for model 7 (newborn homocysteine model)...\n\n")
results_homocysteine_newborn = results_homocysteine_newborn[order(results_homocysteine_newborn$p_val),]
cat("Saving the EWAS results for model 7 (newborn homocysteine model)...\n")

## saves results in .csv file
write.csv(results_homocysteine_newborn, file=paste0(cohort,".",newborn,".",homocysteine.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data:\n\n")
str(results_homocysteine_newborn)
cat("\n")
summary(results_homocysteine_newborn) # Some summary data to check everything looks good

## Calculate the epigenetic inflation factor lambda and save it in a txt file:
lambda_homocysteine_newborn = median(qchisq(results_homocysteine_newborn$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for model 7 (homocysteine newborn model) is: ", lambda_homocysteine_newborn, "\n")
write.table(lambda_homocysteine_newborn, file=paste0(cohort,".",newborn,".",homocysteine.model,".lambda.",analysis.date,".txt"), quote=FALSE, row.names=F)

## Summarise pheno data and save summary as .csv file
homocysteine_newborn_B12.tableone <- as.data.frame(print(CreateTableOne(data=dataB12_homocysteine_newborn[,-1],factorVars=c("mat.ses", "mat.smoking", "parity", "sex", "batch")),stringsAsFactors=FALSE))
write.csv(homocysteine_newborn_B12.tableone,file=paste0(cohort,".",newborn,".",homocysteine.model,".summary.",analysis.date,".csv"))

print("Model 7: newborn homocysteine b12 analysis completed") ## This completes the analysis of model 7. If you are able to run model 8, please do not delete any code from here.



############################################################## MODEL 8: NEWBORN ACTIVE MODEL #######################################################

### IF NOT AVAILABLE IN YOUR COHORT, YOU CAN DELETE THIS PART OF THE SCRIPT (UNTIL ##### COMPLETE ANALYSES ####)


### newborn active B12 model: remove cases with NA for any of newborn active B12 and main covariates
intersecting.samples <- intersect(dataB12$sample.id,colnames(betaFINAL))
dataB12_active_newborn <- na.omit(dataB12[which(dataB12$sample.id %in% intersecting.samples),unique(c(sample.id, "newborn.active.b12", batch, cell.names, main.covariates))])
k_active_newborn = which(is.na(betaFINAL), arr.ind=TRUE)
betaFINAL[k_active_newborn] = rowMedians(betaFINAL, na.rm=TRUE)[k_active_newborn[,1]]

summary(dataB12_active_newborn)## Summarizes data for newborn active b12 model
dim(dataB12_active_newborn) 


## newborn active b12 model 8: Order the files according to ID
index<-which(colnames(betaFINAL) %in% dataB12_active_newborn$sample.id)
length(index)
beta_active_newborn<-betaFINAL[,index]
beta_active_newborn<-beta_active_newborn[,order(colnames(beta_active_newborn))]
dataB12_active_newborn<-dataB12_active_newborn[order(dataB12_active_newborn$sample.id),]
ncol(beta_active_newborn)
nrow(beta_active_newborn)
all.equal(as.character(dataB12_active_newborn$sample.id),as.character(colnames(beta_active_newborn))) #Return must be TRUE, not FALSE!

betaFINAL_active_newborn<-beta_active_newborn


## Transpose the trimmed methylation data (rows are newborns and columns are CpGs)
beta_matrix_active_newborn<-t(betaFINAL_active_newborn)
dim(beta_matrix_active_newborn)
all.equal(as.character(dataB12_active_newborn$sample.id),rownames(beta_matrix_active_newborn))  #Return must be TRUE, not FALSE!


#Put the IDs in phenotype file and in methylation file into the same order
dataB12_active_newborn <- dataB12_active_newborn[dataB12_active_newborn$sample.id %in% intersect(dataB12_active_newborn$sample.id, rownames(beta_matrix_active_newborn)),]
dataB12_active_newborn <- dataB12_active_newborn[match(rownames(beta_matrix_active_newborn), dataB12_active_newborn$sample.id),]
sum(dataB12_active_newborn$sample.id == rownames(beta_matrix_active_newborn))/nrow(beta_matrix_active_newborn) #check that the IDs are in the same order, should return 1
nprobes_active_newborn = ncol(beta_matrix_active_newborn)


### run EWA for MODEL 8: NEWBORN ACTIVE MODEL
results_active_newborn = data.frame(probeID_active_newborn=colnames(beta_matrix_active_newborn),
                 	beta=rep(0, times=nprobes_active_newborn),
                 	se=rep(0, times=nprobes_active_newborn),
                 	p_val=rep(0, times=nprobes_active_newborn),
                 	n=rep(0, times=nprobes_active_newborn))

for(i in 1:nprobes_active_newborn){
  tryCatch({
	CpG_active_newborn = beta_matrix_active_newborn[,i]
	rlm.fit = rlm(CpG_active_newborn ~ newborn.active.b12 +
                	batch +
                	bcell +
                	mono + 
                  cd4t +
                	cd8t + 
                  gran + 
                  nk + 
                  nRBC +
                  mat.age + 
                  mat.ses +
                  mat.bmi +
                  mat.smoking +
                  parity +
                  sex,
              	data=dataB12_active_newborn,
              	maxit=200) 
	test = coeftest(rlm.fit, vcovHC(rlm.fit, type="HC0"))
    
	BETA_active_newborn = test[2,"Estimate"]
	SE_active_newborn = test[2,"Std. Error"]
	PVAL_active_newborn = test[2,"Pr(>|z|)"]
	N_active_newborn = length(rlm.fit$residual)
    
	set(results_active_newborn, i, 2L, BETA_active_newborn)
	set(results_active_newborn, i, 3L, SE_active_newborn)
	set(results_active_newborn, i, 4L, PVAL_active_newborn)
	set(results_active_newborn, i, 5L, N_active_newborn)
  }, error = function(err) {
	message("Error in ", colnames(beta_matrix_active_newborn)[i])
	set(results_active_newborn, i, 2L, NA)
	set(results_active_newborn, i, 3L, NA)
	set(results_active_newborn, i, 4L, NA)
	set(results_active_newborn, i, 5L, NA)
  })
 
  if(i%%4000 == 0) {cat("Progress:", 100*i/nprobes_active_newborn, "%\n")}
}


cat("Sorting and saving result file for newborn active b12 model 8 ...\n\n")
results_active_newborn = results_active_newborn[order(results_active_newborn$p_val),] 
cat("Saving the EWAS results for newborn active b12 model 8...\n")

## saves results in .csv file
write.csv(results_active_newborn, file=paste0(cohort,".",newborn,".",active.model,".results.",analysis.date,".csv"), quote=FALSE, row.names=F)

cat("Done\n\n")
cat("Checking and summary data for newborn active b12 model 8:\n\n")
str(results_active_newborn)
cat("\n")
summary(results_active_newborn) #some summary data to check everything looks good

#Calculate the epigenetic inflation factor lambda for newborn active b12 model 8:
lambda_active_newborn = median(qchisq(results_active_newborn$p_val, df=1, lower.tail = F),na.rm=T)/qchisq(0.5,1)              	 
cat("\n\nLambda for newborn active b12 model 8 is: ", lambda_active_newborn, "\n")
write.table(lambda_active_newborn, file=paste0(cohort,".",newborn,".",active.model,".lambda.",analysis.date,".txt"), col.names=T, row.names=F, quote=F)


# Summarise pheno data and save summary as .csv file
active_newborn_B12.tableone <- as.data.frame(print(CreateTableOne(data=dataB12_active_newborn[,-1],factorVars=c("mat.ses", "mat.smoking", "parity", "sex", "batch")),stringsAsFactors=FALSE))
write.csv(active_newborn_B12.tableone,file=paste0(cohort,".",newborn,".",active.model,".summary.",analysis.date,".csv"))

print("Model 8: newborn active B12 analysis completed") ## This completes the analysis of model 8. 



####################################################################### COMPLETE ANALYSES #############################################################

print("All analyses completed")

