source("new iterative functions.R")

############################################################################################
iterative12 <- 1   ###########################################################################
############################################################################################
# 1 iterative
# 2 non-iterative

indorglobal <- 1
# 1 individual optimization used
# 2 global optimization used

parSetting <- 2
# 1 both fixed + 1 free parameter
# 2 second free parameter + both free

# --------- Specify which functions are used for the new iterative optimization ---------
##############################################################
funcType <- 3 ###############################################
##############################################################
# 1 iterative independent of trial order (half evidence, half prior rate)
# 2 iterative independent of trial order (1-prior rate)
# 3 iterative dependent on trial order 

library("priorinference")

# Data file from Ella
x9data = read.csv(
  "ella_total_allDataCleaned.csv",
  header = TRUE,
  na.strings = c("", " ", "NA")
)

x9data <- x9data[, colnames(x9data)!="X"]

#arguments for the new function simplePragmaticSpeakerWithPrefPriorAll_dependentOnOrder
weights <- c(0.3, 0.4, 0.5, 0.6)
trial <- c(0, 1, 2, 3)

# adding feature property codes (which feature was uttereed, which features were questioned)
uttFeat <- ifelse(x9data$utterance=="green" | x9data$utterance=="red" | x9data$utterance=="blue", 3,
                  ifelse(x9data$utterance=="solid" | x9data$utterance=="striped" | x9data$utterance=="polka-dotted", 2, 1))
x9data$uttFeat <- uttFeat
targetFeat <- x9data$targetFeatureNum

## adding the 1-27 target and object1, object2 & object3 code.
temp <- x9data$simulatedAnswer
temp2 <- (temp - temp %% 10) / 10
temp3 <- (temp2 - temp2 %% 10) / 10
targetOC27 <- temp3 + 3 * ((temp2 %% 10) - 1) + 9 * ((temp %% 10) - 1)
x9data$targetOC27 <- targetOC27

temp <- x9data$obj1
temp2 <- (temp - temp %% 10) / 10
temp3 <- (temp2 - temp2 %% 10) / 10
obj1OC27 <- temp3 + 3 * ((temp2 %% 10) - 1) + 9 * ((temp %% 10) - 1)
x9data$obj1OC27 <- obj1OC27

temp <- x9data$obj2
temp2 <- (temp - temp %% 10) / 10
temp3 <- (temp2 - temp2 %% 10) / 10
obj2OC27 <- temp3 + 3 * ((temp2 %% 10) - 1) + 9 * ((temp %% 10) - 1)
x9data$obj2OC27 <- obj2OC27

temp <- x9data$obj3
temp2 <- (temp - temp %% 10) / 10
temp3 <- (temp2 - temp2 %% 10) / 10
obj3OC27 <- temp3 + 3 * ((temp2 %% 10) - 1) + 9 * ((temp %% 10) - 1)
x9data$obj3OC27 <- obj3OC27

## now determining the recorded subject responses
subjectResponses <- matrix(0,length(x9data$workerid),3)

for(i in c(1:length(x9data$workerid))) {
  subjectResponses[i,1] <- x9data$normResponse0[i] + 1e-100
  subjectResponses[i,2] <- x9data$normResponse1[i] + 1e-100
  subjectResponses[i,3] <- x9data$normResponse2[i] + 1e-100
  #  subjectResponses[i,1:3] <- subjectResponses[i,1:3] / sum(subjectResponses[i,1:3]) # Ella already normalized the data
}

## ordering the recorded subject responses such that they can be compared directly
#   to the model predictions
##             (particularly for visual comparison in the table)
subjectResponsesOrdered <- matrix(NA ,length(x9data$workerid),9)
for(i in c(1:length(x9data$workerid))) {
  for(j in 1:3) {
    subjectResponsesOrdered[i, (j+(targetFeat[i]-1)*3)] <- subjectResponses[i,j]
  }
}
subjectResponsesOrdered <- round(subjectResponsesOrdered, digits=5)

## Reordering objects in input data

targetObject <- rep(NA, length(x9data$workerid))
object2 <- rep(NA, length(x9data$workerid))
object3 <- rep(NA, length(x9data$workerid))

for (i in 1:length(x9data$workerid)){
  if(targetOC27[i] == obj1OC27[i]){
    targetObject[i] <- targetOC27[i]
    object2[i] <- obj2OC27[i]
    object3[i] <- obj3OC27[i]
  } else if (targetOC27[i] == obj2OC27[i])
  {targetObject[i] <- obj2OC27[i]
  object2[i] <- obj1OC27[i]
  object3[i] <- obj3OC27[i]
  } else {
    targetObject[i] <- obj3OC27[i]
    object2[i] <- obj1OC27[i]
    object3[i] <- obj2OC27[i]
  }
}

## recording KL divergence and parameters (base model, 1 param, 2 params)
workerIDs <- x9data$workerid
idMax <- max(workerIDs)
llWorkers12 <- matrix(0,length(unique(workerIDs)), 8)
paramsWorkers12 <- matrix(0,length(unique(workerIDs)), 10)

##########

## recording KL divergence and parameters (base model, 1 param, 2 params)
workerIDs <- x9data$workerid
idMax <- max(workerIDs)

#################################################

if(iterative12 == 1) { #iterative
  if(indorglobal == 1){ #individual 
    if(funcType == 1){
    paramsWorkers12 <- as.matrix(read.csv("optimized/x9params_simpleRSA_indOpt_iterative_indep_half.csv"))
    llWorkers12 <- as.matrix(read.csv("optimized/x9KLDivs_simpleRSA_indOpt_iterative_indep_half.csv"))
    print("used funcType == 1")
    
    } else if (funcType == 2) {
      paramsWorkers12 <- as.matrix(read.csv("optimized/x9params_simpleRSA_indOpt_iterative_indep_pr.csv"))
      llWorkers12 <- as.matrix(read.csv("optimized/x9KLDivs_simpleRSA_indOpt_iterative_indep_pr.csv"))
      print("used funcType == 2")
      
    } else {
      paramsWorkers12 <- as.matrix(read.csv("optimized/x9params_simpleRSA_indOpt_iterative_dep.csv"))
      llWorkers12 <- as.matrix(read.csv("optimized/x9KLDivs_simpleRSA_indOpt_iterative_dep.csv"))
      print("used funcType == 3")
    }
    
  } else { #iterative global
    paramsWorkers12 <- as.matrix(read.csv("x9Params_simpleRSA_globalOpt_iterative.csv"))
    llWorkers12 <- as.matrix(read.csv("x9KLDivs_simpleRSA_globalOpt_iterative.csv"))
  }
  
# not relevant for iterative scenario
} else { # non iterative
  if(indorglobal == 1){
  paramsWorkers12 <- as.matrix(read.csv("x9Params_simpleRSA_indOpt_nonIterative.csv"))
  llWorkers12 <- as.matrix(read.csv("x9KLDivs_simpleRSA_indOpt_nonIterative.csv"))
  }
  else {
    paramsWorkers12 <- as.matrix(read.csv("x9Params_simpleRSA_globalOpt_noniterative.csv"))
    llWorkers12 <- as.matrix(read.csv("x9KLDivs_simpleRSA_globalOpt_noniterative.csv"))
  }
}


llWorkers12 <- llWorkers12[,c(2:ncol(llWorkers12))]
paramsWorkers12 <- paramsWorkers12[,c(2:ncol(paramsWorkers12))]

###########################################
constellationCode <- matrix(0,length(x9data$workerid),6)
uniqueCCode <- rep(0, length(x9data$workerid))
postListMat1Opt <- matrix(0,length(x9data$workerid),9)
postListMat2Opt <- matrix(0,length(x9data$workerid),9)

# added:
postListMat3Opt <- matrix(0,length(x9data$workerid),9)
postListMat45Opt <- matrix(0,length(x9data$workerid),9)

postListMat4Opt <- matrix(0,length(x9data$workerid),9)
postListMat56Opt <- matrix(0,length(x9data$workerid),9)
postListMat78Opt <- matrix(0,length(x9data$workerid),9)

postListMat910Opt <- matrix(0,length(x9data$workerid),9)

###########################################

# when global optimization was used
if(indorglobal == 2){
  params1 <- paramsWorkers12[1,2]
  params2 <- paramsWorkers12[1,3]
  params12 <- paramsWorkers12[1,c(4,5)]
}

logLik <- rep(0,length(x9data$workerid))
workerID <- -1

#weights <- c(0.3, 0.4, 0.5, 0.6)
#trial <- c(0, 1, 2, 3)
#switched x9data$X with x9data$workerid, because we removed X
for(i in c(1:length(x9data$workerid))) { 
  objectConstellation <- c(targetObject[i],object2[i],object3[i])
  print("objectConstellation")
  print(objectConstellation)
  featChoice <- uttFeat[i]
  constellationCode[i,] <- getConstellationCode(objectConstellation, featChoice)[[1]]
  print("constellationCode[i,]")
  print(constellationCode[i,])
  uc <- 0
  for(j in c(1:6)) {
    uc <- (uc * 10) + constellationCode[i,j]
  }
  uniqueCCode[i] <- uc
  if(indorglobal == 1){
  if(workerID != x9data$workerid[i]) {
    workerID <- x9data$workerid[i]
    # all three optimized    params <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(11:13)]
    if(funcType == 1){
    # ---- one parameter optimized ----
    #params1 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(2)] # softness
    #params2 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(3)] # obedience
    #params3 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(4)] # prior rate
    
    # ---- two parameter optimized ----
    #params12 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(4:5)] # softness and obedience
    
    } else if (funcType == 2) {
      #funcType == 3
      
      # ---- one parameter optimized ----
      params1 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(2)] # softness
      params2 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(3)] # obedience
      params3 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(4)] # prior rate
      
      # ---- two parameter optimized ----
      params13 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(5:6)] # softness and prior rate
      params13_1 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(7:8)] # softness and prior rate
      params23 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(9:10)] # obedience and prior rate
      
      # print(params)
    } else {
      #funcType == 3
      
      # ---- one parameter optimized ----
      params1 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(2)] # softness, obed = 0
      params1_1 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(3)] # softness, obed = 0.1
      params2 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(4)] # obedience
      
      # ---- two parameter optimized ----
      params12 <- paramsWorkers12[which(paramsWorkers12[,1]==workerID)[1],c(5:6)] # softness and obedience
      
    } 
      
  }
  }
  validUtterances <- determineValidUtterances(objectConstellation)

  if( (i-1)%%4 == 0) {
    priorPrefAll_1 <- getPreferencesPrior(x9data[i,"targetFeatureNum"])
    priorPrefAll_2 <- getPreferencesPrior(x9data[i,"targetFeatureNum"])
    priorPrefAll_56 <- getPreferencesPrior(x9data[i,"targetFeatureNum"])
    priorPrefAll_78 <- getPreferencesPrior(x9data[i,"targetFeatureNum"])
    priorPrefAll_910 <- getPreferencesPrior(x9data[i,"targetFeatureNum"])
    
  } # uniform focussing on the feature type in question.

  if(iterative12 == 1) { #iterative
    if(parSetting == 1) { # one parameter optimized, two other fixed
      if (funcType == 1){
        # iterative function independent of trial order half, half
        
        # default parameters
        postListMat1Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep(objectConstellation, featChoice,
                                                                                                      0, 0, priorPrefAll_1)
        # 1st parameter optimized: softness                                                                                                                                                                                      0, 0, priorPrefAll_1, 0.5)
        postListMat2Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep(objectConstellation, featChoice,
                                                                                                      abs(params1[1]), 0, priorPrefAll_2)
        # 2nd parameter optimized: obedience                                                                                                                                                                            abs(params1[1]), 0, priorPrefAll_2, 0.5)
        postListMat3Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep(objectConstellation, featChoice,
                                                                                                         0 , abs(params2[1]), priorPrefAll_2)
                                                                                                     
      } else if(funcType == 2){
        # iterative function independent of trial order  (1 - prior rate)
        # default parameters
        postListMat1Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep_pr(objectConstellation, featChoice,
                                                                                                         0, 0, priorPrefAll_1, 0.5)
        # 1st parameter optimized: softness  
        postListMat2Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep_pr(objectConstellation, featChoice,
                                                                                                         abs(params1[1]), 0, priorPrefAll_2, 0.5)
        # 2nd parameter optimized: obedience 
        postListMat3Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep_pr(objectConstellation, featChoice,
                                                                                                         0 , abs(params2[1]), priorPrefAll_2, 0.5)
        # 3rd parameter optimized: prior rate
        postListMat4Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep_pr(objectConstellation, featChoice,
                                                                                                         0 , 0, priorPrefAll_2, abs(params3[1]))
        
      } else {
        # funcType == 3
        postListMat1Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_dep(objectConstellation, featChoice,
                                                                                                         0, 0, priorPrefAll_1)
        # 1st parameter optimized: softness , obed = 0
        postListMat2Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_dep(objectConstellation, featChoice,
                                                                                                         abs(params1[1]), 0, priorPrefAll_2)
        # 1st parameter optimized: softness, obed = 0.1
        postListMat3Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_dep(objectConstellation, featChoice,
                                                                                                    abs(params1_1[1]), 0.1, priorPrefAll_2)
        # 3rd parameter optimized: obedience
        postListMat4Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_dep(objectConstellation, featChoice,
                                                                                                    0 , abs(params2[1]), priorPrefAll_2)
      }
      
    } else if(parSetting == 2) {
      # 2 parameter optimization
      if (funcType == 1){
        # default parameters
        postListMat1Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep(objectConstellation, featChoice,
                                                                                                      0, 0, priorPrefAll_1)
        postListMat45Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep(objectConstellation, featChoice,
                                                                                                       abs(params12[1]), abs(params12[2]), priorPrefAll_2)
                                                                        
      } else if (funcType == 2){
        
        # default parameters
        postListMat1Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep_pr(objectConstellation, featChoice,
                                                                                                      0, 0, priorPrefAll_1, 0.5)
        postListMat56Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep_pr(objectConstellation, featChoice,
                                                                                                       abs(params13[1]), 0, priorPrefAll_56, abs(params13[2]))
        postListMat78Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep_pr(objectConstellation, featChoice,
                                                                                                       abs(params13_1[1]), 0.1, priorPrefAll_78, abs(params13_1[2]))
        postListMat910Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_indep_pr(objectConstellation, featChoice,
                                                                                                        0,  abs(params23[1]), priorPrefAll_910, abs(params23[2]))
        
        
      } else {
        #funcType == 3
        # default parameters
        postListMat1Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_dep(objectConstellation, featChoice,
                                                                                                      0, 0, priorPrefAll_1)
        postListMat56Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions_dep(objectConstellation, featChoice,
                                                                                                       abs(params12[1]), abs(params12[2]), priorPrefAll_1)
                                                                                           
      }
      
      }
    priorPrefAll_1 <- postListMat1Opt[i,]
    priorPrefAll_2 <- postListMat2Opt[i,]
    priorPrefAll_3 <- postListMat3Opt[i,]
    #priorPrefAll_4 <- postListMat4Opt[i,]
    
    #priorPrefAll_5_6 <- postListMat56Opt[i,]
    #priorPrefAll_7_8 <- postListMat78Opt[i,]
    #priorPrefAll_9_10 <- postListMat910Opt[i,]
 
  # not relevant for iterative scenario
  }else if(iterative12 == 2) {
    if(parSetting == 1) {
      postListMat1Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions(objectConstellation, featChoice,
                                                                                 0, 0, priorPrefAll_1, weights, trial)
      postListMat2Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions(objectConstellation, featChoice,
                                                                                 abs(params1[1]), 0, priorPrefAll_2, weights, trial)
    }else if(parSetting == 2) {
      postListMat1Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions(objectConstellation, featChoice,
                                                                                 0,  abs(params2[1]), priorPrefAll_1, weights, trial)
      postListMat2Opt[i,] <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_NewFunctions(objectConstellation, featChoice,
                                                                                 abs(params12[1]), abs(params12[2]), priorPrefAll_2, weights, trial)
    }
  }
}

###########
## adding all those values to the x9data table.
subjectResponsesOrdered <- round(subjectResponsesOrdered, digits=5)
colnames(subjectResponsesOrdered) <- colnames(subjectResponsesOrdered, do.NULL = FALSE, prefix = "DataPost_")
x9data <- data.frame(x9data, as.data.frame(subjectResponsesOrdered))

# default parameters
postListMat1Opt <- round(postListMat1Opt, digits=5)
colnames(postListMat1Opt) <- colnames(postListMat1Opt, do.NULL = FALSE, prefix = "Post1_")
consCodeAndPosteriors <- data.frame(as.data.frame(postListMat1Opt))
x9data <- data.frame(x9data, consCodeAndPosteriors)

# ------ one parameter optimized ------
# softness
postListMat2Opt <- round(postListMat2Opt, digits=5)
colnames(postListMat2Opt) <- colnames(postListMat2Opt, do.NULL = FALSE, prefix = "Post2_")
consCodeAndPosteriorsNO <- data.frame(as.data.frame(postListMat2Opt))
x9data <- data.frame(x9data, consCodeAndPosteriorsNO)

x9data$CCode <- uniqueCCode
x9data$logLik <- logLik

#obedience
postListMat3Opt <- round(postListMat3Opt, digits=5)
colnames(postListMat3Opt) <- colnames(postListMat3Opt, do.NULL = FALSE, prefix = "Post3_")
consCodeAndPosteriorsNO <- data.frame(as.data.frame(postListMat3Opt))
x9data <- data.frame(x9data, consCodeAndPosteriorsNO)

x9data$CCode <- uniqueCCode
x9data$logLik <- logLik

# if func 2 used: prior rate
# if func 3 used: obedience
postListMat4Opt <- round(postListMat4Opt, digits=5)
colnames(postListMat4Opt) <- colnames(postListMat4Opt, do.NULL = FALSE, prefix = "Post4_")
consCodeAndPosteriorsNO <- data.frame(as.data.frame(postListMat4Opt))
x9data <- data.frame(x9data, consCodeAndPosteriorsNO)
x9data$CCode <- uniqueCCode
x9data$logLik <- logLik

# ------- two parameters optimized ------
# ---- if funcType == 1: -----
# softness and obedience (pr = 0.5)
#postListMat45Opt <- round(postListMat45Opt, digits=5)
#colnames(postListMat45Opt) <- colnames(postListMat45Opt, do.NULL = FALSE, prefix = "Post4_")
#consCodeAndPosteriorsNO <- data.frame(as.data.frame(postListMat45Opt))
#x9data <- data.frame(x9data, consCodeAndPosteriorsNO)

#x9data$CCode <- uniqueCCode
#x9data$logLik <- logLik

# ----- if funcType == 2: -----
# softness and prior rate (obed = 0)
# or
# ----- if funcType ==  3: -----
# softness and obedience
postListMat56Opt <- round(postListMat56Opt, digits=5)
colnames(postListMat56Opt) <- colnames(postListMat56Opt, do.NULL = FALSE, prefix = "Post5_")
consCodeAndPosteriorsNO <- data.frame(as.data.frame(postListMat56Opt))
x9data <- data.frame(x9data, consCodeAndPosteriorsNO)

x9data$CCode <- uniqueCCode
x9data$logLik <- logLik

# --- funcType == 2: -----
# softness and prior rate (obed = 0.1)
postListMat78Opt <- round(postListMat78Opt, digits=5)
colnames(postListMat78Opt) <- colnames(postListMat78Opt, do.NULL = FALSE, prefix = "Post6_")
consCodeAndPosteriorsNO <- data.frame(as.data.frame(postListMat78Opt))
x9data <- data.frame(x9data, consCodeAndPosteriorsNO)

x9data$CCode <- uniqueCCode
x9data$logLik <- logLik

# obedience and prior rate (soft = 0)
postListMat910Opt <- round(postListMat910Opt, digits=5)
colnames(postListMat910Opt) <- colnames(postListMat910Opt, do.NULL = FALSE, prefix = "Post7_")
consCodeAndPosteriorsNO <- data.frame(as.data.frame(postListMat910Opt))
x9data <- data.frame(x9data, consCodeAndPosteriorsNO)

x9data$CCode <- uniqueCCode
x9data$logLik <- logLik


if(iterative12 == 1) { #iterative
  if(indorglobal == 1){ #individual
  if(parSetting == 1) { #one parameter optimized
  
    if(funcType == 1){
      
      # optimized softness
      #write.csv(x9data, "data/x9dataAugm_SRSAindOpt1_fixed0_0.5_iterative_indep_half.csv")
      
      # optimized obedience
      #write.csv(x9data, "data/x9dataAugm_SRSAindOpt2_fixed0_0.5_iterative_indep_half.csv")
      
    } else if(funcType == 2){
      # optimized softness
      #write.csv(x9data, "data/x9dataAugm_SRSAindOpt1_fixed0_0.5_iterative_indep_pr.csv")
      
      # optimized obedience
      write.csv(x9data, "data/x9dataAugm_SRSAindOpt2_fixed0_0.5_iterative_indep_pr.csv")
     
      # optimized prior rate
      #write.csv(x9data, "data/x9dataAugm_SRSAindOpt3_fixed0_0.5_iterative_indep_pr.csv")
    
      } else {
      # funcType == 3
      # optimized softness, obed = 0
       #write.csv(x9data, "data/x9dataAugm_SRSAindOpt1_fixed0_iterative_dep.csv")
      
      # optimized softness, obed = 0.1
      # write.csv(x9data, "data/x9dataAugm_SRSAindOpt1_fixed0.1_iterative_dep.csv")
      
      # optimized obedience, soft = 0
       #write.csv(x9data, "data/x9dataAugm_SRSAindOpt2_fixed0_iterative_dep.csv")
    }
    
  }else if(parSetting == 2) {
    
    if (funcType == 1){
      # two parameter optim
      
      #write.csv(x9data, "data/x9dataAugm_SRSAindOpt13_fixed_notObey0_iterative_indep_half.csv")
      #write.csv(x9data, "data/x9dataAugm_SRSAindOpt13_fixed_notObey0.1_iterative_indep_half.csv")
      #write.csv(x9data, "data/x9dataAugm_SRSAindOpt12_fixed_pr0.5_iterative_indep_half.csv")
    
      } else if (funcType == 2){
      
      #write.csv(x9data, "data/x9dataAugm_SRSAindOpt13_fixed_notObey0_iterative_indep_pr.csv")
      #write.csv(x9data, "data/x9dataAugm_SRSAindOpt13_fixed_notObey0.1_iterative_indep_pr.csv")
      #write.csv(x9data, "data/x9dataAugm_SRSAindOpt23_fixed_pref0_iterative_indep_pr.csv")
      
    }else {
      # funcType == 3
      write.csv(x9data, "data/x9dataAugm_SRSAindOpt12_iterative_dep.csv")
      
      
    }
    
  }
  # global
  } else {
    if(parSetting == 1) {
      write.csv(x9data, "data/x9dataAugm_SRSAglobalOpt_fixed00_andOpt1fixed0_iterative_withNewFunctions.csv")
    }else if(parSetting == 2) {
      write.csv(x9data, "data/x9dataAugm_SRSAglobalOpt_fixed0Opt2_andOpt12_iterative_withNewFunctions.csv")
    }
  }
  
# non-iterative
} else if(iterative12 == 2) {
  if(indorglobal == 1){
  if(parSetting == 1) {
    write.csv(x9data, "data/x9dataAugm_SRSAindOpt_fixed00_andOpt1fixed0_nonIterative_withNewFunctions.csv")
  }else if(parSetting == 2) {
    write.csv(x9data, "data/x9dataAugm_SRSAindOpt_fixed0Opt2_andOpt1_nonIterative_withNewFunctions.csv")
  }
  } else {
    if(parSetting == 1) {
      write.csv(x9data, "data/x9dataAugm_SRSAglobalOpt_fixed00_andOpt1fixed0_noniterative_withNewFunctions.csv")
    }else if(parSetting == 2) {
      write.csv(x9data, "data/x9dataAugm_SRSAglobalOpt_fixed0Opt2_andOpt12_noniterative_withNewFunctions.csv")
    }
  }
}

