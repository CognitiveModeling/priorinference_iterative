source("new iterative functions.R")
library("priorinference") # load package

# ---------------------- settings ----------------------------------------
# --------- Which scenario do you want to use for optimization ? ---------
##############################################################
procType <- 3  ###############################################
##############################################################
# 1 iterative
# 2 non-iterative
# 3 iterative new implementation 05.08.21

# if you are using procType == 3:
# --------- Specify which functions are used for the new iterative optimization ---------
##############################################################
funcType <- 3 ###############################################
##############################################################
# 1 iterative independent of trial order (posterior = 0.5 * evidence +  0.5 prior) (NOT INCLUDED IN THE PACKAGE: as it is the same as funcType == 2 with prior rate fixed to 0.5)
# 2 iterative independent of trial order (posterior = (1-prior rate) * evidence + prior rate * prior)
# 3 iterative dependent on trial order


# Data file from Ella
x9data = read.csv(
  "data/ella_total_allDataCleaned.csv",
  header = TRUE,
  na.strings = c("", " ", "NA")
)

# adding feature property codes (which feature was uttered, which features were questioned)
uttFeat <- ifelse(x9data$utterance=="green" | x9data$utterance=="red" | x9data$utterance=="blue", 3,
                  ifelse(x9data$utterance=="solid" | x9data$utterance=="striped" | x9data$utterance==" polka-dotted", 2, 1))
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
subjectResponses <- matrix(0,length(x9data$X),3)

for(i in c(1:length(x9data$X))) {
  subjectResponses[i,1] <- x9data$normResponse0[i] + 1e-100
  subjectResponses[i,2] <- x9data$normResponse1[i] + 1e-100
  subjectResponses[i,3] <- x9data$normResponse2[i] + 1e-100
  #  subjectResponses[i,1:3] <- subjectResponses[i,1:3] / sum(subjectResponses[i,1:3]) # Ella already normalized the data
}

## ordering the recorded subject responses such that they can be compared directly
#   to the model predictions
##             (particularly for visual comparison in the table)
subjectResponsesOrdered <- matrix(NA ,length(x9data$X),9)
for(i in c(1:length(x9data$X))) {
  for(j in 1:3) {
    subjectResponsesOrdered[i, (j+(targetFeat[i]-1)*3)] <- subjectResponses[i,j]
  }
}
subjectResponsesOrdered <- round(subjectResponsesOrdered, digits=5)

## Reordering objects in input data

targetObject <- rep(NA, length(x9data$X))
object2 <- rep(NA, length(x9data$X))
object3 <- rep(NA, length(x9data$X))

for (i in 1:length(x9data$X)){
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
## Starting with simple base model determination:
##
workerIndex <- 1
for(workerID in c(0:idMax)) {
  allIndices <- which(workerIDs == workerID)
  if(length(allIndices)>0) {
    llWorkers12[workerIndex,1] <- workerID
    paramsWorkers12[workerIndex,1] <- workerID
    ## based model -> no change in preferences!
    llWorkers12[workerIndex,2] <- 0 # -2 * length(allIndices) * log(1/3)
    for(i in c(1:length(allIndices))) {
      for(j in c(1:3)) {
        llWorkers12[workerIndex, 2] <- llWorkers12[workerIndex, 2] +
          subjectResponses[allIndices[i],j] *
          (log(subjectResponses[allIndices[i],j]) - log(1/3))
      }
    }
    ## done with this worker -> proceed
    workerIndex <- workerIndex + 1
  }
}

weights <- c(0.3, 0.4, 0.5, 0.6)
trial <- c(0, 1, 2, 3)
#######################
## Optimizing the log likelihoods (maximum log likelihoods for each worker...)
##    1 parameter RSA model optimizations...
workerIndex <- 1
for(workerID in c(0:idMax)) {
  idICases <- which(workerIDs == workerID)
  if(length(idICases)>0) {
    ## generating data matrix for the purpose of optimization
    dataWorker <- matrix(0, length(idICases), 8)
    dataWorker[,1] <- targetObject[idICases]
    dataWorker[,2] <- object2[idICases]
    dataWorker[,3] <- object3[idICases]
    dataWorker[,4] <- uttFeat[idICases]
    dataWorker[,5] <- targetFeat[idICases]
    dataWorker[,6:8] <- subjectResponses[idICases,1:3]

    ######################### 1 parameter optimization ###########################

    if (procType == 1){
      # before optimization:
      optRes1 <- optimize(RSAModelLL1_1simpleRSA4TrialsIterative, c(0,1e+10), dataWorker)
      optRes2 <- optimize(RSAModelLL1_2simpleRSA4TrialsIterative, c(0,1e+10), dataWorker)

    }else if (procType == 2){
      optRes1 <- optimize(RSAModelLL1_1simpleRSA4TrialsIndependent, c(0,1e+10), dataWorker)
      optRes2 <- optimize(RSAModelLL1_2simpleRSA4TrialsIndependent, c(0,1e+10), dataWorker)

    }else{ # procType == 3

      if (funcType == 1){
      # -------- Iterative functions independent of trial order (half evidence, half prior rate) -----------
      # optimizing 1st parameter in iterative model: softness
      optRes1 <- optimize(LL1_1_Iterative_indep_notObey0_pr0.5, c(0,1e+10), dataWorker)

      # optimizing 2nd parameter in iterative model: obedience
      optRes2 <- optimize(LL1_2_Iterative_indep_pref0_pr0.5, c(0,1e+10), dataWorker)

    } else if (funcType == 2){
      # ----------------------------------------------------------------------------------------------------
      # -------- Iterative functions independent of trial order (1-prior rate) -----------------------------
      # optimizing 1st parameter in iterative model: softness
      optRes1 <- optimize(LL1_1_Iterative_pr_notObey0_pr0.5, c(0,1e+10), dataWorker)

      # optimizing 2nd parameter in iterative model: obedience
      optRes2 <- optimize(LL1_2_Iterative_pr_pref0_pr0.5, c(0,1e+10), dataWorker)

      # optimizing 3rd parameter in iterative model: prior rate
      optRes3 <- optimize(LL1_3_Iterative_pr_pref0_notObey0, c(0,1), dataWorker)

    } else {
      # ----------------------------------------------------------------------------------------------------
      # -------- Iterative functions dependent on trial order ----------------------------------
      # optimizing 1st parameter in iterative model: softness (obedience fixed at 0)
      optRes1 <- optimize(LL1_1_Iterative_dep_notObey0, c(0,1e+10), dataWorker)

      # optimizing 1st parameter in iterative model: softness (obedience fixed at 0.1)
      optRes1_1 <- optimize(LL1_1_Iterative_dep_notObey0.1, c(0,1e+10), dataWorker)

      # optimizing 2nd parameter in iterative model: obedience (softness fixed at 0)
      optRes2 <- optimize(LL1_2_Iterative_dep_pref0 , c(0,1e+10), dataWorker)
    }

    }

    ######################### recording results #################################
    if (funcType == 1){
    ## 1 param RSA model
    llWorkers12[workerIndex,3] <- optRes1$objective #softness
    llWorkers12[workerIndex,4] <- optRes2$objective #obedience

    ## resulting parameter choice
    paramsWorkers12[workerIndex,2] <- optRes1$minimum #softness
    paramsWorkers12[workerIndex,3] <- optRes2$minimum #obedience

    } else if (funcType == 2){

      ## 1 param RSA model
      llWorkers12[workerIndex,3] <- optRes1$objective #softness
      llWorkers12[workerIndex,4] <- optRes2$objective #obedience
      llWorkers12[workerIndex,5] <- optRes3$objective #prior rate

      ## resulting parameter choice
      paramsWorkers12[workerIndex,2] <- optRes1$minimum #softness
      paramsWorkers12[workerIndex,3] <- optRes2$minimum #obedience
      paramsWorkers12[workerIndex,4] <- optRes3$minimum #prior rate


    }else {

      ## 1 param RSA model
      llWorkers12[workerIndex,3] <- optRes1$objective #softness (obed = 0)
      llWorkers12[workerIndex,4] <- optRes1_1$objective #softness (obed = 0.1)
      llWorkers12[workerIndex,5] <- optRes2$objective #obedience

      ## resulting parameter choice
      paramsWorkers12[workerIndex,2] <- optRes1$minimum #softness (obed = 0)
      paramsWorkers12[workerIndex,3] <- optRes1_1$minimum #softness (obed = 0.1)
      paramsWorkers12[workerIndex,4] <- optRes2$minimum #obdience

    }

    ####
    print(llWorkers12[workerIndex,])
    print(paramsWorkers12[workerIndex,])
    ####
    workerIndex <- workerIndex + 1
  }
}

######################### 2 parameter optimization ###########################

## Optimizing the log likelihoods (maximum log likelihoods for each worker...)
print("Starting optimization with two free parameters Simple RSA model... ")
workerIDs <- x9data$workerid
idMax <- max(workerIDs)

# llWorkers12 <- matrix(0,length(unique(workerIDs)), 2)
workerIndex <- 1
for(workerID in c(0:idMax)) {
  idICases <- which(workerIDs == workerID)
  if(length(idICases)>0) {
    ## generating data matrix for the purpose of optimization
    dataWorker <- matrix(0, length(idICases), 8)
    dataWorker[,1] <- targetObject[idICases]
    dataWorker[,2] <- object2[idICases]
    dataWorker[,3] <- object3[idICases]
    dataWorker[,4] <- uttFeat[idICases]
    dataWorker[,5] <- targetFeat[idICases]
    dataWorker[,6:8] <- subjectResponses[idICases,1:3]

# before optimization:
    if (procType == 1){
      # before optimization: llWorkers12[1,3] <- RSAModelLL1(c(.2), dataWorker)
      optRes12 <- optim(c(.2, .2), RSAModelLL2_simpleRSA4TrialsIterative, method="L-BFGS-B", gr=NULL, dataWorker,
                         lower = c(0,0), upper = c(1e+10,1e+10))
    } else if (procType == 2){
      optRes12 <- optim(c(.2, .2), RSAModelLL2_simpleRSA4TrialsIndependent, method="L-BFGS-B", gr=NULL, dataWorker,
                         lower = c(0,0), upper = c(1e+10,1e+10))
    } else {
    # procType == 3
    if (funcType == 1){
      # -------- Iterative functions independent of trial order (half evidence, half prior rate) -----------
      #optRes12 <- optim(c(.2, .2), LL2_12_Iterative_indep_pr0.5, method="L-BFGS-B", gr=NULL, dataWorker,
                        #lower = c(0,0), upper = c(1e+10,1e+10))

    } else if (funcType == 2){
      # -------- Iterative functions independent of trial order (1 - prior rate) ---------------------------

      # optimizing 1st and 3rd parameters in iterative model: softness + prior rate
      optRes13 <- optim(c(.2, .2), LL2_13_Iterative_pr_notObey0 , method="L-BFGS-B", gr=NULL, dataWorker,
                        lower = c(0,0), upper = c(1e+10,1))

      # optimizing 1st and 3rd parameters in iterative model: softness + prior rate
      optRes13_1 <- optim(c(.2, .2), LL2_13_Iterative_pr_notObey0.1, method="L-BFGS-B", gr=NULL, dataWorker,
                          lower = c(0,0), upper = c(1e+10,1))

      # optimizing 2nd and 3rd parameters in iterative model: obedience + prior rate
      optRes23 <- optim(c(.2, .2), LL2_23_Iterative_pr_pref0, method="L-BFGS-B", gr=NULL, dataWorker,
                        lower = c(0,0), upper = c(1e+10,1))
    # ----------------------------------------------------------------------------------------------------
    } else {
      # funcType == 3
      # ------------------ Iterative functions dependent on trial order ------------------------------------

      # optimizing 1st and 3rd parameters in iterative model: softness + prior rate
      optRes12 <- optim(c(.2, .2), LL2_12_Iterative_dep, method="L-BFGS-B", gr=NULL, dataWorker,
                        lower = c(0,0), upper = c(1e+10,1e+10))

    }
  }
    ##########
    ## 2 param RSA model2

    if(funcType == 1){
      ## max likelihood parameter choice
      # optimization: softness and obedience (prior rate 0.5)

      #llWorkers12[workerIndex,5] <- optRes12$value

      #paramsWorkers12[workerIndex,4] <- optRes12$par[1]
      #paramsWorkers12[workerIndex,5] <- optRes12$par[2]

    } else if (funcType == 2) {

    # optimization: softness and prior rate (obedience fixed = 0)
    llWorkers12[workerIndex,6] <- optRes13$value

    paramsWorkers12[workerIndex,5] <- optRes13$par[1]
    paramsWorkers12[workerIndex,6] <- optRes13$par[2]

    # optimization: softness and prior rate (obedience fixed = 0.1)
    llWorkers12[workerIndex,7] <- optRes13_1$value

    ## max likelihood parameter choice
    paramsWorkers12[workerIndex,7] <- optRes13_1$par[1]
    paramsWorkers12[workerIndex,8] <- optRes13_1$par[2]

    # optimization: obedience and prior rate
    llWorkers12[workerIndex,8] <- optRes23$value

    ## max likelihood parameter choice
    paramsWorkers12[workerIndex,9] <- optRes23$par[1]
    paramsWorkers12[workerIndex,10] <- optRes23$par[2]

    } else {
      # if funcType == 3
      # optimization: softness and prior rate (obedience fixed = 0)
      llWorkers12[workerIndex,6] <- optRes12$value

      paramsWorkers12[workerIndex,5] <- optRes12$par[1]
      paramsWorkers12[workerIndex,6] <- optRes12$par[2]
    }

    print(llWorkers12[workerIndex,])
    print(paramsWorkers12[workerIndex,])
    workerIndex <- workerIndex + 1
    }
}

## ------ writing out result tables ---------
if(procType == 1) {
write.csv(llWorkers12, "X9_Data/x9KLDivs_simpleRSA_indOpt_iterative.csv")
write.csv(paramsWorkers12, "X9_Data/x9Params_simpleRSA_indOpt_iterative.csv")
# ----------------------------------------------
} else if (procType == 2){
write.csv(llWorkers12, "X9_Data/x9KLDivs_simpleRSA_indOpt_nonIterative.csv")
write.csv(paramsWorkers12, "X9_Data/x9Params_simpleRSA_indOpt_nonIterative.csv")
# ----------------------------------------------
} else{
  # procType == 3: new iterative functions
  # iterative, independent of trial order, half evidence, half prior rate
  if (funcType == 1) {
        write.csv(llWorkers12, "optimized/x9KLDivs_simpleRSA_indOpt_iterative_indep_half.csv")
        write.csv(paramsWorkers12, "optimized/x9params_simpleRSA_indOpt_iterative_indep_half.csv")

    # ----------------------------------------------
    # iterative, independent of trial order, (1 - prior rate)
   } else if (funcType == 2){
        write.csv(llWorkers12, "optimized/x9KLDivs_simpleRSA_indOpt_iterative_indep_pr.csv")
        write.csv(paramsWorkers12, "optimized/x9params_simpleRSA_indOpt_iterative_indep_pr.csv")
    # ----------------------------------------------
    # funcType == 3: iterative, dependent on trial order
  } else {
      write.csv(llWorkers12, "optimized/x9KLDivs_simpleRSA_indOpt_iterative_dep.csv")
      write.csv(paramsWorkers12, "optimized/x9params_simpleRSA_indOpt_iterative_dep.csv")
  }
}

