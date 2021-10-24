#' Simple pragmatic speaker with all prior preferences iterative function
#' Iterative function dependent on trial order (prior rate)
#'
#' @description Simple-RSA (iterative, dependent on trial order)
#'
#' The simple pragmatic speaker considers all "imaginable" (i.e. implemented)
#' preference distributions over objects of the listener.
#'
#' Starting with a prior assumption over the possible listener's preferences. It
#' then infers the posterior over these preferences given the listener makes a
#' particular object choice. P(listener's feature value preferences | utterance,
#' object choice by the listener, prior over listener's feature value
#' preferences).
#'
#' This function takes the evidence from the current trial in consideration and also
#' the prior from the trials before and computes the posterior over the feature preferences of the listener.
#' \code{posterior = (1 - prior rate) x evidence + (prior rate) x prior}.
#'
#' @param utterance The uttered word by the speaker that the listener hears.
#'
#' An index referring to one of the values in the vector validUtterances.
#'
#' @param obj The object chosen by the listener. A value referring to the index
#'   1,2 or 3.
#'
#' @param preferencesPriorAll A vector of length 9.
#'
#' Probability mass over all feature values.
#'
#' Gives a prior preferences distribution over all (nine) feature values.
#'
#' \code{preferencesPriorAll <- rep(1/9, 9)}
#'
#' @param validUtterances A vector of utterances that correspond to all feature
#' values present in the current objects in the scene.
#'
#' For example, it only makes sense to utter \emph{"red"} in a scene if there
#' are \emph{red} objects present.
#'
#' @param currentObjects Vector of three values in \code{{1,...,27}} specifying
#' the target and the other two objects.
#'
#' The target is the first object in the vector \code{(index = 1)}.
#'
#' @param uttToObjProbs A matrix. The rows map each possible utterance that
#' corresponds to each present feature value of the current objects. The
#' columns represent the three objects in the scene.
#'
#' This reflects the obedience-parameter and which objects match the
#' respective utterance. The matrix shows the probability that a certain
#' object is chosen following a certain utterance, that is valid in the scene.
#' The number of rows of the matrix match the length of the validUtterances
#' vector.
#' @param objectPreferenceSoftPriors A list of preference priors for all valid
#' utterances based on the object in the scene.
#'
#' The list has as many rows as the length of the validUtterances vector + 1.
#'
#' Each row in the list contains a vector of length 3, as there are three
#' objects in the scene.
#'
#' The extra row is for the case of no feature preferences whatsoever, i.e.
#' uniform prior over all three objects in the scene.
#'
#' @param weights A vector of length 4 including the weight by which the prior is weighed.
#'
#' weights <- c(0.3, 0.4, 0.5, 0.6)
#'
#' @param trial A vector of length 4 including the number of the current trial.
#'
#' trial <- c(1,2,3,4)
#'
#' @return A vector of length 9. It contains the normalized probability over
#'   preferences (priors).
#'
#' @details
#' This is function is the second of two functions that are used in the iterative setting using the prior rate parameter.
#' The first one is: \code{\link{simplePragmaticSpeakerWithPrefPriorAll_indepOfOrder_pr}}.
#'
#' @examples
#' \donttest{simplePragmaticSpeakerWithPrefPriorAll_depOnOrder(utterance, obj, preferencesPriorAll,
#' validUtterances, currentObjects,
#' uttToObjProbs, objectPreferenceSoftPriors) }
#'
#' output:
#' [1] 0.5333333 0.1333333 0.3333333 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#' [9] 0.0000000
simplePragmaticSpeakerWithPrefPriorAll_depOnOrder <- function(utterance, obj, preferencesPriorAll,
                                                              validUtterances, currentObjects,
                                                              uttToObjProbs, objectPreferenceSoftPriors) {
  #cat("preferencesPriorAll", preferencesPriorAll, "\n")
  preferencesPrior <- preferencesPriorAll[validUtterances]
  prefPost <- rep(0, length(validUtterances))

  weights <- c(0.3, 0.4, 0.5, 0.6)
  trial <- c(1,2,3,4)

  for (pref in c(1:length(validUtterances))) {
    for (trial in trial){
      if (preferencesPrior[pref] > 0) {
        pp <-
          simpleListener(utterance,
                         uttToObjProbs,
                         objectPreferenceSoftPriors[[pref]])
        # trial + 1 because in the csv file we start counting at 0 (trialNum)
        prefPost[pref] <- (1-(weights[trial+1])) * pp[obj] + weights[trial+1] * preferencesPrior[pref]
        #print("weights[trial+1]")
        #print(weights[trial+1])

      }
    } # else{} # prior preference for this preference is zero
  } # done with considering all possible preferences.
  if (sum(prefPost) == 0) { # no evidence for any preferences... -> no inference
    return(preferencesPriorAll)
  }
  # normalizing relevant posterior preferences such that the sum is equal to their prior probability mass
  #   sum(preferencesPrior) is the probability mass of the full prior that we are "entitled" to redistribute
  #           because it concerns the features present in the trial
  #   prefPost / sum(prefPost) is the normalized posterior, so that the updated vector sums up to 1
  prefPost <- sum(preferencesPrior) * (prefPost * (1-weights[trial+1])) / sum(prefPost)
  # replacing the relevant old prior preferences values in preferencesPriorAll with their posteriors (which become the new priors)
  preferencesPriorAll[validUtterances] <- prefPost
  #
  return(preferencesPriorAll / sum(preferencesPriorAll)) # Sidenote: the renormalization here should not be really necessary because it is taken care of by a scaled insertion of new values in the two commands above.
}

# 3rd version (Dependent on trial order)

#' Determine speaker's inference of the posterior listener preferences (iterative setting, dependent on trial order)
#'
#' @description
#' Simple RSA (iterative, dependent on trial order)
#'
#' This function calculates the speaker's posterior guess of the
#' feature value preferences of the listener in the iterative setting.
#' That means how the speaker infers the preferences of the listener based on the object choice.
#'
#' @param currentObjects A vector of three values in \code{{1,...,27}} specifying the target and the other two objects in the scene.
#'
#' The target is the first object in the vector \code{(index = 1)}.
#' @param featureUtt One of the values \code{{1,2,3}} specifying which feature is uttered (i.e. shape = 1 / texture = 2 / or color = 3).
#'
#' @param softPrefValue A parameter value between \code{[0,infinity)} (The larger the higher the tendency towards uniform liking).
#'
#' Value reflects how categorical the listener's preferences are:
#'
#' \strong{0:} The listener always picks her preferred object.
#'
#' If the listener prefers \emph{red} objects, she will always pick the \emph{red} object in the scene.
#'
#' \strong{infinity:} It is as likely for the listener to pick \emph{green}, \emph{blue} or \emph{red} objects.
#' @param notObeyInst Determines the extent to which the instruction of the speaker is obeyed by the listener.
#'
#' (0 = full obedience, infinity = full instruction ignorance).
#'
#' \strong{Example:}
#'
#' \strong{0:} Listener always picks \emph{red} objects following the utterance \emph{"red"}.
#'
#' \strong{infinity:} Listener as likely to pick \emph{green, blue} or \emph{red} objects even if the utterance is \emph{"red"}.
#'
#' @param priorPrefAll A vector of length 9.
#'
#' Probability mass over all feature values.
#'
#' Gives a prior preferences distribution over all (nine) feature values.
#'
#' @return A vector of length 9. It contains the speaker's inference of the feature value preferences of
#' the listener dependent on the trial order.
determineSpeakerPostListPrefsSimpleRSAWithPriorPref_dep <- function(currentObjects, featureUtt, softPrefValue,
                                                                                 notObeyInst, priorPrefAll) {
  validUtterances <- determineValidUtterances(currentObjects)
  mapObjToUtt <- mapObjectToUtterances(currentObjects)
  uttToObjProbs <- determineUttToObjectProbs(validUtterances,
                                             currentObjects,
                                             mapObjToUtt, notObeyInst)
  mapUttToObjDeterministic <- determineUttToObjectProbs(validUtterances,
                                                        currentObjects,
                                                        mapObjToUtt, 0)
  #  print(uttToObjProbs)
  objectPreferenceSoftPriors <- getObjectPreferencePriors(validUtterances, currentObjects,
                                                          softPrefValue, mapUttToObjDeterministic)
  prefPostAll <- simplePragmaticSpeakerWithPrefPriorAll_depOnOrder(which(validUtterances==allObjectsToUtterancesMappings[currentObjects[1],featureUtt]),
                                                                   1, priorPrefAll, validUtterances, currentObjects,
                                                                   uttToObjProbs, objectPreferenceSoftPriors)
  return(prefPostAll)
}
# ------------------

# functions for optimization (dependent on order)

#' Simple RSA model Kullback-Leibler divergence determination
#' (iterative setting, dependent on trial order)
#'
#' @description
#' Simple RSA (iterative, dependent on trial order)
#'
#' The function calculates the optimal parameter values of the free parameters by
#' estimating the log-likelihood of the RSA model given model parameters and data.
#' It also determines the actual RSA model Kullback-Leibler divergence.
#'
#' 2 parameter optimization considering only the available feature values present in the scene, i.e. feature values of shape, texture and color.
#' This function is used in the iterative dependent on the trial scenario.
#'
#' @param data A matrix with data rows.
#'
#' column structure: \code{[1:OC1,OC2,OC3,4:numUttOptions,7-X:TurkerSliderValues]}
#'
#' \strong{1:OC1} Object 1. A value between 1 and 27.
#'
#' \strong{2:OC2} Object 2. A value between 1 and 27.
#'
#' \strong{3:OC3} Object 3. A value between 1 and 27.
#'
#' \strong{4:numUttOptions} The number of valid utterances in the scene.
#'
#' \strong{7-X:TurkerSliderValues} These columns contain the participants' slider values.
#'
#' @param par1 \describe{
#' \item{softness parameter}{The strength of "preferring one entity over others". (The larger the value the higher the tendency towards uniform liking)}
#' }
#'
#' @param par2 \describe{
#' \item{non-obedience parameter}{The extent to which the instruction of the speaker is obeyed by the listener.
#' (0 = full obedience, infinity = full instruction ignorance).}
#' }
#'
#' @return Minimized Kullback-Leibler divergence and the optimal parameter values.
#'
#' @details
#' This function is used in \code{\link{LL1_1_Iterative_dep_notObey0}},
#'
#' \code{\link{LL1_1_Iterative_dep_notObey0.1}},
#'
#' \code{\link{LL1_2_Iterative_dep_pref0}},
#'
#' \code{\link{LL2_12_Iterative_dep}}.
#'
#' @export
#' @return Minimized Kullback-Leibler divergence and the optimal parameters.
RSAModelKLDiv3params_simpleRSA4TrialsIterative_dep <- function(data, par1, par2) {
  #  print(params)
  llRes <- 0
  for(i in c(1:nrow(data))) {
    if( (i-1)%%4 == 0) {
      preferencesPriorAll <- getPreferencesPrior(data[i,5]) # focussing on the feature type in question.
    }
    ## determining the object and utterance
    currentObjects <- c(data[i,1],data[i,2],data[i,3])
    uttFeat <- data[i,4]
    ##
    validUtterances <- determineValidUtterances(currentObjects)
    ## determining the model predictions
    prefPostAll <- determineSpeakerPostListPrefsSimpleRSAWithPriorPref_dep(currentObjects, uttFeat, abs(par1), abs(par2), preferencesPriorAll)
    ## adding the KL Divergence terms of the relevant feature values for the two sets of answers.
    ##
    ## adding the negative log likelihoods
    for(j in c(1:3)) {
      llRes <- llRes + data[i, 5+j] *
        ( log(data[i, 5+j] + 1e-100) - log(prefPostAll[j + (data[i, 5]-1)*3] + 1e-100) )
    }
    preferencesPriorAll <- prefPostAll
  }
  return(llRes)
}

# ---- one parameter optimization ----
# optimizing 1st parameter in iterative model: softness

#' Cost function for one parameter optimization (iterative setting, dependent on trial order).
#' Optimizing softness.
#' Non-obedience fixed at 0.
#'
#' @description
#' Simple RSA
#'
#' 1 parameter optimization; The softness parameter is optimized.
#'
#' The non-obedience parameter is fixed.
#' @param params One value vector, which specifies one of two parameters to be optimized:
#' \enumerate{
#'   \item{softPrefValue is optimized, i.e. The strength of "preferring one entity over others". (The larger the value the higher the tendency towards uniform liking)}
#'   \item{non-obedience (default = 0), i.e. The extent to which the instruction of the speaker is obeyed by the listener.
#'   (0 = full obedience, infinity = full instruction ignorance)}
#'}
#'
#' @param data A Matrix with data rows.
#'
#' column structure:
#'
#' [1:OC1,OC2,OC3,4:UUFeat, 5:Q1Feat,6:Q2Feat]
#'
#' [7:Q1AnswerV1,V2,V3, 10:Q2AnswerV1,V2,V3]
#'
#' \strong{1:OC1} Object 1. A value between 1 and 27.
#'
#' \strong{2:OC2} Object 2. A value between 1 and 27.
#'
#' \strong{3:OC3} Object 3. A value between 1 and 27.
#'
#' \strong{4:UUFeat} Uttered feature. A number between 1 and 3. (1: shape, 2: pattern, 3: color)
#'
#' \strong{5:Q1Feat} Questioned feature 1. A number between 1 and 3. (1: shape, 2: pattern, 3: color).
#'
#' Example: If you utter "blue" (feature: color), then you can learn something about shape and texture preferences.
#'
#' \strong{6:Q2Feat} Questioned feature 2. A number between 1 and 3. (1: shape, 2: pattern, 3: color).
#'
#' Example: If you utter "blue" (feature: color), then you can learn something about shape and texture preferences.
#'
#' \strong{7:Q1AnswerV1, V2, V3} The columns 7-9 contain the participants' slider values for the first questioned feature.
#'
#' \strong{10:Q2AnswerV1, V2, C3} The columns 10-12 contain the participants' slider values for the second questioned feature.
#'
#' @return Minimized Kullback-Leibler divergence and the optimal parameter values.
#' @details
#' This function uses \code{\link{RSAModelKLDiv3params_simpleRSA4TrialsIterative_dep}}.
#' @export
LL1_1_Iterative_dep_notObey0 <- function(params, data) {
  return(RSAModelKLDiv3params_simpleRSA4TrialsIterative_dep(data, abs(params[1]), 0))
}

# ---- optimizing 1st parameter in iterative model: softness ----

#' Cost function for one parameter optimization (iterative setting, dependent on trial order).
#' Optimizing softness.
#' Non-obedience fixed at 0.1.
#'
#' @description
#' Simple RSA
#'
#' 1 parameter optimization; The softness parameter is optimized.
#'
#' The non-obedience parameter is fixed.
#' @param params One value vector, which specifies one of two parameters to be optimized:
#' \enumerate{
#'   \item{softPrefValue is optimized, i.e. The strength of "preferring one entity over others". (The larger the value the higher the tendency towards uniform liking)}
#'   \item{non-obedience fixed at 0.1, i.e. The extent to which the instruction of the speaker is obeyed by the listener.
#'   (0 = full obedience, infinity = full instruction ignorance)}
#'}
#'
#' @param data A Matrix with data rows.
#'
#' column structure:
#'
#' [1:OC1,OC2,OC3,4:UUFeat, 5:Q1Feat,6:Q2Feat]
#'
#' [7:Q1AnswerV1,V2,V3, 10:Q2AnswerV1,V2,V3]
#'
#' \strong{1:OC1} Object 1. A value between 1 and 27.
#'
#' \strong{2:OC2} Object 2. A value between 1 and 27.
#'
#' \strong{3:OC3} Object 3. A value between 1 and 27.
#'
#' \strong{4:UUFeat} Uttered feature. A number between 1 and 3. (1: shape, 2: pattern, 3: color)
#'
#' \strong{5:Q1Feat} Questioned feature 1. A number between 1 and 3. (1: shape, 2: pattern, 3: color).
#'
#' Example: If you utter "blue" (feature: color), then you can learn something about shape and texture preferences.
#'
#' \strong{6:Q2Feat} Questioned feature 2. A number between 1 and 3. (1: shape, 2: pattern, 3: color).
#'
#' Example: If you utter "blue" (feature: color), then you can learn something about shape and texture preferences.
#'
#' \strong{7:Q1AnswerV1, V2, V3} The columns 7-9 contain the participants' slider values for the first questioned feature.
#'
#' \strong{10:Q2AnswerV1, V2, C3} The columns 10-12 contain the participants' slider values for the second questioned feature.
#'
#' @return Minimized Kullback-Leibler divergence and the optimal parameter values.
#' @details
#' This function uses \code{\link{RSAModelKLDiv3params_simpleRSA4TrialsIterative_dep}}.
#' @export
LL1_1_Iterative_dep_notObey0.1 <- function(params, data) {
  return(RSAModelKLDiv3params_simpleRSA4TrialsIterative_dep(data, abs(params[1]), 0.1))
}

# ---- optimizing 2nd parameter in iterative model: obedience ----

#' Cost function for one parameter optimization (iterative setting, dependent on trial order).
#' Optimizing non-obedience.
#' Softness is fixed at 0.
#'
#' @description
#' Simple RSA
#'
#' 1 parameter optimization; The softness parameter is optimized.
#'
#' The non-obedience parameter is fixed.
#' @param params One value vector, which specifies one of two parameters to be optimized:
#' \enumerate{
#'   \item{softPrefValue is fixed to 0, i.e. The strength of "preferring one entity over others". (The larger the value the higher the tendency towards uniform liking)}
#'   \item{non-obedience is optimized i.e. The extent to which the instruction of the speaker is obeyed by the listener.
#'   (0 = full obedience, infinity = full instruction ignorance)}
#'}
#'
#' @param data A Matrix with data rows.
#'
#' column structure:
#'
#' [1:OC1,OC2,OC3,4:UUFeat, 5:Q1Feat,6:Q2Feat]
#'
#' [7:Q1AnswerV1,V2,V3, 10:Q2AnswerV1,V2,V3]
#'
#' \strong{1:OC1} Object 1. A value between 1 and 27.
#'
#' \strong{2:OC2} Object 2. A value between 1 and 27.
#'
#' \strong{3:OC3} Object 3. A value between 1 and 27.
#'
#' \strong{4:UUFeat} Uttered feature. A number between 1 and 3. (1: shape, 2: pattern, 3: color)
#'
#' \strong{5:Q1Feat} Questioned feature 1. A number between 1 and 3. (1: shape, 2: pattern, 3: color).
#'
#' Example: If you utter "blue" (feature: color), then you can learn something about shape and texture preferences.
#'
#' \strong{6:Q2Feat} Questioned feature 2. A number between 1 and 3. (1: shape, 2: pattern, 3: color).
#'
#' Example: If you utter "blue" (feature: color), then you can learn something about shape and texture preferences.
#'
#' \strong{7:Q1AnswerV1, V2, V3} The columns 7-9 contain the participants' slider values for the first questioned feature.
#'
#' \strong{10:Q2AnswerV1, V2, C3} The columns 10-12 contain the participants' slider values for the second questioned feature.
#'
#' @return Minimized Kullback-Leibler divergence and the optimal parameter values.
#' @details
#' This function uses \code{\link{RSAModelKLDiv3params_simpleRSA4TrialsIterative_dep}}.
#' @export
LL1_2_Iterative_dep_pref0 <- function(params, data) {
  return(RSAModelKLDiv3params_simpleRSA4TrialsIterative_dep(data, 0, abs(params[1])))
}

# ---- two parameter optimization ----
# optimizing 1st and 3rd parameters in iterative model: softness + prior rate
# obedience fixed at 0

#' Cost function for one parameter optimization (iterative setting, dependent on trial order).
#' Optimizing softness and non-obedience.
#'
#' @description
#' Simple RSA
#'
#' 2 parameter optimization; The softness and non-obedience parameters are optimized.
#'
#' The non-obedience parameter is fixed.
#' @param params One value vector, which specifies one of two parameters to be optimized:
#' \enumerate{
#'   \item{softPrefValue is optimized, i.e. The strength of "preferring one entity over others". (The larger the value the higher the tendency towards uniform liking)}
#'   \item{non-obedience is optimized, i.e. The extent to which the instruction of the speaker is obeyed by the listener.
#'   (0 = full obedience, infinity = full instruction ignorance)}
#'}
#'
#' @param data A Matrix with data rows.
#'
#' column structure:
#'
#' [1:OC1,OC2,OC3,4:UUFeat, 5:Q1Feat,6:Q2Feat]
#'
#' [7:Q1AnswerV1,V2,V3, 10:Q2AnswerV1,V2,V3]
#'
#' \strong{1:OC1} Object 1. A value between 1 and 27.
#'
#' \strong{2:OC2} Object 2. A value between 1 and 27.
#'
#' \strong{3:OC3} Object 3. A value between 1 and 27.
#'
#' \strong{4:UUFeat} Uttered feature. A number between 1 and 3. (1: shape, 2: pattern, 3: color)
#'
#' \strong{5:Q1Feat} Questioned feature 1. A number between 1 and 3. (1: shape, 2: pattern, 3: color).
#'
#' Example: If you utter "blue" (feature: color), then you can learn something about shape and texture preferences.
#'
#' \strong{6:Q2Feat} Questioned feature 2. A number between 1 and 3. (1: shape, 2: pattern, 3: color).
#'
#' Example: If you utter "blue" (feature: color), then you can learn something about shape and texture preferences.
#'
#' \strong{7:Q1AnswerV1, V2, V3} The columns 7-9 contain the participants' slider values for the first questioned feature.
#'
#' \strong{10:Q2AnswerV1, V2, C3} The columns 10-12 contain the participants' slider values for the second questioned feature.
#'
#' @return Minimized Kullback-Leibler divergence and the optimal parameter values.
#' @details
#' This function uses \code{\link{RSAModelKLDiv3params_simpleRSA4TrialsIterative_dep}}.
#' @export
LL2_12_Iterative_dep <- function(params, data) {
  return(RSAModelKLDiv3params_simpleRSA4TrialsIterative_dep(data, abs(params[1]), abs(params[2])))
}
