% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sRSA_iterative_dep.R
\name{simplePragmaticSpeakerWithPrefPriorAll_depOnOrder}
\alias{simplePragmaticSpeakerWithPrefPriorAll_depOnOrder}
\title{Simple pragmatic speaker with all prior preferences iterative function
Iterative function dependent on trial order (prior rate)}
\usage{
simplePragmaticSpeakerWithPrefPriorAll_depOnOrder(
  utterance,
  obj,
  preferencesPriorAll,
  validUtterances,
  currentObjects,
  uttToObjProbs,
  objectPreferenceSoftPriors
)
}
\arguments{
\item{utterance}{The uttered word by the speaker that the listener hears.

An index referring to one of the values in the vector validUtterances.}

\item{obj}{The object chosen by the listener. A value referring to the index
1,2 or 3.}

\item{preferencesPriorAll}{A vector of length 9.

Probability mass over all feature values.

Gives a prior preferences distribution over all (nine) feature values.

\code{preferencesPriorAll <- rep(1/9, 9)}}

\item{validUtterances}{A vector of utterances that correspond to all feature
values present in the current objects in the scene.

For example, it only makes sense to utter \emph{"red"} in a scene if there
are \emph{red} objects present.}

\item{currentObjects}{Vector of three values in \code{{1,...,27}} specifying
the target and the other two objects.

The target is the first object in the vector \code{(index = 1)}.}

\item{uttToObjProbs}{A matrix. The rows map each possible utterance that
corresponds to each present feature value of the current objects. The
columns represent the three objects in the scene.

This reflects the obedience-parameter and which objects match the
respective utterance. The matrix shows the probability that a certain
object is chosen following a certain utterance, that is valid in the scene.
The number of rows of the matrix match the length of the validUtterances
vector.}

\item{objectPreferenceSoftPriors}{A list of preference priors for all valid
utterances based on the object in the scene.

The list has as many rows as the length of the validUtterances vector + 1.

Each row in the list contains a vector of length 3, as there are three
objects in the scene.

The extra row is for the case of no feature preferences whatsoever, i.e.
uniform prior over all three objects in the scene.}

\item{weights}{A vector of length 4 including the weight by which the prior is weighed.

weights <- c(0.3, 0.4, 0.5, 0.6)}

\item{trial}{A vector of length 4 including the number of the current trial.

trial <- c(1,2,3,4)}
}
\value{
A vector of length 9. It contains the normalized probability over
  preferences (priors).
}
\description{
Simple-RSA (iterative, dependent on trial order)

The simple pragmatic speaker considers all "imaginable" (i.e. implemented)
preference distributions over objects of the listener.

Starting with a prior assumption over the possible listener's preferences. It
then infers the posterior over these preferences given the listener makes a
particular object choice. P(listener's feature value preferences | utterance,
object choice by the listener, prior over listener's feature value
preferences).

This function takes the evidence from the current trial in consideration and also
the prior from the trials before and computes the posterior over the feature preferences of the listener.
\code{posterior = (1 - prior rate) x evidence + (prior rate) x prior}.
}
\details{
This is function is the second of two functions that are used in the iterative setting using the prior rate parameter.
The first one is: \code{\link{simplePragmaticSpeakerWithPrefPriorAll_indepOfOrder_pr}}.
}
\examples{
\donttest{simplePragmaticSpeakerWithPrefPriorAll_depOnOrder(utterance, obj, preferencesPriorAll,
validUtterances, currentObjects,
uttToObjProbs, objectPreferenceSoftPriors) }

output:
[1] 0.5333333 0.1333333 0.3333333 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
[9] 0.0000000
}
