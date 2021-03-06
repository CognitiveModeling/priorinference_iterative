% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sRSA_iterative_indep_pr.R
\name{LL1_1_Iterative_pr_notObey0_pr0.5}
\alias{LL1_1_Iterative_pr_notObey0_pr0.5}
\title{Cost function for one parameter optimization (iterative setting, independent of trial order).
Optimizing softness.
non-obedience is fixed at 0.
prior rate is fixed at 0.5.}
\usage{
LL1_1_Iterative_pr_notObey0_pr0.5(params, data)
}
\arguments{
\item{params}{One value vector, which specifies one of three parameters to be optimized:
\enumerate{
  \item{softPrefValue is optimized, i.e. The strength of "preferring one entity over others". (The larger the value the higher the tendency towards uniform liking)}
  \item{non-obedience is fixed at 0, i.e. The extent to which the instruction of the speaker is obeyed by the listener.
  (0 = full obedience, infinity = full instruction ignorance)}
  \item{priorRate is fixed to 0.5. This parameter specifies how much the prior information is weighed into the decision .}
}}

\item{data}{A Matrix with data rows.

column structure:

[1:OC1,OC2,OC3,4:UUFeat, 5:Q1Feat,6:Q2Feat]

[7:Q1AnswerV1,V2,V3, 10:Q2AnswerV1,V2,V3]

\strong{1:OC1} Object 1. A value between 1 and 27.

\strong{2:OC2} Object 2. A value between 1 and 27.

\strong{3:OC3} Object 3. A value between 1 and 27.

\strong{4:UUFeat} Uttered feature. A number between 1 and 3. (1: shape, 2: pattern, 3: color)

\strong{5:Q1Feat} Questioned feature 1. A number between 1 and 3. (1: shape, 2: pattern, 3: color).

Example: If you utter "blue" (feature: color), then you can learn something about shape and texture preferences.

\strong{6:Q2Feat} Questioned feature 2. A number between 1 and 3. (1: shape, 2: pattern, 3: color).

Example: If you utter "blue" (feature: color), then you can learn something about shape and texture preferences.

\strong{7:Q1AnswerV1, V2, V3} The columns 7-9 contain the participants' slider values for the first questioned feature.

\strong{10:Q2AnswerV1, V2, C3} The columns 10-12 contain the participants' slider values for the second questioned feature.}
}
\value{
Minimized Kullback-Leibler divergence and the optimal parameter values.
}
\description{
Simple RSA

1 parameter optimization; The softness parameter is optimized. (1st)

The non-obedience and prior rate parameter are fixed.
}
\details{
This function uses \code{\link{RSAModelKLDiv3params_simpleRSA4TrialsIterative_pr}}.
}
