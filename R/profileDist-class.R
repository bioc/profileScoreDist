#' ProfileDist
#'
#' This class represents signal and background score distributions for a
#' profile.
#'
#' @section Constructor:
#' \code{ProfileDist(f=numeric, g=numeric, Scores=numeric)}
#'
#' @slot f Signal distribution
#' @slot g Background distribution
#' @slot Scores Scores for the distributions
#'
#' @param x A ProfileDist object.
#' @param object A ProfileDist object for the \code{show} method.
#'
#' @return A ProfileDist object.
#'
#' @export backgroundDist signalDist plotDist ProfileDist
#' @import BiocGenerics methods
#' @importFrom graphics legend lines plot title
#'
# Class definition --------------------------------------------------------

ProfileDist <- setClass("ProfileDist",
                        slots=c(f="numeric", g="numeric", Scores="numeric"),
                        prototype=list(f=1.0, g=1.0, Scores=0),
                        validity=function (object) {
                          lengthtest <- length(object@f) == length(object@g) &&
                            length(object@g) == length(object@Scores)
                          vectest <- all(is.vector(object@f),
                                         is.vector(object@g),
                                         is.vector(object@Scores))
                          sumtest <- all(zapsmall(sum(object@f)) == 1L,
                                         zapsmall(sum(object@g)) == 1L)
                          if (lengthtest && vectest && sumtest)
                            return(TRUE)
                          else if (!lengthtest || !vectest)
                            return(paste("Scores and probability distributions",
                                         "f and g need to be vectors of",
                                         "equal lengths."))
                          else
                            return(paste("Probability distributions f and g",
                                         "should sum to 1, within numerical",
                                         "error as defined by zapsmall()."))
                        })


# Generics ----------------------------------------------------------------

#' Background distribution.
#'
#' \code{backgroundDist} returns the background distribution of a profile
#' object.
#'
#' This is a generic function.
#' @param x A ProfileDist object.
#' @return The background distribution vector.
#' @examples
#' anObject <- ProfileDist()
#' backgroundDist(anObject)
setGeneric("backgroundDist", function (x) {
  standardGeneric("backgroundDist")
})

#' Signal distribution.
#'
#' \code{signalDist} returns the signal distribution of a profile
#' object.
#'
#' This is a generic function.
#' @param x A ProfileDist object.
#' @return The signal distribution vector.
#' @examples
#' anObject <- ProfileDist()
#' backgroundDist(anObject)
setGeneric("signalDist", function (x) {
  standardGeneric("signalDist")
})

#' Plot background and signal distributions.
#'
#' \code{plotDist} creates a rudimentary plot of signals and backgrounds.
#'
#' This is a generic function.
#' @param x A ProfileDist object.
#' @return The scores vector.
#' @examples
#' data(INR)
#' thedist <- computeScoreDist(regularizeMatrix(INR), 0.5)
#' plotDist(thedist)
setGeneric("plotDist", function (x) {
  standardGeneric("plotDist")
})



# Methods -----------------------------------------------------------------

#' @describeIn ProfileDist Shows useful information
setMethod("show", "ProfileDist",
          function (object) {
            cat("\nProfile score distribution\n\n")
            cat(paste("Maximum score:", max(object@Scores), "\n"))
            cat(paste("Minimum score:", min(object@Scores), "\n"))
            cat(paste("Number of discrete probabilities:", length(object@f)))
          })

## Simple accessors:
#' @describeIn ProfileDist Accessor for the scores
setMethod("score", "ProfileDist",
          function (x) {
            x@Scores
          })

#' @describeIn ProfileDist Accessor for the signal distribution
setMethod("signalDist", "ProfileDist",
          function (x) {
            x@f
          })

#' @describeIn ProfileDist Accessor for the background distribution
setMethod("backgroundDist", "ProfileDist",
          function (x) {
            x@g
          })

#' @describeIn ProfileDist Simple plot method for signal and background
#'   distributions
setMethod("plotDist", "ProfileDist",
          function (x) {
            plot(x@Scores, x@f, ylim=c(0, 1.05*max(c(x@f, x@g))),
                 type="l", col="olivedrab", ann=FALSE)
            lines(x@Scores, x@g,
                  type="l", col="red")
            title(xlab="Score", ylab="Probability")
            legend("topleft", legend=c("Background", "Signal"),
                   col=c("red", "olivedrab"), lty=1)
          })

