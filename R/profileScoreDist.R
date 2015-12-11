#' @useDynLib profileScoreDist
#' @importFrom Rcpp sourceCpp
#'
#'
#'

######################################################################

pwmFormatCheck <- function (inputmat) {
  ## checks user input matrix (on the form in MotifDB) and converts to internal
  ## form (Lx4, as in the Rahmann et al. paper)
  if (nrow(inputmat) != 4L)
    stop("First argument must be a matrix with 4 rows, corresponding to
         A, C, G, and T")
  if (ncol(inputmat) < 2)
    stop("First argument must be a matrix with at least 2 columns,
         corresponding to a PWM or PCM of at least length 2")
  if (is.null(rownames(inputmat))) {
    toreturn <- t(inputmat)
    colnames(toreturn) <- c("A", "C", "G", "T")
    toreturn
  } else {
    if (!identical(sort(rownames(inputmat)), c("A", "C", "G", "T")))
      stop("Matrix rownames must be one of each of A, C, G, and T")
    t(inputmat[c("A", "C", "G", "T"), ])
  }
}

computeDelta <- function (w, tau, rho, N) {
  delta <- (1 - w)*tau + w*rho
  relent <- delta*log(delta/rho)
  ## calculate Delta(w): discard NaNs, delta approaches zero faster
  ## than ln(delta) approaches -Inf:
  ## lim_{x->0} x*log(x/c) = 0.
  2*N*sum(relent[!is.nan(relent)])
}

computeDeltadiff <- function (w, tau, rho, N, E, Dzero) {
  computeDelta(w, tau, rho, N) + E - Dzero
}

computew <- function (tau, rho, N, E) {
  ## calculate Delta(0)
  Dzero <- computeDelta(0, tau, rho, N)
  if (Dzero >= E) {
    ## find the w in (0,1) for which Delta(w) = Delta(0) - E
    stats::uniroot(computeDeltadiff, c(0, 1),
                   tau, rho, N, E, Dzero, tol=1e-9)$root
  } else
    ## Delta(0)<E => take w=1
    as.integer(1)
}

regularizeRow <- function (Ci, rho, E) {
  Ni <- sum(Ci)
  taui <- Ci/Ni
  w <- computew(taui, rho, Ni, E)
  if (w == as.integer(1))
    ## in this case, the regularized row is drawn from rho only.
    ## the number of counts can be chosen freely, we choose N.
    rho*Ni
  else {
    ## w = Wi/(Ni + Wi), which leads to:
    Wi <- w*Ni/(1 - w)
    Ni*taui + Wi*rho
  }
}

#' Careful regularization (pseudocount addition) to a position count matrix.
#'
#' Carries out the regularization suggested by Rahmann et al.  This lets each
#' column in the regularized matrix be a linear combination of the column in the
#' non-regularized matrix and rho, the overall base distribution of all
#' positions.  The weighting of the linear combination is determined by the
#' parameter E in a non-trivial way, see Rahmann et al. for more information.  A
#' default value E=1.5 usually works well.
#'
#' @param motif A position count matrix; each column a position and each row a
#'   base corresponding to A, C, G, T.  This order is assumed, unless the rows
#'   are correspondingly named in a different order.
#' @param E Weighting parameter between 0 and 3 for the regularization.
#' @return The regularized matrix
#' @examples
#' data(INR)
#' regularizeMatrix(INR)
#' @references Rahmann, S., Mueller, T., and Vingron, M. (2003). On the power of
#'   profiles for transcription factor binding site detection. Stat Appl Genet
#'   Mol Biol 2, Article7.
#' @export
#'
regularizeMatrix <- function (motif, E=1.5) {
  internalmotif <- pwmFormatCheck(motif)
  if (as.integer(round(sum(internalmotif))) == nrow(internalmotif))
    stop("You have likely supplied a position frequency matrix, not a position
         count matrix.  Consider multiplying your frequency matrix with the
         number of counts for each position.")
  if (as.integer(round(sum(internalmotif)))/nrow(internalmotif) < 10)
    warning("Your position count matrix consists of very few counts.  Results
            may be unreliable.")
  if (E > 3 || E < 0)
    stop("The E weighting parameter has to lie between 0 and 3.")
  ## the regularizing distribution: columnwise division
  rho <- apply(internalmotif, 2, sum)/sum(internalmotif)
  ## for the rare case that a rho is zero, add to numerator and denominator
  if (any(rho <= 1e-15))
    rho <- (apply(internalmotif, 2, sum) + 0.25) / (sum(internalmotif) + 1)
  apply(internalmotif, 1, regularizeRow, rho, E)
}

#' Compute exact position weight/count matrix score distribution.
#'
#' Computes the discretisized score distribution of a position count matrix
#' (PCM) or a position weight matrix (PWM), using the method described by
#' Rahmann et al.
#'
#' @param motif A matrix representing a PCM or PWM; each column a position and
#'   each row a base corresponding to A, C, G, T.  This order is assumed, unless
#'   the rows are correspondingly named in a different order.
#' @param gc A scalar giving the GC fraction to assume.
#' @param granularity The granularity of the discretization, defaults to 0.01.
#' @param unit The logarithm unit of the score computed from the PCM or PWM, can
#'   be "nat" (default, natural logarithm), "bit" (base 2), or "dit" (base 10).
#' @return a ProfileDist object
#' @examples
#' data(INR)
#' thedist <- computeScoreDist(regularizeMatrix(INR), 0.5)
#' plotDist(thedist)
#' @references Rahmann, S., Mueller, T., and Vingron, M. (2003). On the power of
#'   profiles for transcription factor binding site detection. Stat Appl Genet
#'   Mol Biol 2, Article7.
#' @export
#'
computeScoreDist <- function (motif, gc, granularity=0.01, unit="nat") {
  internalmotif <- pwmFormatCheck(motif)
  ## fail if not regularized
  if(any(internalmotif==0L))
    stop(paste("The profile matrix contains zeroes.",
               "Use e.g. regularizeMatrix() to add pseudocounts."))
  ## motif length
  L = nrow(internalmotif)
  ## convert counts to probabilities (check ?/ for explanation)
  pmotif <- internalmotif/rowSums(internalmotif)
  ## vector of background probabilities for ACGT
  bg <- c((1 - gc)/2, gc/2, gc/2, (1 - gc)/2)
  ## replicate: same background for all positions
  pbg <- matrix(rep(bg, L), nrow=L, byrow=TRUE)
  ## log-score matrix
  logpssm <- switch(unit,
                    nat = log(pmotif/pbg),
                    bit = log2(pmotif/pbg),
                    dit = log10(pmotif/pbg),
                    stop("Invalid information unit"))
  ## round S scaled by granularity
  S <- round(logpssm/granularity)
  ## minimal and maximal possible total scores for all words
  Smin <- as.integer(sum(apply(S, 1, min)))
  Smax <- as.integer(sum(apply(S, 1, max)))
  ## initialize vector of probabilities f and g: they correspond to
  ## the probability masses partial sums of scores
  ## (Smin, Smin + granularity, ..., Smax) generated under the matrix
  ## probabilities and background probabilities, respectively
  Scores <- Smin:Smax
  ## first step (calculate f_1):
  f <- profilePosscoreprob(Scores, S[1, ], pmotif[1, ])
  g <- profilePosscoreprob(Scores, S[1, ], bg)
  ## iteratively calculate the partial probability mass sums, using
  ## the convolution method
  for (psum in 2:L) {
    fi <- profilePosscoreprob(Scores, S[psum, ], pmotif[psum, ])
    f  <- profileConvolution(f, fi, abs(Smin), Smax)
    gi <- profilePosscoreprob(Scores, S[psum, ], bg)
    g  <- profileConvolution(g, gi, abs(Smin), Smax)
  }
  ## return results in the list dists
  ProfileDist(f=f, g=g, Scores=Scores*granularity)
}


#' False discovery rate and power for PWM Score distributions.
#'
#' Computes score cutoffs for a PWM or a PCM, given distributions as calculated
#' with \code{computeScoreDist()}.  Cutoffs can be computed for a given false
#' discovery rate (FDR), for a given false negative rate (FNR), and the optimal
#' tradeoff between the two, in the sense that \eqn{c \times FDR = FNR} for some
#' \eqn{c} that the user may choose.
#' @param scoreDist A ProfileDist object, as computed by
#'   \code{computeScoreDist()}
#' @param n The number of scores considered for the given PWM.  If one sequence
#'   is considered and a score is computed for all overlapping windows of the
#'   same length as the PWM, this will be the length of the sequence, minus the
#'   PWM length plus 1.  If scanning a sequence and its reverse complement too,
#'   this number must be further multiplied by two.  The number forms the basis
#'   for the FDR, since this is a multiple testing problem.
#' @param m The number of true positives assumed for computing the FNR.
#' @param c A factor expressing how much more important the FDR is compared to
#'   the FNR, when computing the tradeoff cutoff that considers both FDR and
#'   FNR.  See Rahmann et al. for details.
#' @param cutoff The FDR and FNR considered, typically 0.01 or 0.05.
#' @return a list with elements: \describe{ \item{cutoffa}{Score cutoff for
#'   FDR=\code{cutoff}} \item{cutoffb}{Score cutoff for FNR=\code{cutoff}}
#'   \item{cutoffopt}{Score cutoff for \code{c}*FDR = FNR} }
#' @examples
#' data(INR)
#' thedist <- computeScoreDist(regularizeMatrix(INR), 0.5)
#' scoreDistCutoffs(thedist, n=2000, cutoff=0.05)
#' @references Rahmann, S., Mueller, T., and Vingron, M. (2003). On the power of
#'   profiles for transcription factor binding site detection. Stat Appl Genet
#'   Mol Biol 2, Article7.
#' @export
#'
scoreDistCutoffs <- function (scoreDist, n, m=1, c=1, cutoff=0.01) {
  ## alpha is the right cumulative distribution function of g,
  ## i.e. P(X >= t): the probability that the background nucleotide
  ## distribution gives rise to a score at or above the threshold t
  alphalesseq <- cumsum(backgroundDist(scoreDist))
  ## these are the probabilities of less than
  alphaless <- c(0, alphalesseq[1:(length(alphalesseq) - 1)])
  ## we need the probability of greater than or equal
  alpha  <- 1 - alphaless
  alphan <- 1 - exp(-n*alpha)
  alphaind <- which.min(abs(alphan - cutoff))
  ## beta is the left cumulative distribution function of f,
  ## i.e. P(X < t): the probability that the motif nucleotide
  ## distribution gives rise to a score below the threshold t
  betalesseq <- cumsum(signalDist(scoreDist))
  ## these are the probabilities of less than
  beta <- c(0, betalesseq[1:(length(betalesseq) - 1)])
  betam <- 1 - (1 - beta)^m
  betaind <- which.min(abs(betam - cutoff))
  ## find out the index where c*alphan is closest to betam.
  aboptind <- which.min(abs(c*alphan - betam))
  list(cutoffa=score(scoreDist)[alphaind],
       cutoffb=score(scoreDist)[betaind],
       cutoffopt=score(scoreDist)[aboptind])
}

