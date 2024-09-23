\name{ssize.auc.ci}
\alias{ssize.auc.ci}
\title{Sample Size Calculations for AUC}
\description{
 Compute sample size, power, delta, or significance level of a diagnostic test
 for an expected AUC.
}
\usage{
ssize.auc.ci(AUC = NULL, delta = NULL, n = NULL, sig.level = 0.05,
             power = NULL, prev = NULL, NMAX = 1e4)
}
\arguments{
 \item{AUC}{Expected AUC.}
 \item{n}{Total sample size (number of cases + number of controls).}
 \item{delta}{\code{AUC}-\code{delta} is used as lower confidence limit}
 \item{sig.level}{Significance level (Type I error probability)}
 \item{power}{Assurance probability of confidence interval (1 minus Type II error probability)}
 \item{prev}{Expected prevalence, if \code{NULL} prevalence is ignored which means \code{prev = 0.5}
 is assumed.}
 \item{NMAX}{Maximum sample size considered.}
}

\details{
  Exactly one of the parameters \code{n}, \code{delta}, \code{sig.level}, 
  and \code{power} must be passed as \code{NULL}, and that parameter is determined 
  from the others. Notice that \code{sig.level} has a non-\code{NULL} default 
  so \code{NULL} must be explicitly passed if you want to compute it.

  The computations use the variance of the AUC derived by Hanley and McNeil (1982)
  and incorporate an additional assurance probability (power) as in the approach
  of Flahault et al. (2005).
  
  As noted in Chu and Cole (2007) power is not a monotonically increasing
  function in n but rather saw toothed (see also Chernick and Liu (2002)).
  Hence, in our calculations we use the more conservative approach II); 
  i.e., the minimum sample size \code{n} such that the actual power is 
  larger or equal \code{power} andsuch that for any sample size larger 
  than \code{n} it also holds that the actual power is larger or equal 
  \code{power}.
}
\value{
  Object of class \code{"power.htest"}, a list of the arguments
  (including the computed one) augmented with \code{method} and
  \code{note} elements.
}

\references{ 
  A. Flahault, M. Cadilhac, and G. Thomas (2005). Sample size calculation 
  should be performed for design accuracy in diagnostic test studies. 
  \emph{Journal of Clinical Epidemiology}, \bold{58}(8):859-862.
  
  J.A. Hanley, B.J. McNeil (1982). The meaning and use of the area under a 
  receiver operating characteristic (ROC) curve. \emph{Radiology}, \bold{143}(1):29-36.
}

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}

\note{
  \code{uniroot} is used to solve the equations for unknowns, so
  you may see errors from it, notably about inability to bracket the
  root when invalid arguments are given.
}
\seealso{\code{\link{uniroot}}}

\examples{
## compute n
ssize.auc.ci(AUC = 0.9, delta = 0.05, power = 0.8)
## compute delta
ssize.auc.ci(AUC = 0.9, n = 254, power = 0.8)
## compute power
ssize.auc.ci(AUC = 0.9, n = 254, delta = 0.05)
## compute sig.level
ssize.auc.ci(AUC = 0.9, n = 254, delta = 0.05, power = 0.8, sig.level = NULL)
}
\keyword{htest}


