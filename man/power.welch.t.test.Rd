\name{power.welch.t.test}
\alias{power.welch.t.test}
\alias{power.welch.t.test2}
\alias{ssize.welch.t.test}
\title{Power and Sample Size Calculations for Welch Two-sample t Test}
\usage{
power.welch.t.test(n = NULL, delta = NULL, sd1 = 1, sd2 = 1, sig.level = 0.05,
                   power = NULL, alternative = c("two.sided", "one.sided"),
                   strict = FALSE, tol = .Machine$double.eps^0.25)
power.welch.t.test2(n = NULL, delta = NULL, sd1 = 1, sd2 = 1, sig.level = 0.05,
                    power = NULL, alternative = c("two.sided", "one.sided"),
                    strict = FALSE, tol = .Machine$double.eps^0.25, 
                    abstol = .Machine$double.eps^0.5)
ssize.welch.t.test(delta = 1, sd1 = 1, sd2 = 1, sig.level = 0.05, power = 0.8,
                   alternative = c("two.sided", "one.sided"), strict = FALSE)                    
}
\arguments{
 \item{n}{number of observations (per group)}
 \item{delta}{(expected) true difference in means}
 \item{sd1}{(expected) standard deviation of group 1}
 \item{sd2}{(expected) standard deviation of group 2}
 \item{sig.level}{significance level (Type I error probability)}
 \item{power}{power of test (1 minus Type II error probability)}
 \item{alternative}{one- or two-sided test.  Can be abbreviated.}
 \item{strict}{use strict interpretation in two-sided case}
 \item{tol}{numerical tolerance used in root finding, the default
    providing (at least) four significant digits.}
 \item{abstol}{controls the absolute convergence tolerance of the 
    \code{"Nelder-Mead"} method. Only useful for non-negative functions, 
    as a tolerance for reaching zero; see \code{\link{optim}}.}
}
\description{
  Compute the power of the two-sample Welch t test, or determine parameters
  to obtain a target power.
}

\details{
  In case of \code{power.welch.t.test} and \code{power.welch.t.test2}, exactly 
  one of the parameters \code{n}, \code{delta}, \code{power},
  \code{sd1}, \code{sd2} and \code{sig.level} must be passed as \code{NULL},
  and that parameter is determined from the others. Notice that the last three
  have non-NULL defaults, so NULL must be explicitly passed if you want to
  compute them.
  
  In the case of the \code{power.welch.t.test} function, both groups are assumed 
  to be of equal size (i.e. balanced design), whereas the \code{power.welch.t.test2} 
  function allows for unequal sample sizes per group (i.e. unbalanced design). 
  For unequal variances, the balanced design does not give the minimum total 
  sample size (sum of the group sample sizes).
  The Function \code{ssize.welch.t.test} returns the optimal group sample sizes 
  by minimizing the total sample size (sum of the group sample sizes) that achieves
  the given power. The solution of function \code{power.welch.t.test2} might not 
  represent the global optimum in terms of total sample size.

  If \code{strict = TRUE} is used, the power will include the probability of
  rejection in the opposite direction of the true effect, in the two-sided
  case. Without this the power will be half the significance level if the
  true difference is zero.
}
\value{
  Object of class \code{"power.htest"}, a list of the arguments
  (including the computed one) augmented with \code{method} and
  \code{note} elements.
}

\references{
  S.L. Jan and G. Shieh (2011). Optimal sample sizes for Welch's test under
  various allocation and cost considerations. \emph{Behav Res Methods}, 43,
  4:1014-22.
}

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}

\note{
  The function and its documentation was adapted from \code{power.t.test}
  implemented by Peter Dalgaard and based on previous work by Claus Ekstroem.

  \code{uniroot} is used to solve the power equation for unknowns, so
  you may see errors from it, notably about inability to bracket the
  root when invalid arguments are given.
}
\seealso{\code{\link{power.t.test}}, \code{\link{t.test}}, \code{\link{uniroot}}}

\examples{
 ## identical results as power.t.test, since sd = sd1 = sd2 = 1
 power.welch.t.test(n = 20, delta = 1)
 power.welch.t.test(power = .90, delta = 1)
 power.welch.t.test(power = .90, delta = 1, alternative = "one.sided")
 
 power.welch.t.test2(n = c(20, 20), delta = 1)
 power.welch.t.test2(power = .90, delta = 1)
 power.welch.t.test2(power = .90, delta = 1, alternative = "one.sided")
 
 ssize.welch.t.test(power = 0.90, delta = 1)
 power.welch.t.test2(n = c(25, 20), delta = 1)
 power.welch.t.test2(n = c(20, 25), delta = 1)
 power.welch.t.test(n = 23, delta = 1)
 power.welch.t.test2(n = c(23, 23), delta = 1)
 
 ssize.welch.t.test(power = .90, delta = 1, alternative = "one.sided")
 power.welch.t.test2(n = c(19, 17), delta = 1, alternative = "one.sided")
 power.welch.t.test2(n = c(17, 19), delta = 1, alternative = "one.sided")
 power.welch.t.test(n = 18, delta = 1, alternative = "one.sided")
 power.welch.t.test2(n = c(18, 18), delta = 1, alternative = "one.sided")

 ###########################################################
 ## see Table 1 of Jan and Shieh (2011)
 ## delta = 1, sig.level = 0.05, power = 0.9
 ###########################################################

 #######################################
 ## sd1 = 1/3, sd2 = 1
 #######################################
 power.welch.t.test(delta = 1, sd1 = 1/3, sd2 = 1, power = 0.9)
 power.welch.t.test2(delta = 1, sd1 = 1/3, sd2 = 1, power = 0.9)
 ssize.welch.t.test(delta = 1, sd1 = 1/3, sd2 = 1, power = 0.9)
 
 ## n1 = 14, n2 = 14
 power.welch.t.test(n = 14, delta = 1, sd1 = 1/3, sd2 = 1, power = NULL)
 power.welch.t.test2(n = c(14, 14), delta = 1, sd1 = 1/3, sd2 = 1, power = NULL)
 
 ## n1 = 8, n2 = 14
 power.welch.t.test2(n = c(8, 16), delta = 1, sd1 = 1/3, sd2 = 1, power = NULL)
 
 ## n1 = 6, n2 = 18
 power.welch.t.test2(n = c(6, 18), delta = 1, sd1 = 1/3, sd2 = 1, power = NULL)

 #######################################
 ## sd1 = 1/2, sd2 = 1
 #######################################
 power.welch.t.test(delta = 1, sd1 = 1/2, sd2 = 1, power = 0.9)
 power.welch.t.test2(delta = 1, sd1 = 1/2, sd2 = 1, power = 0.9)
 ssize.welch.t.test(delta = 1, sd1 = 1/2, sd2 = 1, power = 0.9)
 
 ## n1 = 15, n2 = 15
 power.welch.t.test(n = 15, delta = 1, sd1 = 1/2, sd2 = 1, power = NULL)
 power.welch.t.test2(n = c(15, 15), delta = 1, sd1 = 1/2, sd2 = 1, power = NULL)
 
 ## n1 = 9, n2 = 18
 power.welch.t.test2(n = c(9, 18), delta = 1, sd1 = 1/2, sd2 = 1, power = NULL)
 
 ## n1 = 7, n2 = 21
 power.welch.t.test2(n = c(7, 21), delta = 1, sd1 = 1/2, sd2 = 1, power = NULL)

 #######################################
 ## sd1 = 1, sd2 = 1
 #######################################
 power.t.test(delta = 1, sd = 1, power = 0.9)
 power.welch.t.test(delta = 1, sd1 = 1, sd2 = 1, power = 0.9)
 power.welch.t.test2(delta = 1, sd1 = 1, sd2 = 1, power = 0.9)
 ssize.welch.t.test(delta = 1, sd1 = 1, sd2 = 1, power = 0.9)
 
 ## n1 = 23, n2 = 23
 power.t.test(n = 23, delta = 1, sd = 1, power = NULL)
 power.welch.t.test(n = 23, delta = 1, sd1 = 1, sd2 = 1, power = NULL)
 power.welch.t.test2(n = c(23, 23), delta = 1, sd1 = 1, sd2 = 1, power = NULL)
 
 ## n1 = 17, n2 = 34
 power.welch.t.test2(n = c(17, 34), delta = 1, sd1 = 1, sd2 = 1, power = NULL)
 
 ## n1 = 16, n2 = 48
 power.welch.t.test2(n = c(16, 48), delta = 1, sd1 = 1, sd2 = 1, power = NULL)

 #######################################
 ## sd1 = 2, sd2 = 1
 #######################################
 power.welch.t.test(delta = 1, sd1 = 2, sd2 = 1, power = 0.9)
 power.welch.t.test2(delta = 1, sd1 = 2, sd2 = 1, power = 0.9)
 ssize.welch.t.test(delta = 1, sd1 = 2, sd2 = 1, power = 0.9)
 
 ## n1 = 54, n2 = 54
 power.welch.t.test(n = 54, delta = 1, sd1 = 2, sd2 = 1, power = NULL)
 power.welch.t.test2(n = c(54, 54), delta = 1, sd1 = 2, sd2 = 1, power = NULL)
 
 ## n1 = 49, n2 = 98 
 power.welch.t.test2(n = c(49, 98), delta = 1, sd1 = 2, sd2 = 1, power = NULL)
 
 ## n1 = 48, n2 = 144
 power.welch.t.test2(n = c(48, 144), delta = 1, sd1 = 2, sd2 = 1, power = NULL)

 #######################################
 ## sd1 = 3, sd2 = 1
 #######################################
 power.welch.t.test(delta = 1, sd1 = 3, sd2 = 1, power = 0.9)
 power.welch.t.test2(delta = 1, sd1 = 3, sd2 = 1, power = 0.9)
 ssize.welch.t.test(delta = 1, sd1 = 3, sd2 = 1, power = 0.9)
 
 ## n1 = 107, n2 = 107
 power.welch.t.test(n = 107, delta = 1, sd1 = 3, sd2 = 1, power = NULL)
 power.welch.t.test2(n = c(107, 107), delta = 1, sd1 = 3, sd2 = 1, power = NULL)
 
 ## n1 = 102, n2 = 204
 power.welch.t.test2(n = c(102, 204), delta = 1, sd1 = 3, sd2 = 1, power = NULL)
 
 ## n1 = 100, n2 = 300
 power.welch.t.test2(n = c(100, 300), delta = 1, sd1 = 3, sd2 = 1, power = NULL)
}
\keyword{htest}


