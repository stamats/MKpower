\name{power.hsu.t.test}
\alias{power.hsu.t.test}
\alias{power.hsu.t.test2}
\alias{ssize.hsu.t.test}
\title{Power Calculations for Two-sample Hsu t Test}
\usage{
power.hsu.t.test(n = NULL, delta = NULL, sd1 = 1, sd2 = 1, sig.level = 0.05,
                   power = NULL, alternative = c("two.sided", "one.sided"),
                   strict = FALSE, tol = .Machine$double.eps^0.25)
power.hsu.t.test2(n = NULL, delta = NULL, sd1 = 1, sd2 = 1, sig.level = 0.05,
                   power = NULL, alternative = c("two.sided", "one.sided"),
                   strict = FALSE, tol = .Machine$double.eps^0.25, 
                    abstol = .Machine$double.eps^0.5)
ssize.hsu.t.test(delta = 1, sd1 = 1, sd2 = 1, sig.level = 0.05, power = 0.8,
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
  Compute the power of the two-sample Hsu t test, or determine parameters
  to obtain a target power; see Section 7.4.4 in Hedderich and Sachs (2016).
}

\details{
  In case of \code{power.hsu.t.test} and \code{power.hsu.t.test2}, 
  exactly one of the parameters \code{n}, \code{delta}, \code{power},
  \code{sd1}, \code{sd2} and \code{sig.level} must be passed as \code{NULL},
  and that parameter is determined from the others. Notice that the last three
  have non-NULL defaults, so NULL must be explicitly passed if you want to
  compute them.

  In the case of the \code{power.hsu.t.test} function, both groups are assumed 
  to be of equal size (i.e. balanced design), whereas the \code{power.hsu.t.test2} 
  function allows for unequal sample sizes per group (i.e. unbalanced design). 
  For unequal variances, the balanced design does not give the minimum total 
  sample size (sum of the group sample sizes).
  The Function \code{ssize.hsu.t.test} returns the optimal group sample sizes 
  by minimizing the total sample size (sum of the group sample sizes) that achieves
  the given power. The solution of function \code{power.hsu.t.test2} might not 
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
  J. Hedderich, L. Sachs. \emph{Angewandte Statistik: Methodensammlung mit R}.
  Springer 2016.
}

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}

\note{
  The function and its documentation was adapted from \code{power.t.test}
  implemented by Peter Dalgaard and based on previous work by Claus Ekstroem.

  \code{uniroot} is used to solve the power equation for unknowns, so
  you may see errors from it, notably about inability to bracket the
  root when invalid arguments are given.
}
\seealso{\code{\link{power.welch.t.test}}, \code{\link{power.t.test}},
\code{\link{t.test}}, \code{\link{uniroot}}}

\examples{
 ## more conservative than classical or Welch t-test
 power.t.test(n = 20, delta = 1)
 power.welch.t.test(n = 20, delta = 1)
 power.hsu.t.test(n = 20, delta = 1)
 power.hsu.t.test(power = .90, delta = 1)
 power.hsu.t.test(power = .90, delta = 1, alternative = "one.sided")

 ## sd1 = 0.5, sd2 = 1
 power.welch.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9)
 power.hsu.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9)
 power.welch.t.test2(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9) # 22+87 = 109
 power.hsu.t.test2(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9) # 44+59 = 103
 power.hsu.t.test2(n = c(44, 59), delta = 0.5, sd1 = 0.5, sd2 = 1)
 ssize.welch.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9) # 35+62 = 97
 ssize.hsu.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9) # 39+62 = 101
 power.hsu.t.test2(n = c(39, 62), delta = 0.5, sd1 = 0.5, sd2 = 1)
 power.welch.t.test2(n = c(39, 62), delta = 0.5, sd1 = 0.5, sd2 = 1)
}
\keyword{htest}


