\name{sim.ssize.perm.t.test}
\alias{sim.ssize.perm.t.test}
\title{Sample Size for Permutation t Tests}
\usage{
sim.ssize.perm.t.test(rx, ry = NULL, mu = 0, sig.level = 0.05, power = 0.8, 
                      type = c("two.sample", "one.sample", "paired"), 
                      alternative = c("two.sided", "less", "greater"),
                      var.equal = FALSE, R = 9999, symmetric = TRUE,
                      useCombn = FALSE, n.min = 10, n.max = 200, 
                      step.size = 10, iter = 10000, BREAK = TRUE,
                      parallel = FALSE, cl = NULL)
}
\arguments{
 \item{rx}{ function to simulate the values of x, respectively x-y in the paired case.}
 \item{ry}{ function to simulate the values of y in the two-sample case}
 \item{mu}{ true values of the location shift for the null hypothesis.}
 \item{sig.level}{ significance level (Type I error probability)}
 \item{power}{ two-sample, one-sample or paired test }
 \item{type}{ one- or two-sided test.  Can be abbreviated.}
 \item{alternative}{one- or two-sided test.  Can be abbreviated.}
 \item{var.equal}{a logical variable indicating whether to treat the
  two variances as being equal. If \code{TRUE} then the pooled
  variance is used to estimate the variance otherwise the Welch
  (or Satterthwaite) approximation to the degrees of freedom is used.}
 \item{R}{ number of (Monte-Carlo) permutations. }
 \item{symmetric}{a logical variable indicating whether to assume symmetry
  in the two-sided test. If \code{TRUE} then the symmetric permutation p value 
  otherwise the equal-tail permutation p value is computed.}
 \item{useCombn}{ logical; use all combinations instead of permutations in the
  case of the two-sample t-test; see \code{\link[MKinfer]{perm.t.test}}. }
 \item{n.min}{ integer, start value of grid search.}
 \item{n.max}{ integer, stop value of grid search.}
 \item{step.size}{ integer, step size used in the grid search.}
 \item{iter}{ integer, number of interations of the simulations.}
 \item{BREAK}{ logical, grid search stops when the emperical power is larger
  than the requested power.}
 \item{parallel}{logical; use parellel computing via function \code{\link{parRapply}}.}
 \item{cl}{a cluster object, created by package \code{\link[parallel:makeCluster]{parallel}} 
  or by package \code{\link[snow:makeCluster]{snow}}. If \code{NULL}, a new cluster 
  will be generated with function \code{\link[parallel:makeCluster]{makeCluster}} 
  using the maximum number of available cores - 1.}
}
\description{
  Simulate the empirical power of permutation t tests for computing the required 
  sample size.
}

\details{
  Functions \code{rx} and \code{ry} are used to simulate the data.
  
  We recommend a two steps procedure: In the first step, start with a wide grid 
  and find out in which range of sample size values the intended power will
  be achieved. In the second step, the interval identified in the first step 
  is used to find the sample size that leads to the required power setting
  \code{step.size = 1} and \code{BREAK = FALSE}. This approach is applied 
  in the examples below.
}
\value{
  Object of class \code{"power.htest"}, a list of the arguments
  (including the computed one) augmented with \code{method} and
  \code{note} elements.
}

\references{
  B. Efron, R.J. Tibshirani. \emph{An Introduction to the Bootstrap}.
  Chapman and Hall/CRC 1993.
  
  A. Janssen (1997). Studentized permutation tests for non-i.i.d, hypotheses 
  and the generalized Behrens-Fisher problem. 
  \emph{Statistics and Probability Letters}, \bold{36}, 9-21.
  
  A. Janssen, T. Pauls. How do bootstrap and permutation tests work?. 
  \emph{Ann. Statist.}, \bold{31}(3), 768-806.
  
  E. Chung, J.P. Romano (2013). Exact and asymptotically robust permutation
  tests. \emph{The Annals of Statistics}, \bold{41}(2), 484-507.
}

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}

%\note{}
\seealso{\code{\link[MKinfer]{perm.t.test}}}

\examples{
 \dontrun{
  ###############################################################################
  ## two-sample
  ## iter = 100 to reduce computation time
  ###############################################################################
  rx <- function(n) rnorm(n, mean = 0, sd = 1) 
  ry <- function(n) rnorm(n, mean = 0.5, sd = 1) 
  sim.ssize.perm.t.test(rx = rx, ry = ry, n.max = 100, iter = 100, 
                        var.equal = TRUE)
  sim.ssize.perm.t.test(rx = rx, ry = ry, n.max = 100, iter = 100, 
                        var.equal = TRUE, parallel = TRUE)
  sim.ssize.perm.t.test(rx = rx, ry = ry, n.min = 60, n.max = 65, step.size = 1, 
                        iter = 1000, BREAK = FALSE, var.equal = TRUE,
                        parallel = TRUE)
  ## compared to
  power.t.test(delta = 0.5, power = 0.8)
  
  rx <- function(n) rnorm(n, mean = 0, sd = 1) 
  ry <- function(n) rnorm(n, mean = 0.5, sd = 1.5) 
  sim.ssize.perm.t.test(rx = rx, ry = ry, n.max = 100, iter = 1000, alternative = "less")
  sim.ssize.perm.t.test(rx = rx, ry = ry, n.min = 85, n.max = 90, step.size = 1, 
                        iter = 1000, BREAK = FALSE, alternative = "less")
  ## compared to
  power.welch.t.test(delta = 0.5, sd = 1, sd2 = 1.5, power = 0.8, alternative = "one.sided")
  
  rx <- function(n) rnorm(n, mean = 0.5, sd = 1)
  ry <- function(n) rnorm(n, mean = 0, sd = 1)
  sim.ssize.perm.t.test(rx = rx, ry = ry, n.max = 100, iter = 1000, alternative = "greater")
  sim.ssize.perm.t.test(rx = rx, ry = ry, n.min = 50, n.max = 55, step.size = 1, 
                        iter = 1000, BREAK = FALSE, alternative = "greater")
  ## compared to
  power.t.test(delta = 0.5, power = 0.8, alternative = "one.sided")
  
  rx <- function(n) rgamma(n, scale = 10, shape = 1)
  ry <- function(n) rgamma(n, scale = 15, shape = 1)
  sim.ssize.perm.t.test(rx = rx, ry = ry, n.max = 200, iter = 1000)
  sim.ssize.perm.t.test(rx = rx, ry = ry, n.min = 125, n.max = 135, step.size = 1, 
                        iter = 1000, BREAK = FALSE)
  
  ###############################################################################
  ## one-sample
  ## iter = 1000 to reduce computation time
  ###############################################################################
  rx <- function(n) rnorm(n, mean = 0.5, sd = 1)
  sim.ssize.perm.t.test(rx = rx, mu = 0, type = "one.sample", n.max = 100, iter = 1000)
  sim.ssize.perm.t.test(rx = rx, mu = 0, type = "one.sample", n.min = 33, n.max = 38, 
                        step.size = 1, iter = 1000, BREAK = FALSE)
  ## compared to
  power.t.test(delta = 0.5, power = 0.8, type = "one.sample")
  
  sim.ssize.perm.t.test(rx = rx, mu = 0, type = "one.sample", n.max = 100, iter = 1000,
                        alternative = "greater")
  sim.ssize.perm.t.test(rx = rx, mu = 0, type = "one.sample", n.min = 25, n.max = 30, 
                        step.size = 1, iter = 1000, BREAK = FALSE, alternative = "greater")
  ## compared to
  power.t.test(delta = 0.5, power = 0.8, type = "one.sample", alternative = "one.sided")
  
  sim.ssize.perm.t.test(rx = rx, mu = 1, type = "one.sample", n.max = 100, iter = 1000,
                        alternative = "less")
  sim.ssize.perm.t.test(rx = rx, mu = 1, type = "one.sample", n.min = 20, n.max = 30, 
                        step.size = 1, iter = 1000, BREAK = FALSE, alternative = "less")
  ## compared to
  power.t.test(delta = 0.5, power = 0.8, type = "one.sample", alternative = "one.sided")
  
  rx <- function(n) rgamma(n, scale = 10, shape = 1)
  sim.ssize.perm.t.test(rx = rx, mu = 5, type = "one.sample", n.max = 200, iter = 1000)
  sim.ssize.perm.t.test(rx = rx, mu = 5, type = "one.sample", n.min = 40, n.max = 50, 
                        step.size = 1, iter = 1000, BREAK = FALSE)
  
  ###############################################################################
  ## paired
  ## identical to one-sample, requires random number generating function 
  ## that simulates the difference x-y
  ## iter = 1000 to reduce computation time
  ###############################################################################
  rxy <- function(n) rnorm(n, mean = 0.5, sd = 1)
  sim.ssize.perm.t.test(rx = rxy, mu = 0, type = "paired", n.max = 100, 
                        iter = 1000)
  sim.ssize.perm.t.test(rx = rxy, mu = 0, type = "paired", n.min = 33, 
                        n.max = 38, step.size = 1, iter = 1000, BREAK = FALSE)
  ## compared to
  power.t.test(delta = 0.5, power = 0.8, type = "paired")
}
}
\keyword{htest}
