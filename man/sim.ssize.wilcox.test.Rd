\name{sim.ssize.wilcox.test}
\alias{sim.ssize.wilcox.test}
\title{Sample Size for Wilcoxon Rank Sum and Signed Rank Tests}
\usage{
sim.ssize.wilcox.test(rx, ry = NULL, mu = 0, sig.level = 0.05, power = 0.8, 
                      type = c("two.sample", "one.sample", "paired"), 
                      alternative = c("two.sided", "less", "greater"),
                      n.min = 10, n.max = 200, step.size = 10, 
                      iter = 10000, BREAK = TRUE)
}
\arguments{
 \item{rx}{ function to simulate the values of x, respectively x-y in the paired case.}
 \item{ry}{ function to simulate the values of y in the two-sample case}
 \item{mu}{ true values of the location shift for the null hypothesis.}
 \item{sig.level}{ significance level (Type I error probability)}
 \item{power}{ two-sample, one-sample or paired test }
 \item{type}{ one- or two-sided test.  Can be abbreviated.}
 \item{alternative}{one- or two-sided test.  Can be abbreviated.}
 \item{n.min}{ integer, start value of grid search.}
 \item{n.max}{ integer, stop value of grid search.}
 \item{step.size}{ integer, step size used in the grid search.}
 \item{iter}{ integer, number of interations of the simulations.}
 \item{BREAK}{ logical, grid search stops when the emperical power is larger
 than the requested power.}
}
\description{
  Simulate the empirical power of Wilcoxon rank sum and signed rank tests for
  computing the required sample size.
}

\details{
  Functions \code{rx} and \code{ry} are used to simulate the data and 
  functions \code{row_wilcoxon_twosample} and \code{row_wilcoxon_onesample}
  are used to efficiently compute the p values of the respective test.
  
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
 Wilcoxon, F (1945). Individual Comparisons by Ranking Methods. 
 \emph{Biometrics Bulletin}, \bold{1}, 80-83.
}

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}

%\note{}
\seealso{\code{\link[stats]{wilcox.test}}, \code{\link[matrixTests]{wilcoxon}}}

\examples{
 \donttest{
  ###############################################################################
  ## two-sample
  ## iter = 1000 to reduce check time
  ###############################################################################
  rx <- function(n) rnorm(n, mean = 0, sd = 1) 
  ry <- function(n) rnorm(n, mean = 0.5, sd = 1) 
  sim.ssize.wilcox.test(rx = rx, ry = ry, n.max = 100, iter = 1000)
  sim.ssize.wilcox.test(rx = rx, ry = ry, n.min = 65, n.max = 70, step.size = 1, 
                        iter = 1000, BREAK = FALSE)
  ## compared to
  power.t.test(delta = 0.5, power = 0.8)
  
  rx <- function(n) rnorm(n, mean = 0, sd = 1) 
  ry <- function(n) rnorm(n, mean = 0.5, sd = 1.5) 
  sim.ssize.wilcox.test(rx = rx, ry = ry, n.max = 100, iter = 1000, alternative = "less")
  sim.ssize.wilcox.test(rx = rx, ry = ry, n.min = 85, n.max = 90, step.size = 1, 
                        iter = 1000, BREAK = FALSE, alternative = "less")
  ## compared to
  power.welch.t.test(delta = 0.5, sd = 1, sd2 = 1.5, power = 0.8, alternative = "one.sided")
  
  rx <- function(n) rnorm(n, mean = 0.5, sd = 1)
  ry <- function(n) rnorm(n, mean = 0, sd = 1)
  sim.ssize.wilcox.test(rx = rx, ry = ry, n.max = 100, iter = 1000, alternative = "greater")
  sim.ssize.wilcox.test(rx = rx, ry = ry, n.min = 50, n.max = 55, step.size = 1, 
                        iter = 1000, BREAK = FALSE, alternative = "greater")
  ## compared to
  power.t.test(delta = 0.5, power = 0.8, alternative = "one.sided")
  
  rx <- function(n) rgamma(n, scale = 10, shape = 1)
  ry <- function(n) rgamma(n, scale = 15, shape = 1)
  sim.ssize.wilcox.test(rx = rx, ry = ry, n.max = 200, iter = 1000)
  sim.ssize.wilcox.test(rx = rx, ry = ry, n.min = 125, n.max = 135, step.size = 1, 
                        iter = 1000, BREAK = FALSE)
  
  ###############################################################################
  ## one-sample
  ## iter = 1000 to reduce check time
  ###############################################################################
  rx <- function(n) rnorm(n, mean = 0.5, sd = 1)
  sim.ssize.wilcox.test(rx = rx, mu = 0, type = "one.sample", n.max = 100, iter = 1000)
  sim.ssize.wilcox.test(rx = rx, mu = 0, type = "one.sample", n.min = 33, n.max = 38, 
                        step.size = 1, iter = 1000, BREAK = FALSE)
  ## compared to
  power.t.test(delta = 0.5, power = 0.8, type = "one.sample")
  
  sim.ssize.wilcox.test(rx = rx, mu = 0, type = "one.sample", n.max = 100, iter = 1000,
                        alternative = "greater")
  sim.ssize.wilcox.test(rx = rx, mu = 0, type = "one.sample", n.min = 25, n.max = 30, 
                        step.size = 1, iter = 1000, BREAK = FALSE, alternative = "greater")
  ## compared to
  power.t.test(delta = 0.5, power = 0.8, type = "one.sample", alternative = "one.sided")
  
  sim.ssize.wilcox.test(rx = rx, mu = 1, type = "one.sample", n.max = 100, iter = 1000,
                        alternative = "less")
  sim.ssize.wilcox.test(rx = rx, mu = 1, type = "one.sample", n.min = 20, n.max = 30, 
                        step.size = 1, iter = 1000, BREAK = FALSE, alternative = "less")
  ## compared to
  power.t.test(delta = 0.5, power = 0.8, type = "one.sample", alternative = "one.sided")
  
  rx <- function(n) rgamma(n, scale = 10, shape = 1)
  sim.ssize.wilcox.test(rx = rx, mu = 5, type = "one.sample", n.max = 200, iter = 1000)
  sim.ssize.wilcox.test(rx = rx, mu = 5, type = "one.sample", n.min = 40, n.max = 50, 
                        step.size = 1, iter = 1000, BREAK = FALSE)
  
  ###############################################################################
  ## paired
  ## identical to one-sample, requires random number generating function 
  ## that simulates the difference x-y
  ## iter = 1000 to reduce check time
  ###############################################################################
  rxy <- function(n) rnorm(n, mean = 0.5, sd = 1)
  sim.ssize.wilcox.test(rx = rxy, mu = 0, type = "paired", n.max = 100, 
                        iter = 1000)
  sim.ssize.wilcox.test(rx = rxy, mu = 0, type = "paired", n.min = 33, 
                        n.max = 38, step.size = 1, iter = 1000, BREAK = FALSE)
  ## compared to
  power.t.test(delta = 0.5, power = 0.8, type = "paired")
}
}
\keyword{htest}
