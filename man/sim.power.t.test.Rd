\name{sim.power.t.test}
\alias{sim.power.t.test}
\title{Monte Carlo Simulations for Empirical Power of Two-sample t-Tests}
\description{
  Simulate the empirical power and type-I-error of two-sample t-tests; i.e.,
  classical (equal variances), Welch and Hsu t-tests.
}
\usage{
sim.power.t.test(nx, rx, rx.H0 = NULL, ny, ry, ry.H0 = NULL, 
                 sig.level = 0.05, mu = 0, 
                 alternative = c("two.sided", "less", "greater"), 
                 iter = 10000)
}
\arguments{
 \item{nx}{ single numeric, sample size of first group.}
 \item{rx}{ function to simulate the values of first group (assuming H1).}
 \item{rx.H0}{ \code{NULL} or function to simulate the values of first group (assuming H0).}
 \item{ny}{ single numeric, sample size of second group.}
 \item{ry}{ function to simulate the values of second group (assuming H1).}
 \item{ry.H0}{ \code{NULL} or function to simulate the values of second group (assuming H0).}
 \item{sig.level}{ significance level (type I error probability)}
 \item{mu}{ true value of the location shift for the null hypothesis.}
 \item{alternative}{one- or two-sided test.  Can be abbreviated.}
 \item{iter}{ single integer, number of interations of the simulations.}
}
\details{
  Functions \code{rx} and \code{ry} are used to simulate the data under the
  alternative hypothesis H1. If specified, functions \code{rx.H0} and \code{ry.H0}
  simulte the data unter the null hypothesis H0.
  
  For fast computations functions from package \code{matrixTests} are used.
}
\value{
  Object of class \code{"sim.power.ttest"} with the results of the three t-tests
  in the list elements \code{Classical}, \code{Welch} and \code{Hsu}. In addition,
  the simulation setup is saved in element \code{SetUp}.
}
\references{
  J. Hedderich, L. Sachs. \emph{Angewandte Statistik: Methodensammlung mit R}.
  Springer 2018.
  
  Hsu, P. (1938). Contribution to the theory of ``student's'' t-test as 
  applied to the problem of two samples. \emph{Statistical Research Memoirs} \bold{2}: 
  1-24.
  
  Student (1908). The Probable Error of a Mean. \emph{Biometrika}, \bold{6}(1): 1-25.
  
  Welch, B. L. (1947). The generalization of ``Student's'' problem when several 
  different population variances are involved. \emph{Biometrika}, \bold{34} 
  (1-2): 28-35.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\seealso{\code{\link[stats]{t.test}}, \code{\link[MKinfer]{hsu.t.test}}, 
\code{\link[matrixTests]{ttest}}}
\examples{
## Equal variance, small sample size
power.t.test(n = 5, delta = 2)
power.welch.t.test(n = 5, delta = 2)
power.hsu.t.test(n = 5, delta = 2)
sim.power.t.test(nx = 5, rx = rnorm, rx.H0 = rnorm,
                 ny = 5, ry = function(x) rnorm(x, mean = 2), ry.H0 = rnorm)

## Equal variance, moderate sample size
power.t.test(n = 25, delta = 0.8)
power.welch.t.test(n = 25, delta = 0.8)
power.hsu.t.test(n = 25, delta = 0.8)
sim.power.t.test(nx = 25, rx = rnorm, rx.H0 = rnorm,
                 ny = 25, ry = function(x) rnorm(x, mean = 0.8), ry.H0 = rnorm)

## Equal variance, high sample size
power.t.test(n = 100, delta = 0.4)
power.welch.t.test(n = 100, delta = 0.4)
power.hsu.t.test(n = 100, delta = 0.4)
sim.power.t.test(nx = 100, rx = rnorm, rx.H0 = rnorm,
                 ny = 100, ry = function(x) rnorm(x, mean = 0.4), ry.H0 = rnorm)

## Unequal variance, small sample size
power.welch.t.test(n = 5, delta = 5, sd1 = 1, sd2 = 3)
power.hsu.t.test(n = 5, delta = 5, sd1 = 1, sd2 = 3)
sim.power.t.test(nx = 5, rx = rnorm, rx.H0 = rnorm, 
                 ny = 5, ry = function(x) rnorm(x, mean = 5, sd = 3), 
                 ry.H0 = function(x) rnorm(x, sd = 3))
                 
## Unequal variance, moderate sample size
power.welch.t.test(n = 25, delta = 1.8, sd1 = 1, sd2 = 3)
power.hsu.t.test(n = 25, delta = 1.8, sd1 = 1, sd2 = 3)
sim.power.t.test(nx = 25, rx = rnorm, rx.H0 = rnorm, 
                 ny = 25, ry = function(x) rnorm(x, mean = 1.8, sd = 3), 
                 ry.H0 = function(x) rnorm(x, sd = 3))
                 
## Unequal variance, high sample size
power.welch.t.test(n = 100, delta = 0.9, sd1 = 1, sd2 = 3)
power.hsu.t.test(n = 100, delta = 0.9, sd1 = 1, sd2 = 3)
sim.power.t.test(nx = 100, rx = rnorm, rx.H0 = rnorm, 
                 ny = 100, ry = function(x) rnorm(x, mean = 0.9, sd = 3), 
                 ry.H0 = function(x) rnorm(x, sd = 3))
                 
## Unequal variance, unequal sample sizes
## small sample sizes
sim.power.t.test(nx = 10, rx = rnorm, rx.H0 = rnorm, 
                 ny = 5, ry = function(x) rnorm(x, mean = 5, sd = 3), 
                 ry.H0 = function(x) rnorm(x, sd = 3))
sim.power.t.test(nx = 5, rx = rnorm, rx.H0 = rnorm, 
                 ny = 10, ry = function(x) rnorm(x, mean = 3, sd = 3), 
                 ry.H0 = function(x) rnorm(x, sd = 3))

## Unequal variance, unequal sample sizes
## moderate sample sizes
sim.power.t.test(nx = 25, rx = rnorm, rx.H0 = rnorm, 
                 ny = 50, ry = function(x) rnorm(x, mean = 1.5, sd = 3), 
                 ry.H0 = function(x) rnorm(x, sd = 3))

## Unequal variance, unequal sample sizes
## high sample sizes
sim.power.t.test(nx = 100, rx = rnorm, rx.H0 = rnorm, 
                 ny = 200, ry = function(x) rnorm(x, mean = 0.6, sd = 3), 
                 ry.H0 = function(x) rnorm(x, sd = 3))
}
\keyword{htest}
