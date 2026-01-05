\name{sim.power.cond.test}
\alias{sim.power.cond.test}
\title{Monte Carlo Simulations for Empirical Power of Conditional Two-sample Test}
\description{
  Simulate the empirical power and type-I-error of the two-sample test that 
  follows pre-testing of assumptions.
}
\usage{
sim.power.cond.test(nx, rx = NULL, rx.H0 = NULL, ny, ry = NULL, ry.H0 = NULL, 
                    sig.level = 0.05, conf.int = FALSE, mu = 0, 
                    alternative = c("two.sided", "less", "greater"), 
                    iter = 10000, normality = TRUE, variance = TRUE, 
                    var.equal = FALSE, exact = TRUE)
}
\arguments{
 \item{nx}{ single numeric, sample size of first group.}
 \item{rx}{ function to simulate the values of first group (assuming H1).}
 \item{rx.H0}{ \code{NULL} or function to simulate the values of first group (assuming H0).}
 \item{ny}{ single numeric, sample size of second group.}
 \item{ry}{ function to simulate the values of second group (assuming H1).}
 \item{ry.H0}{ \code{NULL} or function to simulate the values of second group (assuming H0).}
 \item{alternative}{one- or two-sided test.  Can be abbreviated.}
 \item{sig.level}{ significance level (type I error probability)}
 \item{conf.int}{ logical, shall confidence intervals be computed. 
 Strongly increases computation time!}
 \item{mu}{ true value of the location shift for the null hypothesis.}
 \item{iter}{ single positive integer, number of interations of the simulations.}
 \item{normality}{ logical, group-wise pre-testing of the assumption of normality 
  applying the Shapiro-Wilk normality test; see \code{\link[stats]{shapiro.test}}.}
 \item{variance}{ logical, pre-testing of the assumption of equal variances 
  applying Levene's test; see \code{\link[matrixTests]{levene}}.}
 \item{var.equal}{ logical, only relevant if \code{variance} is \code{FALSE}. 
  If \code{TRUE}, Student's t-test is used, otherwise Welch's t-test.}
 \item{exact}{ logical, compute an exact p-value; see \code{\link[matrixTests]{wilcoxon}}.}
}
\details{
  If specified, functions \code{rx} and \code{ry} are used to simulate the data under the
  alternative hypothesis H1. If specified, functions \code{rx.H0} and \code{ry.H0}
  simulte the data unter the null hypothesis H0.
  
  At least one pair of functions has to be provided, either \code{rx} and \code{ry}
  or \code{rx.H0} and \code{ry.H0}.
  
  For fast computations functions from package \code{matrixTests} are used.
  
  If the test for normality is significant, the Wilcoxon-Mann-Whitney test is computed.
  
  If the test for equal variances is significant, the Welch t-test is used - otherwise,
  Student't t-test.
  
  If \code{normality} is \code{FALSE}, the Shapiro-Wilk normality test and 
  the Wilcoxon-Mann-Whitney test are skipped.
  
  If \code{variance} is \code{FALSE}, Levene's test is skipped and Student's or
  Welch t-test is computed depending on \code{var.equal}.
}
\value{
  Object of class \code{"sim.power.cond.test"} with the results of the 
  conditional two-sample tests. A list with elements \code{H1} (optional), 
  \code{H0} (optional) and \code{SetUp}.
}
\references{
  Levene, H. (1960). Robust tests for equality of variances. In Olkin, I.,
  Hotelling, H. et al. (eds.). \emph{Contributions to Probability and Statistics: 
  Essays in Honor of Harold Hotelling}. Stanford University Press. pp. 278-292.

  Mann, H and Withney, D (1947). On a test of whether one of two random variables 
  is stochastically larger than the other. \emph{Annals of mathematical Statistics},
  \bold{18}, 50-60.
  
  Shapiro, S. S., Wilk, M. B. (1965). An Analysis of Variance Test for Normality 
  (Complete Samples). \emph{Biometrika}, \bold{52}(3/4), 591-611.
  
  Student (1908). The Probable Error of a Mean. \emph{Biometrika}, \bold{6}(1): 1-25.
  
  Welch, B. L. (1947). The generalization of ``Student's'' problem when several 
  different population variances are involved. \emph{Biometrika}, \bold{34} 
  (1-2): 28-35.
  
  Wilcoxon, F (1945). Individual Comparisons by Ranking Methods. 
  \emph{Biometrics Bulletin}, \bold{1}, 80-83.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\seealso{\code{\link[stats]{wilcox.test}}, \code{\link[coin]{LocationTests}}, 
\code{\link[matrixTests]{wilcoxon}}}
\examples{
  ## Equal variance, small sample size
  power.t.test(n = 6, power = 0.8)
  (res <- sim.power.cond.test(nx = 6, rx = rnorm, rx.H0 = rnorm,
                              ny = 6, ry = function(x) rnorm(x, mean = 2), 
                              ry.H0 = rnorm))
  ## test of normality not significant?
  table(res$H1$normality)
  ## test of equal variance not significant? 
  ## (only if test for normality was not significant)
  table(res$H1$variance, useNA = "ifany")
  table(res$H0$normality)
  table(res$H0$variance, useNA = "ifany")
  
  ## only power
  (res <- sim.power.cond.test(nx = 6, rx = rnorm, 
                             ny = 6, ry = function(x) rnorm(x, mean = 2)))
  ## only type I error
  (res <- sim.power.cond.test(nx = 6, rx.H0 = rnorm, 
                              ny = 6, ry.H0 = rnorm))
                              
\dontrun{
  ## only test for equal variances
  ## assume equal variances -> Student's t-test
  (res <- sim.power.cond.test(nx = 6, rx = rnorm, rx.H0 = rnorm,
                              ny = 6, ry = function(x) rnorm(x, mean = 2), 
                              ry.H0 = rnorm, normality = FALSE, var.equal = TRUE))
  ## test of normality not significant?
  table(res$H1$normality, useNA = "ifany")
  ## test of equal variance not significant? 
  ## (only if test for normality was not significant)
  table(res$H1$variance, useNA = "ifany")
  table(res$H0$normality)
  table(res$H0$variance, useNA = "ifany")
  
  ## only test for equal variances
  ## assume unequal variances -> Welch's t-test
  (res <- sim.power.cond.test(nx = 6, rx = rnorm, rx.H0 = rnorm,
                              ny = 6, ry = function(x) rnorm(x, mean = 2), 
                              ry.H0 = rnorm, normality = FALSE, var.equal = FALSE))
  ## test of normality not significant?
  table(res$H1$normality, useNA = "ifany")
  ## test of equal variance not significant? 
  ## (only if test for normality was not significant)
  table(res$H1$variance, useNA = "ifany")
  table(res$H0$normality)
  table(res$H0$variance, useNA = "ifany")
  
  ## only test of normality
  ## assume equal variances -> Student's t-test
  (res <- sim.power.cond.test(nx = 6, rx = rnorm, rx.H0 = rnorm,
                              ny = 6, ry = function(x) rnorm(x, mean = 2), 
                              ry.H0 = rnorm, variance = FALSE, var.equal = TRUE))
  ## test of normality not significant?
  table(res$H1$normality, useNA = "ifany")
  ## test of equal variance not significant? 
  ## (only if test for normality was not significant)
  table(res$H1$variance, useNA = "ifany")
  table(res$H0$normality)
  table(res$H0$variance, useNA = "ifany")
  
  ## only test of normality
  ## assume unequal variances -> Welch's t-test
  (res <- sim.power.cond.test(nx = 6, rx = rnorm, rx.H0 = rnorm,
                              ny = 6, ry = function(x) rnorm(x, mean = 2), 
                              ry.H0 = rnorm, variance = FALSE, var.equal = FALSE))
  ## test of normality not significant?
  table(res$H1$normality, useNA = "ifany")
  ## test of equal variance not significant? 
  ## (only if test for normality was not significant)
  table(res$H1$variance, useNA = "ifany")
  table(res$H0$normality)
  table(res$H0$variance, useNA = "ifany")
}                    
}
\keyword{htest}
