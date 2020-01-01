\name{sim.power.wilcox.test}
\alias{sim.power.wilcox.test}
\title{Monte Carlo Simulations for Empirical Power of Wilcoxon-Mann-Whitney Tests}
\description{
  Simulate the empirical power and type-I-error of Wilcoxon-Mann-Whitney tests.
}
\usage{
sim.power.wilcox.test(nx, rx, rx.H0 = NULL, ny, ry, ry.H0 = NULL, 
                      alternative = c("two.sided", "less", "greater"), 
                      sig.level = 0.05, conf.int = FALSE, approximate = FALSE,
                      ties = FALSE, iter = 10000, nresample = 10000,
                      parallel = "no", ncpus = 1L, cl = NULL)
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
 \item{approximate}{ logical, shall an approximate test be computed; 
 see \code{\link[coin]{LocationTests}}. Increases computation time!}
 \item{ties}{ logical, indicating whether ties may occur. Increases computation time!}
 \item{iter}{ single positive integer, number of interations of the simulations.}
 \item{nresample}{ single positive integer, the number of Monte Carlo replicates 
 used for the computation of the approximative reference distribution; 
 see \code{\link[coin]{NullDistribution}}.}
 \item{parallel}{ a character, the type of parallel operation: either \code{"no"} 
 (default), \code{"multicore"} or \code{"snow"}; see \code{\link[coin]{NullDistribution}}.}
 \item{ncpus}{ a single integer, the number of processes to be used in parallel operation. 
 Defaults to 1L; see \code{\link[coin]{NullDistribution}}.}
 \item{cl}{ an object inheriting from class \code{"cluster"}, specifying an optional 
 parallel or snow cluster if \code{parallel = "snow"}. Defaults to NULL; see 
 \code{\link[coin]{NullDistribution}}.}
}
\details{
  Functions \code{rx} and \code{ry} are used to simulate the data under the
  alternative hypothesis H1. If specified, functions \code{rx.H0} and \code{ry.H0}
  simulte the data unter the null hypothesis H0.
  
  For fast computations functions from package \code{matrixTests} and package
  \code{coin} are used.
}
\value{
  Object of class \code{"sim.power.wtest"} with the results of the 
  Wilcoxon-Mann-Whitney tests. A list elements \code{Exact}, \code{Asymptotic} 
  and \code{Approximate}. In addition, the simulation setup is saved in 
  element \code{SetUp}.
}
\references{
  Mann, H and Withney, D (1947). On a test of whether one of two random variables 
  is stochastically larger than the other. \emph{Annals of mathematical Statistics},
  \bold{18}, 50-60.
  
  Wilcoxon, F (1945). Individual Comparisons by Ranking Methods. 
  \emph{Biometrics Bulletin}, \bold{1}, 80-83.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\seealso{\code{\link[stats]{wilcox.test}}, \code{\link[coin]{LocationTests}}, 
\code{\link[matrixTests]{wilcoxon}}}
\examples{
## Equal variance, small sample size
power.t.test(n = 5, power = 0.8)
sim.ssize.wilcox.test(rx = rnorm, ry = function(x) rnorm(x, mean = 2), 
                      power = 0.8, n.min = 3, n.max = 10, step.size = 1)
sim.power.wilcox.test(nx = 6, rx = rnorm, rx.H0 = rnorm,
                      ny = 6, ry = function(x) rnorm(x, mean = 2), 
                      ry.H0 = rnorm)
}
\keyword{htest}
