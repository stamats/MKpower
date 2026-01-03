\name{sim.power.t.test}
\alias{sim.power.t.test}
\title{Monte Carlo Simulations for Empirical Power of Two-sample t-Tests}
\description{
  Simulate the empirical power and type-I-error of two-sample t-tests; i.e.,
  Student (equal variances), Welch, Hsu, Xiao, permutation and boostrap t-tests.
}
\usage{
sim.power.t.test(nx, rx, rx.H0 = NULL, ny, ry, ry.H0 = NULL, 
                 sig.level = 0.05, conf.int = FALSE, mu = 0, 
                 alternative = c("two.sided", "less", "greater"), 
                 methods = c("student", "welch", "hsu", "xiao"),
                 iter = 10000, R = 9999, useCombn = FALSE, 
                 parallel = TRUE, cl = NULL)
}
\arguments{
 \item{nx}{ single numeric, sample size of first group.}
 \item{rx}{ function to simulate the values of first group (assuming H1).}
 \item{rx.H0}{ \code{NULL} or function to simulate the values of first group (assuming H0).}
 \item{ny}{ single numeric, sample size of second group.}
 \item{ry}{ function to simulate the values of second group (assuming H1).}
 \item{ry.H0}{ \code{NULL} or function to simulate the values of second group (assuming H0).}
 \item{sig.level}{ significance level (type I error probability)}
 \item{conf.int}{ logical, shall confidence intervals be computed. 
   Increases computation time!}
 \item{mu}{ true value of the location shift for the null hypothesis.}
 \item{alternative}{one- or two-sided test.  Can be abbreviated.}
 \item{methods}{ character, one or more of the following names: \code{"student"},
   \code{"welch"}, \code{"hsu"}, \code{"xiao"}, \code{"perm.stud"}, 
   \code{"perm.welch"}, \code{"boot.stud"}, \code{"boot.welch"}.}
 \item{iter}{ single integer, number of interations of the simulations.}
 \item{R}{ number of (Monte-Carlo) permutations for permutation and bootstrap tests.}
 \item{useCombn}{logical; use all combinations instead of permutations in the
    case of the two-sample t-test; see \code{\link[MKinfer:perm.t.test]{perm.t.test}}.}
  \item{parallel}{logical; use parellel computing via function \code{\link{parRapply}}
    in the case of method \code{"xiao"}.}
  \item{cl}{a cluster object, created by package \code{\link[parallel:makeCluster]{parallel}} 
    or by package \code{\link[snow:makeCluster]{snow}}. If \code{NULL}, a new cluster 
    will be generated with function \code{\link[parallel:makeCluster]{makeCluster}} 
    using the maximum number of available cores - 1.}
}
\details{
  Functions \code{rx} and \code{ry} are used to simulate the data under the
  alternative hypothesis H1. If specified, functions \code{rx.H0} and \code{ry.H0}
  simulte the data unter the null hypothesis H0.
  
  For fast computations functions from package \code{matrixTests} are used.
}
\value{
  Object of class \code{"sim.power.ttest"} with the results of the selected t-tests
  in the list elements \code{Student}, \code{Welch}, \code{Hsu}, \code{Xiao}, 
  \code{Perm}, \code{Boot}. In addition, the simulation setup is saved in 
  element \code{SetUp}.
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
  
  Y. Xiao (2018). On the Solution of a Generalized Behrens-Fisher Problem. 
  \emph{Far East Journal of Theoretical Statistics}, \bold{54} (1), 21-140.
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
                   ny = 5, ry = function(n) rnorm(n, mean = 2), 
                   ry.H0 = rnorm, methods = c("student", "welch"))
  
  ## Unequal variance, small sample size
  power.welch.t.test(n = 5, delta = 5, sd1 = 1, sd2 = 3)
  power.hsu.t.test(n = 5, delta = 5, sd1 = 1, sd2 = 3)
  sim.power.t.test(nx = 5, rx = rnorm, rx.H0 = rnorm, 
                   ny = 5, ry = function(n) rnorm(n, mean = 5, sd = 3), 
                   ry.H0 = function(n) rnorm(n, sd = 3), 
                   methods = c("student", "welch"))
                   
  ## Unequal variance, small unequal sample sizes
  power.welch.t.test2(n = c(5, 10), delta = 3, sd1 = 1, sd2 = 3)
  sim.power.t.test(nx = 5, rx = rnorm, rx.H0 = rnorm, 
                   ny = 10, ry = function(n) rnorm(n, mean = 3, sd = 3), 
                   ry.H0 = function(n) rnorm(n, sd = 3), 
                   methods = c("student", "welch"))
  \dontrun{
    ## with parallel computing in case of xiao
    sim.power.t.test(nx = 5, rx = rnorm, rx.H0 = rnorm, 
                     ny = 10, ry = function(n) rnorm(n, mean = 3, sd = 3), 
                     ry.H0 = function(n) rnorm(n, sd = 3), 
                     methods = c("student", "welch", "hsu", "xiao"))
  }
}
\keyword{htest}
