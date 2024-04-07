\name{power.mpe.unknown.var}
\alias{power.mpe.unknown.var}
\title{Multiple Co-Primary Endpoints with Unknown Covariance}
\description{
  The function calculates either sample size or power for continuous multiple
  co-primary endpoints with unknown covariance.
}
\usage{
  power.mpe.unknown.var(K, n = NULL, delta = NULL, Sigma, SD, rho, sig.level = 0.05,
                        power = NULL, M = 10000, n.min = NULL, n.max = NULL,
                        tol = .Machine$double.eps^0.25, use.uniroot = TRUE)
}
\arguments{
  \item{K}{number of co-primary endpoints}
  \item{n}{optional: sample size}
  \item{delta}{expected effect size (length \code{K})}
  \item{Sigma}{unknown covariance matrix (dimension \code{K} x \code{K})}
  \item{SD}{unknown standard deviations (length \code{K})}
  \item{rho}{unknown correlations (length \code{0.5*K*(K-1)})}
  \item{sig.level}{significance level (Type I error probability)}
  \item{power}{optional: power of test (1 minus Type II error probability)}
  \item{M}{Number of replications for the required simulations.}
  \item{n.min}{Starting point of search interval for sample size}
  \item{n.max}{End point of search interval for sample size, must be larger than \code{n.min}}
  \item{tol}{the desired accuracy for \code{\link[stats]{uniroot}}}
  \item{use.uniroot}{Finds one root of one equation}
}
\details{
  The function can be used to either compute sample size or power for continuous
  multiple co-primary endpoints with unknown covariance. The implementation is
  based on the formulas given in the references below.

  The null hypothesis reads \eqn{\mu_{Tk}-\mu_{Ck}\le 0}{mu_Tk-mu_Ck <= 0} for
  at least one \eqn{k\in\{1,\ldots,K\}}{k in {1,...,K}} where Tk is treatment k,
  Ck is control k and K is the number of co-primary endpoints.

  One has to specify either \code{n} or \code{power}, the other parameter is
  determined. An approach to calculate sample size \code{n}, is to first call
  \code{\link{power.mpe.known.var}} and use the result as \code{n.min}. The input for
  \code{n.max} must be larger then \code{n.min}. Moreover, either covariance
  matrix \code{Sigma} or standard deviations \code{SD} and correlations \code{rho}
  must be given.

  The sample size is calculated by simulating Wishart distributed random matrices,
  hence the results include a certain random variation.
}
\value{
  Object of class \code{power.mpe.test}, a list of arguments (including the
  computed one) augmented with method and note elements.
}
\references{
  Sugimoto, T. and Sozu, T. and Hamasaki, T. (2012). A convenient formula for sample
  size calculations in clinical trials with multiple co-primary continuous endpoints.
  \emph{Pharmaceut. Statist.}, \bold{11}: 118-128. doi:10.1002/pst.505

  Sozu, T. and Sugimoto, T. and Hamasaki, T. and Evans, S.R. (2015). \emph{Sample
  Size Determination in Clinical Trials with Multiple Endpoints}. Springer Briefs in
  Statistics, ISBN 978-3-319-22005-5.
}
\author{Srinath Kolampally, Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\note{The function first appeared in package \pkg{mpe}, which is now archived on CRAN.}
\seealso{ \code{\link{power.mpe.known.var}} }
\examples{
## compute power
\dontrun{
power.mpe.unknown.var(K = 2, n = 20, delta = c(1,1), Sigma = diag(c(1,1)))

## To compute sample size, first assume covariance as known
power.mpe.known.var(K = 2, delta = c(1,1), Sigma = diag(c(2,2)), power = 0.9,
                  sig.level = 0.025)

## The value of n, which is 51, is used as n.min and n.max must be larger
## then n.min so we try 60.
power.mpe.unknown.var(K = 2, delta = c(1,1), Sigma = diag(c(2,2)), power = 0.9,
                  sig.level = 0.025, n.min = 51, n.max = 60)

## More complex example with unknown covariance matrix assumed to be
Sigma <- matrix(c(1.440, 0.840, 1.296, 0.840,
                  0.840, 1.960, 0.168, 1.568,
                  1.296, 0.168, 1.440, 0.420,
                  0.840, 1.568, 0.420, 1.960), ncol = 4)
## compute power
power.mpe.unknown.var(K = 4, n = 90, delta = c(0.5, 0.75, 0.5, 0.75), Sigma = Sigma)
## equivalent: unknown SDs and correlation rho
power.mpe.unknown.var(K = 4, n = 90, delta = c(0.5, 0.75, 0.5, 0.75),
                      SD = c(1.2, 1.4, 1.2, 1.4),
                      rho = c(0.5, 0.9, 0.5, 0.1, 0.8, 0.25))
}}
\keyword{htest}
\keyword{multivariate}
