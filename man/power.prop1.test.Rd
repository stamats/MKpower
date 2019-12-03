\name{power.prop1.test}
\alias{power.prop1.test}
\encoding{UTF-8}
\title{Power Calculations for One-Sample Test for Proportions}
\description{
  Compute the power of the one-sample test for proportions, or determine
  parameters to obtain a target power.
}
\usage{
power.prop1.test(n = NULL, p1 = NULL, p0 = 0.5, sig.level = 0.05, 
                             power = NULL, 
                             alternative = c("two.sided", "less", "greater"),
                             cont.corr = TRUE, tol = .Machine$double.eps^0.25)
}
\arguments{
  \item{n}{number of observations (per group)}
  \item{p1}{expected probability}
  \item{p0}{probability under the null hypothesis}
  \item{sig.level}{significance level (Type I error probability)}
  \item{power}{power of test (1 minus Type II error probability)}
  \item{alternative}{one- or two-sided test.  Can be abbreviated.}
  \item{cont.corr}{use continuity correction}
  \item{tol}{numerical tolerance used in root finding, the default
    providing (at least) four significant digits.}
}
\details{
  Exactly one of the parameters \code{n}, \code{p1}, \code{power}, and 
  \code{sig.level} must be passed as NULL, and that parameter is determined 
  from the others.  Notice that \code{sig.level} has a non-NULL default 
  so \code{NULL} must be explicitly passed if you want it computed.

  The computation is based on the asymptotic formulas provided in Section 2.5.1
  of Fleiss et al. (2003). If \code{cont.corr = TRUE} a continuity correction 
  is applied, which may lead to better approximations of the finite-sample
  values.
}
\value{
  Object of class \code{"power.htest"}, a list of the arguments
  (including the computed one) augmented with \code{method} and
  \code{note} elements.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\note{
  The documentation was adapted from \code{\link{power.prop.test}}.
}

\references{
  J.L. Fleiss, B. Levin and M.C. Paik (2003). 
  \emph{Statistical Methods for Rates and Proportions}. 
  Wiley Series in Probability and Statistics.
}

\seealso{\code{\link{power.prop.test}}, \code{\link{prop.test}}}

\examples{
power.prop1.test(p1 = 0.4, power = 0.8)
power.prop1.test(p1 = 0.4, power = 0.8, cont.corr = FALSE)
power.prop1.test(p1 = 0.6, power = 0.8)
power.prop1.test(n = 204, power = 0.8)
power.prop1.test(n = 204, p1 = 0.4, power = 0.8, sig.level = NULL)
power.prop1.test(n = 194, p1 = 0.4, power = 0.8, sig.level = NULL, 
                 cont.corr = FALSE)

power.prop1.test(p1 = 0.1, p0 = 0.3, power = 0.8, alternative = "less")
power.prop1.test(p1 = 0.1, p0 = 0.3, power = 0.8, alternative = "less", 
                 cont.corr = FALSE)
power.prop1.test(n = 31, p0 = 0.3, power = 0.8, alternative = "less")
power.prop1.test(n = 31, p1 = 0.1, p0 = 0.3, power = 0.8, sig.level = NULL, 
                 alternative = "less")


power.prop1.test(p1 = 0.5, p0 = 0.3, power = 0.8, alternative = "greater")
power.prop1.test(p1 = 0.5, p0 = 0.3, power = 0.8, alternative = "greater", 
                 cont.corr = FALSE)
power.prop1.test(n = 40, p0 = 0.3, power = 0.8, alternative = "greater")
power.prop1.test(n = 40, p1 = 0.5, p0 = 0.3, power = 0.8, sig.level = NULL, 
                 alternative = "greater")
}
\keyword{ htest }
