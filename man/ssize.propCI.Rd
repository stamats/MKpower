\name{ssize.propCI}
\alias{ssize.propCI}
\encoding{UTF-8}
\title{Sample Size Calculation for Confidence Interval of a Proportion}
\description{
  Compute the sample size for the two-sided confidence interval of a single proportion.
}
\usage{
ssize.propCI(prop, width, conf.level = 0.95,  method = "wald-cc")
}
\arguments{
  \item{prop}{expected proportion}
  \item{width}{width of the confidence interval}
  \item{conf.level}{confidence level}
  \item{method}{ method used to compute the confidence interval; see Details.}
}
\details{
  The computation is based on the asymptotic formulas provided in Section 2.5.2
  of Fleiss et al. (2003). If \code{method = "wald-cc"} a continuity correction 
  is applied.
}
\value{
  Object of class \code{"power.htest"}, a list of the arguments
  (including the computed one) augmented with \code{method} and
  \code{note} elements.
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}

\references{
  J.L. Fleiss, B. Levin and M.C. Paik (2003). 
  \emph{Statistical Methods for Rates and Proportions}. 
  Wiley Series in Probability and Statistics.
}

\seealso{\code{\link{power.prop1.test}}, \code{\link[MKinfer]{binomCI}}}

\examples{
ssize.propCI(prop = 0.3, width = 0.1)
ssize.propCI(prop = 0.1, width = 0.1)
}
\keyword{htest}
