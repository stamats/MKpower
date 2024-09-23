\name{ssize.propCI}
\alias{ssize.propCI}
\alias{ssize.prop.ci}
\encoding{UTF-8}
\title{Sample Size Calculation for Confidence Interval of a Proportion}
\description{
  Compute the sample size for the two-sided confidence interval of a single proportion.
}
\usage{
ssize.prop.ci(prop, width, conf.level = 0.95,  method = "wald-cc")

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
  
  There are also methods for Jeffreys, Clopper-Pearson, Wilson and the 
  Agresti-Coull interval; see also \code{\link[MKinfer]{binomCI}}.
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
  
  W.W. Piegorsch (2004). Sample sizes for improved binomial confidence intervals. 
  \emph{Computational Statistics & Data Analysis}, \bold{46}, 309-316.
  
  M. Thulin (2014). The cost of using exact confidence intervals for a binomial 
  proportion. \emph{Electronic Journal of Statistics}, \bold{8}(1), 817-840.
}

\seealso{\code{\link{power.prop1.test}}, \code{\link[MKinfer]{binomCI}}}

\examples{
ssize.propCI(prop = 0.1, width = 0.1)
ssize.propCI(prop = 0.3, width = 0.1)
ssize.propCI(prop = 0.3, width = 0.1, method = "wald")
ssize.propCI(prop = 0.3, width = 0.1, method = "jeffreys")
ssize.propCI(prop = 0.3, width = 0.1, method = "clopper-pearson")
ssize.propCI(prop = 0.3, width = 0.1, method = "wilson")
ssize.propCI(prop = 0.3, width = 0.1, method = "agresti-coull")
}
\keyword{htest}
