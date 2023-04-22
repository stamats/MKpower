\name{qqunif}
\alias{qqunif}
\alias{qqunif.default}
\alias{qqunif.sim.power.ttest}
\alias{qqunif.sim.power.wtest}
\title{ qq-Plots for Uniform Distribution}
\description{
  Produce qq-plot(s) of the given effect size and p values assuming a uniform distribution. 
}
\usage{
qqunif(x, \dots)

\S3method{qqunif}{default}(x, min = 0, max = 1, \dots)

\S3method{qqunif}{sim.power.ttest}(x, color.line = "orange", shape = 19, size = 1, 
                           alpha = 1, \dots)

\S3method{qqunif}{sim.power.wtest}(x, color.line = "orange", shape = 19, size = 1, 
                           alpha = 1, \dots)
}
\arguments{
  \item{x}{ numeric vector or data (object).}
  \item{min}{ single numeric, lower limit of the distribution.}
  \item{max}{ single numeric, upper limit of the distribution.}
  \item{color.line}{ color of the line indicating the uniform distribution.}
  \item{shape}{ point shape.}
  \item{size}{ point size.}
  \item{alpha}{ bleding factor (default: no blending.}
  \item{\dots}{ further arguments that may be passed through).}
}
\details{
  The plot generates a \code{ggplot2} object that is shown. 

  Missing values are handled by the \code{ggplot2} functions.
}
\value{
  Object of class \code{gg} and \code{ggplot}.
}
%\references{}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
%\seealso{}
\examples{
## default
qqunif(runif(100))

## visualization of empirical power and type-I-error
res1 <- sim.power.t.test(nx = 5, rx = rnorm, rx.H0 = rnorm, 
                        ny = 10, ry = function(x) rnorm(x, mean = 3, sd = 3), 
                        ry.H0 = function(x) rnorm(x, sd = 3))
qqunif(res1, alpha = 0.1)

res2 <- sim.power.wilcox.test(nx = 6, rx = rnorm, rx.H0 = rnorm,
                      ny = 6, ry = function(x) rnorm(x, mean = 2), 
                      ry.H0 = rnorm)
qqunif(res2)
}
\keyword{hplot}
