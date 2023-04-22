\name{hist}
\alias{hist}
\alias{hist.sim.power.ttest}
\alias{hist.sim.power.wtest}
\title{Histograms}
\description{Produce histograms for simulations of power and type-I-error of tests.}
\usage{
\S3method{hist}{sim.power.ttest}(x, color.hline = "orange", \dots)

\S3method{hist}{sim.power.wtest}(x, color.hline = "orange", \dots)
}
\arguments{
  \item{x}{ object of class \code{sim.power.ttest}.}
  \item{color.hline}{ color of horizontal line indicating uniform distribution of p values.}
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
\seealso{\code{\link[graphics]{hist}}}
\examples{
res1 <- sim.power.t.test(nx = 5, rx = rnorm, rx.H0 = rnorm, 
                        ny = 10, ry = function(x) rnorm(x, mean = 3, sd = 3), 
                        ry.H0 = function(x) rnorm(x, sd = 3))
hist(res1)
res2 <- sim.power.wilcox.test(nx = 6, rx = rnorm, rx.H0 = rnorm,
                      ny = 6, ry = function(x) rnorm(x, mean = 2), 
                      ry.H0 = rnorm)
hist(res2)
}
\keyword{hplot}
