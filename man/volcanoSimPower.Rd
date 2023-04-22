\name{volcano}
\alias{volcano}
\alias{volcano.sim.power.ttest}
\alias{volcano.sim.power.wtest}
\title{ Volcano Plots }
\description{
  Produce volcano plot(s) for simulations of power and type-I-error of tests. 
}
\usage{
\S3method{volcano}{sim.power.ttest}(x, alpha = 1, shape = 19, 
                                    hex = FALSE, bins = 50, \dots)

\S3method{volcano}{sim.power.wtest}(x, alpha = 1, shape = 19, 
                                    hex = FALSE, bins = 50, \dots)
}
\arguments{
  \item{x}{ object of class \code{sim.power.ttest}.}
  \item{alpha}{ bleding factor (default: no blending.}
  \item{shape}{ point shape used.}
  \item{hex}{ logical, should hexagonal binning be used.}
  \item{bins}{ number of bins used for hexagonal binning.}
  \item{\dots}{ further arguments that may be passed through).}
}
\details{
  The plot generates a \code{ggplot2} object that is shown. 

  Missing values are handled by the \code{ggplot2} functions.
}
\value{
  Object of class \code{gg} and \code{ggplot}.
}
\references{
  Wikipedia contributors, \emph{Volcano plot (statistics)}, Wikipedia, The Free Encyclopedia, 
  \url{https://en.wikipedia.org/w/index.php?title=Volcano_plot_(statistics)&oldid=900217316} 
  (accessed December 25, 2019). 

  For more sophisticated and flexible volcano plots see for instance:
  Blighe K, Rana S, Lewis M (2019). EnhancedVolcano: Publication-ready volcano 
  plots with enhanced colouring and labeling. R/Bioconductor package.
  \url{https://github.com/kevinblighe/EnhancedVolcano}.
}
\author{ Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link[MKinfer]{volcano}}}
\examples{
res1 <- sim.power.t.test(nx = 5, rx = rnorm, rx.H0 = rnorm, 
                        ny = 10, ry = function(x) rnorm(x, mean = 3, sd = 3), 
                        ry.H0 = function(x) rnorm(x, sd = 3))
volcano(res1)
\donttest{
## low number of iterations to reduce computation time
res2 <- sim.power.wilcox.test(nx = 6, rx = rnorm, rx.H0 = rnorm,
                      ny = 6, ry = function(x) rnorm(x, mean = 2), 
                      ry.H0 = rnorm, iter = 100, conf.int = TRUE)
volcano(res2)
}
}
\keyword{hplot}
