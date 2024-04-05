\name{print.power.mpe.test}
\alias{print.power.mpe.test}
\title{Print Methods for Hypothesis Tests, Sample size and Power Calculations}
\description{
  Printing objects of class \code{"power.mpe.test"} by simple \code{\link{print}} methods.
}
\usage{
  \method{print}{power.mpe.test}(x, digits = getOption("digits"), ...)
}
\arguments{
  \item{x}{object of class \code{"power.mpe.test"}.}
  \item{digits}{number of significant digits to be used.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  The \code{print} method is based on the respective method 
  \code{print.power.htest} of package \pkg{stats}.

  A \code{power.mpe.test} object is just a named list of numbers and
  character strings, supplemented with \code{method} and \code{note}
  elements.  The \code{method} is displayed as a title, the \code{note}
  as a footnote, and the remaining elements are given in an aligned
  \sQuote{name = value} format.
}
\value{
  the argument \code{x}, invisibly, as for all \code{\link{print}}
  methods.
}
\author{Srinath Kolampally, Matthias Kohl \email{Matthias.Kohl@stamats.de}}
\note{The function first appeared in package \pkg{mpe}, which is now archived on CRAN.}
\seealso{\code{\link[stats]{print.power.htest}}, \code{\link{power.mpe.known.var}}, 
\code{\link{power.mpe.unknown.var}}
}
\examples{
(pkv <- power.mpe.known.var(K = 2, delta = c(1,1), Sigma = diag(c(2,2)), power = 0.9,
                            sig.level = 0.025))
print(pkv, digits =  4) # using less digits than default
print(pkv, digits = 12) # using more digits than default
}
\keyword{htest}
\keyword{power.htest}
