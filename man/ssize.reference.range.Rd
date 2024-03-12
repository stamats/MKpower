\name{ssize.reference.range}
\alias{ssize.reference.range}
\title{Power Calculations for Two-sample Hsu t Test}
\usage{
ssize.reference.range(n = NULL, delta = NULL, ref.prob = 0.95, conf.prob = NULL,
                      alternative = c("two.sided", "one.sided"),
                      method = "parametric", exact = TRUE, 
                      tol = .Machine$double.eps^0.5)
}
\arguments{
 \item{n}{number of observations}
 \item{delta}{ difference between empirical and target coverage of reference range}
 \item{ref.prob}{target coverage of reference range}
 \item{conf.prob}{confidence probability to acchieve given difference between empirical and target coverage}
 \item{alternative}{a character string specifying \code{"two.sided"} (default),
    or one-sided reference ranges. You can specify just the initial letter.}
 \item{method}{either \code{"parametric"} or \code{"nonparametric"}; see details}
 \item{exact}{use exact or approximate method}
 \item{tol}{numerical tolerance used in root finding, the default
    providing (at least) eight significant digits.}
}
\description{
  Compute the sample size for reference range studies, or determine parameters
  for a given sample size; see Jennen-Steinmetz and Wellek (2005).
}

\details{
  Exactly one of the parameters \code{n}, \code{delta}, \code{ref.prob} 
  and \code{conf.prob} must be passed as \code{NULL}, and that parameter is 
  determined from the others. In case of \code{ref.prob} \code{NULL} must be explicitly 
  passed if you want to compute it.
  
  If method \code{"parametric"} a normal distribution is assumed for the investigated
  quantity.
  
  If method \code{"nonparametric"} an arbitrary continuous probability distribution
  is assumed.

  If \code{exact = TRUE} is used, the computations use the exact formulas (5) 
  and (9) of Jennen-Steinmetz and Wellek (2005).
  
  If \code{exact = FALSE} is used, the computations use the approximate formulas (6) 
  and (10) of Jennen-Steinmetz and Wellek (2005).
}
\value{
  Object of class \code{"power.htest"}, a list of the arguments
  (including the computed one) augmented with \code{method} and
  \code{note} elements.
}

\references{
  C. Jennen-Steinmetz, S. Wellek (2005). A new approach to sample size calculation
  for reference interval studies. \emph{Statistics in Medicine} 24:3199-3212.
}

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}

\note{
  \code{uniroot} is used to solve the equations for unknowns, so
  you may see errors from it, notably about inability to bracket the
  root when invalid arguments are given.
}
\seealso{\code{\link{uniroot}}}

\examples{
  ## see Table 1 in Jennen-Steinmetz and Wellek (2005)
  ssize.reference.range(delta = 0.03, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE)
  ## 135 vs 125 (error in Table 1)
  ssize.reference.range(delta = 0.03, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE)
  ssize.reference.range(delta = 0.03, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE)
  ssize.reference.range(delta = 0.03, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE)
  
  ssize.reference.range(delta = 0.025, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE)
  ssize.reference.range(delta = 0.025, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE)
  ssize.reference.range(delta = 0.025, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE)
  ssize.reference.range(delta = 0.025, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE)
  
  ssize.reference.range(delta = 0.02, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE)
  ## 314 vs. 305 (error Table 1?)
  ssize.reference.range(delta = 0.02, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE)
  ssize.reference.range(delta = 0.02, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE)
  ssize.reference.range(delta = 0.02, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE)
  
  ssize.reference.range(delta = 0.015, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE)
  ssize.reference.range(delta = 0.015, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE)
  ssize.reference.range(delta = 0.015, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE)
  ssize.reference.range(delta = 0.015, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE)
  
  ssize.reference.range(delta = 0.01, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE)
  ssize.reference.range(delta = 0.01, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE)
  ssize.reference.range(delta = 0.01, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE)
  ssize.reference.range(delta = 0.01, ref.prob = 0.9, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE)
  
  ssize.reference.range(delta = 0.015, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE)
  ssize.reference.range(delta = 0.015, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE)
  ssize.reference.range(delta = 0.015, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE)
  ssize.reference.range(delta = 0.015, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE)
  
  ssize.reference.range(delta = 0.0125, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE)
  ssize.reference.range(delta = 0.0125, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE)
  ssize.reference.range(delta = 0.0125, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE)
  ssize.reference.range(delta = 0.0125, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE)
  
  ssize.reference.range(delta = 0.01, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE)
  ssize.reference.range(delta = 0.01, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE)
  ssize.reference.range(delta = 0.01, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE)
  ssize.reference.range(delta = 0.01, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE)
  
  ssize.reference.range(delta = 0.0075, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE)
  ssize.reference.range(delta = 0.0075, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE)
  ssize.reference.range(delta = 0.0075, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE)
  ssize.reference.range(delta = 0.0075, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE)
  
  ssize.reference.range(delta = 0.005, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE)
  ssize.reference.range(delta = 0.005, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE)
  ssize.reference.range(delta = 0.005, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE)
  ssize.reference.range(delta = 0.005, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE)
  
  
  ## results are equivalent to one-sided reference range with coverage of 
  ## 95 percent instead of 90 percent; for example
  ssize.reference.range(delta = 0.03, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = TRUE, alternative = "one.sided")
  ## 135 vs 125 (error in Table 1)
  ssize.reference.range(delta = 0.03, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = TRUE, alternative = "one.sided")
  ssize.reference.range(delta = 0.03, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "parametric", exact = FALSE, alternative = "one.sided")
  ssize.reference.range(delta = 0.03, ref.prob = 0.95, conf.prob = 0.9, 
                        method = "nonparametric", exact = FALSE, alternative = "one.sided")

}
\keyword{htest}


