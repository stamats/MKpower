\name{sim.power.xiao.t.test}
\alias{sim.power.xiao.t.test}
\alias{sim.ssize.xiao.t.test}
\title{Power Simulations for Xiao Two-sample t Test}
\usage{
sim.power.xiao.t.test(n1, n2, delta = 1, sd1 = 1, sd2 = 2, mu = 0, sig.level = 0.05, 
                      alternative = c("two.sided", "less", "greater"),
                      iter = 10000, parallel = FALSE, cl = NULL)
sim.ssize.xiao.t.test(delta = 1, sd1 = 1, sd2 = 2, mu = 0, sig.level = 0.05, 
                      power = 0.8, alternative = c("two.sided", "one.sided"),
                      n.delta = 5, step.size = 1, iter = 10000,
                      parallel = FALSE, cl = NULL)
}
\arguments{
 \item{n1}{number of observations of group 1}
 \item{n2}{number of observations of group 2}
 \item{delta}{(expected) true difference in means}
 \item{sd1}{(expected) standard deviation of group 1}
 \item{sd2}{(expected) standard deviation of group 2}
 \item{mu}{ true value of the location shift for the null hypothesis.}
 \item{sig.level}{significance level (Type I error probability)}
 \item{power}{power of test (1 minus Type II error probability)}
 \item{alternative}{one- or two-sided test.}
 \item{iter}{ single integer, number of interations of the simulations.}
 \item{parallel}{logical; use parellel computing via function \code{\link{parRapply}}
  in the case of method \code{"xiao"} as well as in the case of permutation and 
  bootstrap tests.}
 \item{cl}{a cluster object, created by package \code{\link[parallel:makeCluster]{parallel}} 
  or by package \code{\link[snow:makeCluster]{snow}}. If \code{NULL}, a new cluster 
  will be generated with function \code{\link[parallel:makeCluster]{makeCluster}} 
  using the maximum number of available cores - 1.}
 \item{n.delta}{ single integer, neighborhood around the values obtained by 
  \code{\link{ssize.welch.t.test}}.}
 \item{step.size}{ step size for the neighborhood around the values obtained 
  by \code{\link{ssize.welch.t.test}}.}
}

\description{
  Use Monte-Carlo simulations to determine the empirical power and sample 
  size of the Xiao t-test. 
}

\details{
  To find the (optimal) group sample sizes a grid search in a (small) neighborhood 
  around the (optimal) values of the Welch t-test is performed.
}
\value{
  Object of class \code{"power.htest"}, a list of the arguments
  (including the computed one) augmented with \code{method} and
  \code{note} elements.
}

\references{
  Y. Xiao (2018). On the Solution of a Generalized Behrens-Fisher Problem. 
  \emph{Far East Journal of Theoretical Statistics}, \bold{54} (1), 21-140.
}

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}

%\note{}
\seealso{\code{\link{ssize.welch.t.test}}, \code{\link[MKinfer]{xiao.t.test}}}

\examples{
 \dontrun{
   ## delta = 0.5, sd1 = 0.5, sd2 = 1,
   power.welch.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9)
   power.hsu.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9)
   power.welch.t.test2(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9) # 22+87 = 109
   power.hsu.t.test2(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9) # 44+59 = 103
   power.hsu.t.test2(n = c(44, 59), delta = 0.5, sd1 = 0.5, sd2 = 1)
   sim.power.xiao.t.test(n1 = 22, n2 = 87, delta = 0.5, sd1 = 0.5, sd2 = 1, parallel = TRUE)
   ssize.welch.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9) # 35+62 = 97
   ssize.hsu.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9) # 39+62 = 101
   sim.ssize.xiao.t.test(delta = 0.5, sd1 = 0.5, sd2 = 1, power = 0.9, parallel = TRUE)
 }
}
\keyword{htest}


