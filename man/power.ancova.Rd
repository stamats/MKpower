\name{power.ancova}
\alias{power.ancova}
\title{Power Calculation for ANCOVA}
\description{
 Compute sample size for ANCOVA.
}
\usage{
power.ancova(n = NULL, mu = NULL, var = 1, nr.covs = 1L, group.ratio = NULL, 
             contr.mat = NULL, sig.level = 0.05, power = NULL, n.max = 1000L,
             rel.tol = .Machine$double.eps^0.25)
}
\arguments{
 \item{n}{vector of sample sizes per groups.}
 \item{mu}{vector of mean values of the groups.}
 \item{var}{error variance.}
 \item{nr.covs}{number of covariates (larger or equal than 1).}
 \item{group.ratio}{vector of group sizes relative to group 1; i.e., first 
 entry should always be one. If \code{NULL}, a balanced design is used.}
 \item{contr.mat}{matrix of contrasts (number of columns must be idential to 
 number of groups). If \code{NULL}, standard ANCOVA contrasts are used; see
 examples below.}
 \item{sig.level}{significance level (type I error probability)}
 \item{power}{power of test (1 minus type II error probability)}
 \item{n.max}{maximum sample size considered in the computations.}
 \item{rel.tol}{relative tolerance passed to function \code{integrate}.}
}

\details{
  Exactly one of the parameters \code{n} and \code{power} must be passed as
  \code{NULL}, and that parameter is determined from the other.
  
  The function includes an implementation of the exact approach of Shieh (2020).
  It is based on the code provided in the supplement of Shieh (2020), but
  uses \code{integrate} instead of the trapezoid rule and \code{uniroot} for
  finding the required sample size.
}
\value{
  Object of class \code{"power.htest"}, a list of the arguments
  (including the computed one) augmented with a \code{note} element.
}

\references{
  G. Shieh (2020). Power Analysis and Sample Size Planning in ANCOVA Designs. 
  \emph{Psychometrika} \bold{85}:101-120. \doi{10.1007/s11336-019-09692-3}.
  
  S.E. Maxwell and H.D. Delaney (2004). \emph{Designing experiments and analyzing 
  data: A model comparison perspective} (2nd ed.). Mahwah, 
  NJ: Lawrence Erlbaum Associates.
}

\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}

%\note{}
\seealso{\code{\link[stats]{power.anova.test}}, \code{\link[stats]{power.t.test}}}

\examples{
## Default matrix of contrasts
## 3 groups
cbind(rep(1,2), -diag(2))
## 4 groups
cbind(rep(1,3), -diag(3))

## Table 1 in Shieh (2020)
power.ancova(mu=c(400, 450, 500), var = 9900, power = 0.8)
power.ancova(n = rep(63/3, 3), mu=c(400, 450, 500), var = 9900)
power.ancova(mu=c(400, 450, 500), var = 9900, power = 0.8, nr.covs = 10)
power.ancova(n = rep(72/3, 3), mu=c(400, 450, 500), var = 9900, nr.covs = 10)

## Table 2 in Shieh (2020)
power.ancova(mu=c(400, 450, 500), var = 7500, power = 0.8)
power.ancova(n = rep(48/3, 3), mu=c(400, 450, 500), var = 7500)
power.ancova(mu=c(400, 450, 500), var = 7500, power = 0.8, nr.covs = 10)
power.ancova(n = rep(60/3, 3), mu=c(400, 450, 500), var = 7500, nr.covs = 10)

## Table 3 in Shieh (2020)
power.ancova(mu=c(400, 450, 500), var = 1900, power = 0.8)
power.ancova(n = rep(18/3, 3), mu=c(400, 450, 500), var = 1900)
power.ancova(mu=c(400, 450, 500), var = 1900, power = 0.8, nr.covs = 10)
power.ancova(n = rep(27/3, 3), mu=c(400, 450, 500), var = 1900, nr.covs = 10)

## ANOVA approach for Table 1-3
power.anova.test(groups = 3, between.var = var(c(400, 450, 500)), 
                 within.var = 10000, power = 0.8)
power.anova.test(n = 63/3, groups = 3, between.var = var(c(400, 450, 500)), 
                 within.var = 10000)

## Table 4 in Shieh (2020)
power.ancova(mu=c(410, 450, 490), var = 9900, power = 0.8)
power.ancova(n = rep(96/3, 3), mu=c(410, 450, 490), var = 9900)
power.ancova(mu=c(410, 450, 490), var = 9900, power = 0.8, nr.covs = 10)
power.ancova(n = rep(105/3, 3), mu=c(410, 450, 490), var = 9900, nr.covs = 10)

## Table 5 in Shieh (2020)
power.ancova(mu=c(410, 450, 490), var = 7500, power = 0.8)
power.ancova(n = rep(72/3, 3), mu=c(410, 450, 490), var = 7500)
power.ancova(mu=c(410, 450, 490), var = 7500, power = 0.8, nr.covs = 10)
power.ancova(n = rep(84/3, 3), mu=c(410, 450, 490), var = 7500, nr.covs = 10)

## Table 6 in Shieh (2020)
power.ancova(mu=c(410, 450, 490), var = 1900, power = 0.8)
power.ancova(n = rep(24/3, 3), mu=c(410, 450, 490), var = 1900)
power.ancova(mu=c(410, 450, 490), var = 1900, power = 0.8, nr.covs = 10)
power.ancova(n = rep(33/3, 3), mu=c(410, 450, 490), var = 1900, nr.covs = 10)

## ANOVA approach for Table 4-6
power.anova.test(groups = 3, between.var = var(c(410, 450, 490)), 
                 within.var = 10000, power = 0.8)
power.anova.test(n = 96/3, groups = 3, between.var = var(c(410, 450, 490)), 
                 within.var = 10000)

###############################################################################
## Example from Maxwell and Delaney (2004) according to Shieh (2020)
###############################################################################
## ANCOVA (balanced design)
power.ancova(n = rep(30/3, 3), mu=c(7.5366, 11.9849, 13.9785), var = 29.0898)
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.8)
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.9)

## ANOVA
power.anova.test(n = 30/3, groups = 3, between.var = var(c(7.5366, 11.9849, 13.9785)), 
                 within.var = 29.0898)
power.anova.test(groups = 3, between.var = var(c(7.5366, 11.9849, 13.9785)), 
                 within.var = 29.0898, power = 0.8)
power.anova.test(groups = 3, between.var = var(c(7.5366, 11.9849, 13.9785)), 
                 within.var = 29.0898, power = 0.9)
                 
## ANCOVA - imbalanced design
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.8, 
             group.ratio = c(1, 1.25, 1.5))
power.ancova(n = c(13, 16, 19), mu=c(7.5366, 11.9849, 13.9785), var = 29.0898,  
             group.ratio = c(1, 1.25, 1.5))
power.ancova(mu=c(7.5366, 11.9849, 13.9785), var = 29.0898, power = 0.8, 
             group.ratio = c(1, 0.8, 2/3))
power.ancova(n = c(17, 14, 12), mu=c(7.5366, 11.9849, 13.9785), var = 29.0898,  
             group.ratio = c(1, 0.8, 2/3))
}
\keyword{htest}


