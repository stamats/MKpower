\name{MKpower-package}
\alias{MKpower-package}
\alias{MKpower}
\docType{package}
\title{
Power Analysis and Sample Size Calculation.
}
\description{
Power analysis and sample size calculation for Welch and Hsu 
(Hedderich and Sachs (2018), ISBN:978-3-662-56657-2) t-tests including 
Monte-Carlo simulations of empirical power and type-I-error. 
Power and sample size calculation for Wilcoxon rank sum and signed rank tests 
via Monte-Carlo simulations. Power and sample size required for the evaluation 
of a diagnostic test(-system) (Flahault et al. (2005), 
<doi:10.1016/j.jclinepi.2004.12.009>; Dobbin and Simon (2007), 
<doi:10.1093/biostatistics/kxj036>) as well as for a single proportion 
(Fleiss et al. (2003), ISBN:978-0-471-52629-2) and comparing 
two negative binomial rates (Zhu and Lakkis (2014), <doi:10.1002/sim.5947>).
}
\details{
\tabular{ll}{
Package: \tab MKpower\cr
Type: \tab Package\cr
Version: \tab 0.4\cr
Date: \tab 2020-01-23\cr
Depends: \tab R(>= 3.5.0)\cr
Imports: \tab stats, matrixTests, ggplot2, MKdescr, MKinfer(>= 0.4), qqplotr, coin\cr
Suggests: \tab knitr, rmarkdown\cr
License: \tab LGPL-3\cr
URL: \tab http://www.stamats.de/\cr
}
library(MKpower)
}
\author{ Matthias Kohl \url{http://www.stamats.de}

Maintainer: Matthias Kohl  \email{matthias.kohl@stamats.de}}
%\references{
%~~ Literature or other references for background information ~~
%}
%~~ Optionally other keywords from doc/KEYWORDS, one per line
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
\keyword{package}
