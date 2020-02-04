# MKpower
The repository includes the development version of R package MKpower

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MKpower)](http://cran.r-project.org/package=MKpower)
[![cran checks](https://cranchecks.info/badges/summary/MKpower)](https://cran.r-project.org/web/checks/check_results_MKpower.html)

## Description
Power analysis and sample size calculation for Welch and Hsu 
(Hedderich and Sachs (2018), ISBN:978-3-662-56657-2) t-tests including 
Monte-Carlo simulations of empirical power and type-I-error. Power and sample 
size calculation for Wilcoxon rank sum and signed rank tests via Monte-Carlo 
simulations. Power and sample size required for the evaluation of a diagnostic 
test(-system) (Flahault et al. (2005), <doi:10.1016/j.jclinepi.2004.12.009>; 
Dobbin and Simon (2007), <doi:10.1093/biostatistics/kxj036>) as well as for a 
single proportion (Fleiss et al. (2003), ISBN:978-0-471-52629-2) and 
comparing two negative binomial rates (Zhu and Lakkis (2014), <doi:10.1002/sim.5947>).

## Installation

```{r, eval = FALSE}
## Installation of CRAN version
install.packages("MKpower")

## Development version from GitHub
# install.packages("remotes")
remotes::install_github("stamats/MKpower")
```

## Getting started

```{r}
library(MKpower)
```
