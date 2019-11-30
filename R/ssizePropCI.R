ssize.propCI <- function(prop, width, conf.level = 0.95,  method = "wald-cc"){
  METHODS <- c("wald", "wald-cc")
  method <- pmatch(method, METHODS)
  
  if (is.na(method))
    stop("invalid method")
  if (method == -1)
    stop("ambiguous method")
  
  alpha <- 1-conf.level
  if(method == 1){ # wald
    n <- 4*qnorm(1-alpha/2)^2*prop*(1-prop)/width^2
  }
  if(method == 2){ # wald-cc
    z.alpha <- qnorm(1-alpha/2)
    TERM <- 2*sqrt(width*prop*(1-prop)*z.alpha^2 + prop^2*(1-prop)^2*z.alpha^4)/width^2
    n <- 2*prop*(1-prop)*z.alpha^2/width^2 + 1/width + TERM
  }
  METHOD <- paste("Sample size calculation by method of", METHODS[method])
  NOTE <- "Two-sided confidence interval"
  structure(list(n = n, prop = prop, width = width, conf.level = conf.level,
                 note = NOTE, method = METHOD), class = "power.htest")
}
