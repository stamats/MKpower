power.welch.t.test2 <- function (n = NULL, delta = NULL, sd1 = 1, sd2 = 1,
                                sig.level = 0.05, power = NULL,
                                alternative = c("two.sided", "one.sided"),
                                strict = FALSE, tol = .Machine$double.eps^0.25, 
                                abstol = .Machine$double.eps^0.5) {
  if (sum(sapply(list(n, delta, sd1, sd2, power, sig.level), is.null)) != 1)
    stop("exactly one of 'n', 'delta', 'sd1', 'sd2', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
    stop("'sig.level' must be numeric in [0, 1]")
  if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 1))
    stop("'power' must be numeric in [0, 1]")
  if (!is.null(sd1) && !is.numeric(sd1) || any(0 > sd1))
    stop("'sd1' must be a positive numeric")
  if (!is.null(sd2) && !is.numeric(sd2) || any(0 > sd2))
    stop("'sd2' must be a positive numeric")
  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 2 && !is.null(delta))
    delta <- abs(delta)
  p.body <- if (strict && tside == 2)
    quote({
      n1 <- n[1]
      n2 <- n[2]
      nu <- (sd1^2/n1+sd2^2/n2)^2/(sd1^4/(n1^2*(n1-1)) + sd2^4/(n2^2*(n2-1)))
      qu <- qt(sig.level/tside, nu, lower.tail = FALSE)
      sd <- sqrt(sd1^2/n1 + sd2^2/n2)
      pt(qu, nu, ncp = delta/sd, lower.tail = FALSE) +
        pt(-qu, nu, ncp = delta/sd, lower.tail = TRUE)
    })
  else quote({
    n1 <- n[1]
    n2 <- n[2]
    nu <- (sd1^2/n1+sd2^2/n2)^2/(sd1^4/(n1^2*(n1-1)) + sd2^4/(n2^2*(n2-1)))
    sd <- sqrt(sd1^2/n1 + sd2^2/n2)
    pt(qt(sig.level/tside, nu, lower.tail = FALSE), nu, ncp = delta/sd,
       lower.tail = FALSE)
  })
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- optim(c(5,5), function(n){ abs(eval(p.body) - power) }, 
               method = "Nelder-Mead", control = list(abstol = abstol))$par
  else if (is.null(sd1))
    sd1 <- uniroot(function(sd1) eval(p.body) - power, delta *
                    c(1e-07, 1e+07), tol = tol, extendInt = "downX")$root
  else if (is.null(sd2))
    sd2 <- uniroot(function(sd2) eval(p.body) - power, delta *
                    c(1e-07, 1e+07), tol = tol, extendInt = "downX")$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power,
                     sqrt(sd1^2 + sd2^2) * c(1e-07, 1e+07), tol = tol, 
                     extendInt = "upX")$root
  else if (is.null(sig.level))
    sig.level <- uniroot(function(sig.level) eval(p.body) -
                           power, c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root
  else stop("internal error", domain = NA)
  NOTE <- "unequal sample sizes per group (allowed)"
  METHOD <- "Two-sample Welch t test power calculation"
  structure(list(n1 = n[1], n2 = n[2], delta = delta, sd1 = sd1, sd2 = sd2, 
                 sig.level = sig.level, power = power, alternative = alternative, 
                 note = NOTE, method = METHOD), class = "power.htest")
}
