ssize.welch.t.test <- function (delta = 1, sd1 = 1, sd2 = 1,
                                sig.level = 0.05, power = 0.8,
                                alternative = c("two.sided", "one.sided"),
                                strict = FALSE) {
  if (!is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
    stop("'sig.level' must be numeric in [0, 1]")
  if (!is.numeric(power) || any(0 > power | power > 1))
    stop("'power' must be numeric in [0, 1]")
  if (!is.numeric(sd1) || any(0 > sd1))
    stop("'sd1' must be a positive numeric")
  if (!is.numeric(sd2) || any(0 > sd2))
    stop("'sd2' must be a positive numeric")
  alternative <- match.arg(alternative)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 2 && !is.null(delta))
    delta <- abs(delta)
  p.body <- if (strict && tside == 2)
    quote({
      nu <- (sd1^2/n1+sd2^2/n2)^2/(sd1^4/(n1^2*(n1-1)) + sd2^4/(n2^2*(n2-1)))
      qu <- qt(sig.level/tside, nu, lower.tail = FALSE)
      sd <- sqrt(sd1^2/n1 + sd2^2/n2)
      pt(qu, nu, ncp = delta/sd, lower.tail = FALSE) +
        pt(-qu, nu, ncp = delta/sd, lower.tail = TRUE)
    })
  else quote({
    nu <- (sd1^2/n1+sd2^2/n2)^2/(sd1^4/(n1^2*(n1-1)) + sd2^4/(n2^2*(n2-1)))
    sd <- sqrt(sd1^2/n1 + sd2^2/n2)
    pt(qt(sig.level/tside, nu, lower.tail = FALSE), nu, ncp = delta/sd,
       lower.tail = FALSE)
  })
  n0 <- power.welch.t.test(delta = delta, sd1 = sd1, sd2 = sd2, 
                           sig.level = sig.level, power = power, 
                           alternative = alternative, strict = strict)$n
  n0 <- ceiling(n0)
  fun <- function(n1, n2){ eval(p.body) }
  N1 <- seq.int(from = 2, to = 2*n0, by = 1)
  N2 <- N1
  res <- outer(N1, N2, fun)
  SUM <- outer(N1, N2, "+")
  SSIZE <- outer(N1, N2, function(n1, n2) paste0(n1, ",", n2))
  ind <- which.min(SUM[res > power])
  ssize <- SSIZE[res > power][ind]
  n <- as.integer(unlist(strsplit(ssize, "\\,")))
  power <- power.welch.t.test2(n = n, delta = delta, sd1 = sd1, sd2 = sd2, 
                               sig.level = sig.level, alternative = alternative,
                               strict = strict)$power

  NOTE <- "minimum total sample size (sum of group sample sizes)"
  METHOD <- "Two-sample Welch t test sample size calculation"
  structure(list(n1 = n[1], n2 = n[2], delta = delta, sd1 = sd1, sd2 = sd2, 
                 sig.level = sig.level, power = power, alternative = alternative, 
                 note = NOTE, method = METHOD), class = "power.htest")
}
