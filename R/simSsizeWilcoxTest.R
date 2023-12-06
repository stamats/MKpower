sim.ssize.wilcox.test <- function(rx, ry = NULL, mu = 0, sig.level = 0.05, power = 0.8, 
                                  type = c("two.sample", "one.sample", "paired"), 
                                  alternative = c("two.sided", "less", "greater"),
                                  n.min = 10, n.max = 200, step.size = 10, 
                                  iter = 10000, BREAK = TRUE, exact = NA, 
                                  correct = TRUE){
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  if(length(sig.level) != 1)
    stop("'sig.level' has to be of length 1 (type I error)")
  if(sig.level <= 0 | sig.level >= 1)
    stop("'sig.level' has to be in (0, 1)")
  if(length(power) != 1)
    stop("'power' has to be of length 1")
  if(power <= 0 | power >= 1)
    stop("'power' has to be in (0, 1)")
  stopifnot(is.function(rx))
  stopifnot(is.numeric(mu), length(mu) == 1)
  stopifnot(is.numeric(n.min), length(n.min) == 1, n.min >= 1)
  stopifnot(is.numeric(n.max), length(n.max) == 1, n.max > n.min)
  stopifnot(is.numeric(step.size), length(step.size) == 1, step.size < n.max-n.min)
  stopifnot(is.numeric(iter), length(iter) == 1, iter >= 1)
  n.min <- trunc(n.min)
  n.max <- trunc(n.max)
  step.size <- trunc(step.size)
  iter <- trunc(iter)
  
  if(type == "two.sample"){
    stopifnot(!is.null(ry))
    stopifnot(is.function(ry))
    ns <- seq.int(from = n.min, to = n.max, by = step.size)
    empPower <- numeric(length(ns))
    for(i in seq_len(length(ns))){
      data.x <- matrix(rx(ns[i]*iter), nrow = iter)
      data.y <- matrix(ry(ns[i]*iter), nrow = iter)
      res <- row_wilcoxon_twosample(data.x, data.y, 
                                    alternative = alternative, null = mu,
                                    exact = exact, correct = correct)[,"pvalue"]
      empPower[i] <- sum(res < sig.level)/iter
      if(empPower[i] > power) if(BREAK) break
    }
    METHOD <- "Wilcoxon rank sum test"
  }
  if(type == "one.sample"){
    ns <- seq.int(from = n.min, to = n.max, by = step.size)
    empPower <- numeric(length(ns))
    for(i in seq_len(length(ns))){
      data.x <- matrix(rx(ns[i]*iter), nrow = iter)
      res <- row_wilcoxon_onesample(data.x, alternative = alternative, null = mu,
                                    exact = exact, correct = correct)[,"pvalue"]
      empPower[i] <- sum(res < sig.level)/iter
      if(empPower[i] > power) if(BREAK) break
    }
    METHOD <- "Wilcoxon signed rank test"
  }
  if(type == "paired"){
    ns <- seq.int(from = n.min, to = n.max, by = step.size)
    empPower <- numeric(length(ns))
    for(i in seq_len(length(ns))){
      data.xy <- matrix(rx(ns[i]*iter), nrow = iter)
      res <- row_wilcoxon_onesample(data.xy, 
                                    alternative = alternative, null = mu,
                                    exact = exact, correct = correct)[,"pvalue"]
      empPower[i] <- sum(res < sig.level)/iter
      if(empPower[i] > power) if(BREAK) break
    }
    METHOD <- "Wilcoxon signed rank test"
  }
  empPower[empPower == 0] <- NA
  names(empPower) <- ns
  empPower
  NOTE <- "n is number in *each* group"
  if(is.null(ry))
    res <- structure(list(n = ns[!is.na(empPower)], rx = body(rx), 
                          sig.level = sig.level, emp.power = empPower[!is.na(empPower)],
                          alternative = alternative, note = NOTE, 
                          method = METHOD), class = "power.htest")
  else
    res <- structure(list(n = ns[!is.na(empPower)], rx = body(rx), ry = body(ry), 
                          sig.level = sig.level, emp.power = empPower[!is.na(empPower)],
                          alternative = alternative, note = NOTE, 
                          method = METHOD), class = "power.htest")
  res
}

