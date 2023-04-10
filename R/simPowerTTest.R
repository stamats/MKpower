sim.power.t.test <- function(nx, rx, rx.H0 = NULL, ny, ry, ry.H0 = NULL, 
                             sig.level = 0.05, conf.int = FALSE, mu = 0, 
                             alternative = c("two.sided", "less", "greater"), 
                             iter = 10000){
  alternative <- match.arg(alternative)
  if(length(sig.level) != 1)
    stop("'sig.level' has to be of length 1 (type I error)")
  if(sig.level <= 0 | sig.level >= 1)
    stop("'sig.level' has to be in (0, 1)")
  stopifnot(is.function(rx), is.function(ry))
  stopifnot(is.numeric(nx), length(nx) == 1, nx >= 1)
  nx <- trunc(nx)
  stopifnot(is.numeric(ny), length(ny) == 1, ny >= 1)
  ny <- trunc(ny)
  stopifnot(is.numeric(mu), length(mu) == 1)
  stopifnot(is.numeric(iter), length(iter) == 1, iter >= 1)
  iter <- trunc(iter)
  
  if(conf.int){
    conf.level <- 1 - sig.level
    alpha <- sig.level
  }else{
    conf.level <- NA
  }
  data.x <- matrix(rx(nx*iter), nrow = iter)
  data.y <- matrix(ry(ny*iter), nrow = iter)
  if(!is.null(rx.H0) & !is.null(ry.H0)){
    data.x.H0 <- matrix(rx.H0(nx*iter), nrow = iter)
    data.y.H0 <- matrix(ry.H0(ny*iter), nrow = iter)
  }
  ## classical 2-sample t-test
  res <- row_t_equalvar(data.x, data.y, 
                        alternative = alternative, null = mu, 
                        conf.level = conf.level)
  if(!is.null(rx.H0) & !is.null(ry.H0)){
    res.H0 <- row_t_equalvar(data.x.H0, data.y.H0, 
                             alternative = alternative, 
                             null = mu, conf.level = conf.level)
    CLASSICAL <- list("H1" = res, "H0" = res.H0)
  }else{
    CLASSICAL <- list("H1" = res)  
  }
  ## Welch t-test
  res <- row_t_welch(data.x, data.y, 
                     alternative = alternative, null = mu, 
                     conf.level = conf.level)
  if(!is.null(rx.H0) & !is.null(ry.H0)){
    res.H0 <- row_t_welch(data.x.H0, data.y.H0, 
                          alternative = alternative, 
                          null = mu, conf.level = conf.level)
    WELCH <- list("H1" = res, "H0" = res.H0)
  }else{
    WELCH <- list("H1" = res)
  }
  ## Hsu t-test
  df <- min(nx, ny) - 1
  res$df <- df
  if(alternative == "less"){
    res$pvalue <- pt(res$statistic, df = df)
    if(conf.int){
      res$conf.high <- mu + (res$statistic + qt(conf.level, df = df))*res$stderr
    }
  }else if(alternative == "greater") {
    res$pvalue <- pt(res$statistic, df = df, lower.tail = FALSE)
    if(conf.int){
      res$conf.low <- mu + (res$statistic - qt(conf.level, df = df))*res$stderr
    }
  }
  else{
    res$pvalue <- 2 * pt(-abs(res$statistic), df = df)
    if(conf.int){
      crit <- qt(1 - alpha/2, df)
      cint <- mu + (res$statistic + c(-crit, crit))*res$stderr
      res$conf.low <- cint[1]
      res$conf.high <- cint[2]
    }
  }
  if(!is.null(rx.H0) & !is.null(ry.H0)){
    if(alternative == "less"){
      res.H0$pvalue <- pt(res.H0$statistic, df = df)
      if(conf.int){
        res.H0$conf.high <- mu + (res.H0$statistic + qt(conf.level, df = df))*res.H0$stderr
      }
    }else if(alternative == "greater") {
      res.H0$pvalue <- pt(res.H0$statistic, df = df, lower.tail = FALSE)
      if(conf.int){
        res.H0$conf.low <- mu + (res.H0$statistic - qt(conf.level, df = df))*res.H0$stderr
      }
    }
    else{
      res.H0$pvalue <- 2 * pt(-abs(res.H0$statistic), df = df)
      if(conf.int){
        crit <- qt(1 - alpha/2, df)
        cint <- mu + (res.H0$statistic + c(-crit, crit))*res.H0$stderr
        res.H0$conf.low <- cint[1]
        res.H0$conf.high <- cint[2]
      }
    }
    HSU <- list("H1" = res, "H0" = res.H0)
  }else{
    HSU <- list("H1" = res)
  }
  SetUp <- c("nx" = nx, "rx" = rx, "rx.H0" = rx.H0, 
             "ny" = ny, "ry" = ry, "ry.H0" = ry.H0, 
             "sig.level" = sig.level, "mu" = mu,
             "alternative" = alternative, "iter" = iter)
  res <- list(Classical = CLASSICAL, Welch = WELCH, Hsu = HSU, SetUp = SetUp)
  class(res) <- "sim.power.ttest"
  res
}

print.sim.power.ttest <- function(x, digits = getOption("digits"), ...){
  cat("\n    Simulation Set-up\n")
  alpha <- x$SetUp$sig.level
  y0 <- x$SetUp
  iter <- x$SetUp$iter
  cat(paste(format(names(y0), width = 15L, justify = "right"), 
            format(y0, digits = digits), sep = " = "), sep = "\n")
  cat("\n    Classical Two-sample t-Test\n")
  y1 <- c("emp.power" = sum(x$Classical$H1$pvalue < alpha)/iter)
  if(!is.null(x$Classical$H0)){
    y1 <- c(y1, "emp.type.I.error" = sum(x$Classical$H0$pvalue < alpha)/iter)
  }
  cat(paste(format(names(y1), width = 15L, justify = "right"), 
            format(y1, digits = digits), sep = " = "), sep = "\n")
  cat("\n    Welch Two-sample t-Test\n")
  y2 <- c("emp.power" = sum(x$Welch$H1$pvalue < alpha)/iter)
  if(!is.null(x$Welch$H0)){
    y2 <- c(y2, "emp.type.I.error" = sum(x$Welch$H0$pvalue < alpha)/iter)
  }
  cat(paste(format(names(y2), width = 15L, justify = "right"), 
            format(y2, digits = digits), sep = " = "), sep = "\n")
  cat("\n    Hsu Two-sample t-Test\n")
  y3 <- c("emp.power" = sum(x$Hsu$H1$pvalue < alpha)/iter)
  if(!is.null(x$Hsu$H0)){
    y3 <- c(y3, "emp.type.I.error" = sum(x$Hsu$H0$pvalue < alpha)/iter)
  }
  cat(paste(format(names(y3), width = 15L, justify = "right"), 
            format(y3, digits = digits), sep = " = "), sep = "\n")
  cat("\n")
  invisible(c(y1, y2, y3))
}
