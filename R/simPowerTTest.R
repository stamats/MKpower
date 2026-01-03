sim.power.t.test <- function(nx, rx, rx.H0 = NULL, ny, ry, ry.H0 = NULL, 
                             sig.level = 0.05, conf.int = FALSE, mu = 0, 
                             alternative = c("two.sided", "less", "greater"), 
                             methods = c("student", "welch", "hsu", "xiao"), 
                             iter = 10000, R = 9999, useCombn = FALSE, 
                             parallel = TRUE, cl = NULL){
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
  
  ## generate the data
  data.x <- matrix(rx(nx*iter), nrow = iter)
  data.y <- matrix(ry(ny*iter), nrow = iter)
  if(!is.null(rx.H0) & !is.null(ry.H0)){
    data.x.H0 <- matrix(rx.H0(nx*iter), nrow = iter)
    data.y.H0 <- matrix(ry.H0(ny*iter), nrow = iter)
  }
  
  ## set-up for parallel computing
  PARA <- any(c("xiao", "perm.student", "perm.welch", "boot.student", 
                "boot.welch") %in% methods)
  if(parallel && PARA){
    if(!("package:MKinfer" %in% search())){
      DETACH <- TRUE
      attachNamespace("MKinfer")
    }else{
      DETACH <- FALSE
    }
    if(is.null(cl)){
      STOP <- TRUE
      ncores <- detectCores()-1
      cl <- makePSOCKcluster(rep("localhost", ncores))
    }else{
      STOP <- FALSE
    }
  }
  
  RESULT <- list()
  
  #########################################################
  ## Student 2-sample t-test
  #########################################################
  if("student" %in% methods){
    res <- row_t_equalvar(data.x, data.y, 
                          alternative = alternative, null = mu, 
                          conf.level = conf.level)
    if(!is.null(rx.H0) & !is.null(ry.H0)){
      res.H0 <- row_t_equalvar(data.x.H0, data.y.H0, 
                               alternative = alternative, 
                               null = mu, conf.level = conf.level)
      STUDENT <- list("H1" = res, "H0" = res.H0)
    }else{
      STUDENT <- list("H1" = res)  
    }
    RESULT <- c(RESULT, list(Student = STUDENT))
  }
  
  #########################################################
  ## Welch 2-sample t-test
  #########################################################
  if("welch" %in% methods){
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
    RESULT <- c(RESULT, list(Welch = WELCH))
  }
  
  #########################################################
  ## Hsu 2-sample t-test
  #########################################################
  if("hsu" %in% methods){
    df <- min(nx, ny) - 1
    if("welch" %in% methods){
      res <- WELCH$H1 
    }else{
      res <- row_t_welch(data.x, data.y, 
                         alternative = alternative, null = mu, 
                         conf.level = conf.level)
    }
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
        res$conf.low <- mu + (res$statistic - crit)*res$stderr
        res$conf.high <- mu + (res$statistic + crit)*res$stderr
      }
    }
    if(!is.null(rx.H0) & !is.null(ry.H0)){
      if("welch" %in% methods){
        res.H0 <- WELCH$H0
      }else{
        res.H0 <- row_t_welch(data.x.H0, data.y.H0, 
                              alternative = alternative, 
                              null = mu, conf.level = conf.level)
      }
      res.H0$df <- df
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
          res.H0$conf.low <- mu + (res.H0$statistic - crit)*res.H0$stderr
          res.H0$conf.high <- mu + (res.H0$statistic + crit)*res.H0$stderr
        }
      }
      HSU <- list("H1" = res, "H0" = res.H0)
    }else{
      HSU <- list("H1" = res)
    }
    RESULT <- c(RESULT, list(Hsu = HSU))
  }
  
  #########################################################
  ## Xiao two-sample t-test
  #########################################################
  if("xiao" %in% methods){
    if("welch" %in% methods || "hsu" %in% methods){
      if("welch" %in% methods){
        res <- WELCH$H1
      }else{
        res <- HSU$H1
      }
    }else{
      res <- row_t_welch(data.x, data.y, 
                         alternative = alternative, null = mu, 
                         conf.level = conf.level)
    }
    
    vx <- res$var.x
    vy <- res$var.y
    df <- ifelse(nx*(nx-1)/vx <= ny*(ny-1)/vy, 
                 (ny-1)*(1 + ny/nx*vx/vy), (nx-1)*(1 + nx/ny*vy/vx))
    res$df <- df
    if(alternative == "less"){
      res$pvalue <- pgt(res$statistic, n1 = nx, n2 = ny, sd1 = sqrt(vx), 
                        sd2 = sqrt(vy), parallel = parallel, cl = cl)
      if(conf.int){
        res$conf.high <- mu + (res$statistic + qgt(conf.level, n1 = nx, n2 = ny, 
                                                   sd1 = sqrt(vx), sd2 = sqrt(vy), 
                                                   parallel = parallel, cl = cl))*res$stderr
      }
    }else if(alternative == "greater") {
      res$pvalue <- pgt(res$statistic, n1 = nx, n2 = ny, sd1 = sqrt(vx), 
                        sd2 = sqrt(vy), lower.tail = FALSE, parallel = parallel, 
                        cl = cl)
      if(conf.int){
        res$conf.low <- mu + (res$statistic - qgt(conf.level, n1 = nx, n2 = ny, 
                                                  sd1 = sqrt(vx), sd2 = sqrt(vy), 
                                                  parallel = parallel, cl = cl))*res$stderr
      }
    }
    else{
      res$pvalue <- 2 * pgt(-abs(res$statistic), n1 = nx, n2 = ny, 
                            sd1 = sqrt(vx), sd2 = sqrt(vy), 
                            parallel = parallel, cl = cl)
      if(conf.int){
        crit <- qgt(1 - alpha/2, n1 = nx, n2 = ny, sd1 = sqrt(vx), sd2 = sqrt(vy), 
                    parallel = parallel, cl = cl)
        res$conf.low <- mu + (res$statistic - crit)*res$stderr
        res$conf.high <- mu + (res$statistic + crit)*res$stderr
      }
    }
    if(!is.null(rx.H0) & !is.null(ry.H0)){
      if("welch" %in% methods || "hsu" %in% methods){
        if("welch" %in% methods){
          res.H0 <- WELCH$H0
        }else{
          res.H0 <- HSU$H0
        }
      }else{
        res.H0 <- row_t_welch(data.x.H0, data.y.H0, 
                              alternative = alternative, 
                              null = mu, conf.level = conf.level)
      }
      vx <- res.H0$var.x
      vy <- res.H0$var.y
      df <- ifelse(nx*(nx-1)/vx <= ny*(ny-1)/vy, 
                   (ny-1)*(1 + ny/nx*vx/vy), (nx-1)*(1 + nx/ny*vy/vx))
      res.H0$df <- df
      if(alternative == "less"){
        res.H0$pvalue <- pgt(res.H0$statistic, n1 = nx, n2 = ny, 
                             sd1 = sqrt(vx), sd2 = sqrt(vy), 
                             parallel = parallel, cl = cl)
        if(conf.int){
          res.H0$conf.high <- mu + (res.H0$statistic + 
                                      qgt(conf.level, n1 = nx, n2 = ny, 
                                          sd1 = sqrt(vx), sd2 = sqrt(vy), 
                                          parallel = parallel, cl = cl))*res.H0$stderr
        }
      }else if(alternative == "greater") {
        res.H0$pvalue <- pgt(res.H0$statistic, n1 = nx, n2 = ny, 
                             sd1 = sqrt(vx), sd2 = sqrt(vy), lower.tail = FALSE, 
                             parallel = parallel, cl = cl)
        if(conf.int){
          res.H0$conf.low <- mu + (res.H0$statistic - 
                                     qgt(conf.level, n1 = nx, n2 = ny, 
                                         sd1 = sqrt(vx), sd2 = sqrt(vy), 
                                         parallel = parallel, cl = cl))*res.H0$stderr
        }
      }
      else{
        res.H0$pvalue <- 2 * pgt(-abs(res.H0$statistic), n1 = nx, n2 = ny, 
                                 sd1 = sqrt(vx), sd2 = sqrt(vy), 
                                 parallel = parallel, cl = cl)
        if(conf.int){
          crit <- qgt(1 - alpha/2, n1 = nx, n2 = ny, sd1 = sqrt(vx), sd2 = sqrt(vy), 
                      parallel = parallel, cl = cl)
          res.H0$conf.low <- mu + (res.H0$statistic - crit)*res.H0$stderr
          res.H0$conf.high <- mu + (res.H0$statistic + crit)*res.H0$stderr
        }
      }
      XIAO <- list("H1" = res, "H0" = res.H0)
    }else{
      XIAO <- list("H1" = res)
    }
    RESULT <- c(RESULT, list(Xiao = XIAO))
  }
  
  #########################################################
  ## Permutation Student 2-sample t-test
  #########################################################
  if("perm.student" %in% methods){
    if(parallel){
      row.perm.stud <- function(xy, nx, ny, conf.level, mu, alternative, R, useCombn){
        tmp <- perm.t.test(x = xy[1:nx], 
                           y = xy[(nx+1):(nx+ny)], var.equal = TRUE, 
                           conf.level = conf.level, mu = mu, 
                           alternative = alternative, R = R, 
                           useCombn = useCombn)
        res <- c(tmp$statistic, tmp$perm.p.value, 
                 tmp$R.true, tmp$perm.conf.int[1], tmp$perm.conf.int[2],
                 tmp$perm.estimate, tmp$perm.stderr)
        names(res) <- c("statistic", "pvalue", "R", "conf.low",
                        "conf.high", "mean.diff", "stderr")
        res
      }
      data.comb <- cbind(data.x, data.y)
      res <- parRapply(cl = cl, x = data.comb, FUN = row.perm.stud, 
                       nx = nx, ny = ny, conf.level = 1-sig.level, mu = mu, 
                       R = R, useCombn = useCombn)
    }else{
      res <- matrix(NA, nrow = iter, ncol = 7)
      for(i in 1:iter){
        tmp <- perm.t.test(x = data.x[i,], y = data.y[i,], var.equal = TRUE, 
                           conf.level = 1-sig.level, mu = mu, 
                           alternative = alternative, R = R, useCombn = useCombn)
        res[,i] <- c(tmp$statistic, tmp$perm.p.value, 
                     tmp$R.true, tmp$perm.conf.int[1], tmp$perm.conf.int[2],
                     tmp$perm.estimate, tmp$perm.stderr)
      }
      colnames(res) <- c("statistic", "pvalue", "R", "conf.low",
                        "conf.high", "mean.diff", "stderr")
    }
    if(!is.null(rx.H0) & !is.null(ry.H0)){
      if(parallel){
        data.comb.H0 <- cbind(data.x.H0, data.y.H0)
        res.H0 <- parRapply(cl = cl, x = data.comb.H0, FUN = row.perm.stud, 
                            nx = nx, ny = ny, conf.level = 1-sig.level, mu = mu, 
                            R = R, useCombn = useCombn)
      }else{
        res.H0 <- matrix(NA, nrow = iter, ncol = 7)
        for(i in 1:iter){
          tmp <- perm.t.test(x = data.x.H0[i,], y = data.y.H0[i,],
                             var.equal = TRUE, conf.level = 1-sig.level, 
                             conf.int = conf.int, mu = mu, 
                             alternative = alternative, R = R, 
                             useCombn = useCombn)
          res.H0[,i] <- c(tmp$statistic, tmp$perm.p.value, 
                          tmp$R.true, tmp$perm.conf.int[1], tmp$perm.conf.int[2],
                          tmp$perm.estimate, tmp$perm.stderr)
        }
        colnames(res.H0) <- c("statistic", "pvalue", "R", "conf.low",
                             "conf.high", "mean.diff", "stderr")
      }
      PERM.STUDENT <- list("H1" = res, "H0" = res.H0)
    }else{
      PERM.STUDENT <- list("H1" = res)  
    }
    RESULT <- c(RESULT, list(Perm.Student = PERM.STUDENT))
  }
  
  #########################################################
  ## Permutation Welch 2-sample t-test
  #########################################################
  if("perm.student" %in% methods){
    if(parallel){
      row.perm.welch <- function(xy, nx, ny, conf.level, mu, alternative, R, useCombn){
        tmp <- perm.t.test(x = xy[1:nx], 
                           y = xy[(nx+1):(nx+ny)], var.equal = FALSE, 
                           conf.level = conf.level, mu = mu, 
                           alternative = alternative, R = R, 
                           useCombn = useCombn)
        res <- c(tmp$statistic, tmp$perm.p.value, 
                 tmp$R.true, tmp$perm.conf.int[1], tmp$perm.conf.int[2],
                 tmp$perm.estimate, tmp$perm.stderr)
        names(res) <- c("statistic", "pvalue", "R", "conf.low",
                        "conf.high", "mean.diff", "stderr")
        res
      }
      data.comb <- cbind(data.x, data.y)
      res <- parRapply(cl = cl, x = data.comb, FUN = row.perm.welch, 
                       nx = nx, ny = ny, conf.level = 1-sig.level, mu = mu, 
                       R = R, useCombn = useCombn)
    }else{
      res <- matrix(NA, nrow = iter, ncol = 7)
      for(i in 1:iter){
        tmp <- perm.t.test(x = data.x[i,], y = data.y[i,], var.equal = FALSE, 
                           conf.level = 1-sig.level, conf.int = conf.int, mu = mu, 
                           alternative = alternative, R = R, useCombn = useCombn)
        res[,i] <- c(tmp$statistic, tmp$perm.p.value, 
                     tmp$R.true, tmp$perm.conf.int[1], tmp$perm.conf.int[2],
                     tmp$perm.estimate, tmp$perm.stderr)
      }
      colnames(res) <- c("statistic", "pvalue", "R", "conf.low",
                        "conf.high", "mean.diff", "stderr")
    }
    if(!is.null(rx.H0) & !is.null(ry.H0)){
      if(parallel){
        data.comb.H0 <- cbind(data.x.H0, data.y.H0)
        res.H0 <- parRapply(cl = cl, x = data.comb.H0, FUN = row.perm.welch, 
                         nx = nx, ny = ny, conf.level = 1-sig.level, mu = mu, 
                         R = R, useCombn = useCombn)
      }else{
        res.H0 <- matrix(NA, nrow = iter, ncol = 7)
        for(i in 1:iter){
          tmp <- perm.t.test(x = data.x.H0[i,], y = data.y.H0[i,],
                             var.equal = FALSE, conf.level = 1-sig.level, 
                             conf.int = conf.int, mu = mu, 
                             alternative = alternative, R = R, 
                             useCombn = useCombn)
          res.H0[,i] <- c(tmp$statistic, tmp$perm.p.value, 
                          tmp$R.true, tmp$perm.conf.int[1], tmp$perm.conf.int[2],
                          tmp$perm.estimate, tmp$perm.stderr)
        }
        colnames(res.H0) <- c("statistic", "pvalue", "R", "conf.low",
                             "conf.high", "mean.diff", "stderr")
      }
      PERM.WELCH <- list("H1" = res, "H0" = res.H0)
    }else{
      PERM.WELCH <- list("H1" = res)  
    }
    RESULT <- c(RESULT, list(Perm.Welch = PERM.WELCH))
  }
  
  #########################################################
  ## Bootstrap Student 2-sample t-test
  #########################################################
  if("boot.student" %in% methods){
    if(parallel){
      row.boot.stud <- function(xy, nx, ny, conf.level, mu, alternative, R){
        tmp <- boot.t.test(x = xy[1:nx], 
                           y = xy[(nx+1):(nx+ny)], var.equal = TRUE, 
                           conf.level = conf.level, mu = mu, 
                           alternative = alternative, R = R)
        res <- c(tmp$statistic, tmp$boot.p.value, 
                 tmp$R, tmp$boot.conf.int[1], tmp$boot.conf.int[2],
                 tmp$boot.estimate, tmp$boot.stderr)
        names(res) <- c("statistic", "pvalue", "R", "conf.low",
                        "conf.high", "mean.diff", "stderr")
        res
      }
      data.comb <- cbind(data.x, data.y)
      res <- parRapply(cl = cl, x = data.comb, FUN = row.boot.stud, 
                       nx = nx, ny = ny, conf.level = 1-sig.level, mu = mu, 
                       R = R)
    }else{
      res <- matrix(NA, nrow = iter, ncol = 7)
      for(i in 1:iter){
        tmp <- boot.t.test(x = data.x[i,], y = data.y[i,], var.equal = TRUE, 
                           conf.level = 1-sig.level, conf.int = conf.int, mu = mu, 
                           alternative = alternative, R = R)
        res[,i] <- c(tmp$statistic, tmp$boot.p.value, 
                     tmp$R, tmp$boot.conf.int[1], tmp$boot.conf.int[2],
                     tmp$boot.estimate, tmp$boot.stderr)
      }
      colnames(res) <- c("statistic", "pvalue", "R", "conf.low",
                        "conf.high", "mean.diff", "stderr")
    }
    if(!is.null(rx.H0) & !is.null(ry.H0)){
      if(parallel){
        data.comb.H0 <- cbind(data.x.H0, data.y.H0)
        res.H0 <- parRapply(cl = cl, x = data.comb.H0, FUN = row.boot.stud, 
                            nx = nx, ny = ny, conf.level = 1-sig.level, mu = mu, 
                            R = R)
      }else{
        res.H0 <- matrix(NA, nrow = iter, ncol = 7)
        for(i in 1:iter){
          tmp <- boot.t.test(x = data.x.H0[i,], y = data.y.H0[i,],
                             var.equal = TRUE, conf.level = 1-sig.level, 
                             conf.int = conf.int, mu = mu, 
                             alternative = alternative, R = R)
          res.H0[,i] <- c(tmp$statistic, tmp$boot.p.value, 
                          tmp$R, tmp$boot.conf.int[1], tmp$boot.conf.int[2],
                          tmp$boot.estimate, tmp$boot.stderr)
        }
        colnames(res.H0) <- c("statistic", "pvalue", "R", "conf.low",
                             "conf.high", "mean.diff", "stderr")
      }
      BOOT.STUDENT <- list("H1" = res, "H0" = res.H0)
    }else{
      BOOT.STUDENT <- list("H1" = res)  
    }
    RESULT <- c(RESULT, list(Boot.Student = BOOT.STUDENT))
  }
  
  #########################################################
  ## Bootstrap Welch 2-sample t-test
  #########################################################
  if("boot.welch" %in% methods){
    if(parallel){
      row.boot.welch <- function(xy, nx, ny, conf.level, mu, alternative, R){
        tmp <- boot.t.test(x = xy[1:nx], 
                           y = xy[(nx+1):(nx+ny)], var.equal = FALSE, 
                           conf.level = conf.level, mu = mu, 
                           alternative = alternative, R = R)
        res <- c(tmp$statistic, tmp$boot.p.value, 
                 tmp$R, tmp$boot.conf.int[1], tmp$boot.conf.int[2],
                 tmp$boot.estimate, tmp$boot.stderr)
        names(res) <- c("statistic", "pvalue", "R", "conf.low",
                        "conf.high", "mean.diff", "stderr")
        res
      }
      data.comb <- cbind(data.x, data.y)
      res <- parRapply(cl = cl, x = data.comb, FUN = row.boot.welch, 
                       nx = nx, ny = ny, conf.level = 1-sig.level, mu = mu, 
                       R = R)
    }else{
      res <- matrix(NA, nrow = iter, ncol = 7)
      for(i in 1:iter){
        tmp <- boot.t.test(x = data.x[i,], y = data.y[i,], var.equal = FALSE, 
                           conf.level = 1-sig.level, conf.int = conf.int, mu = mu, 
                           alternative = alternative, R = R)
        res[,i] <- c(tmp$statistic, tmp$boot.p.value, 
                     tmp$R, tmp$boot.conf.int[1], tmp$boot.conf.int[2],
                     tmp$boot.estimate, tmp$boot.stderr)
      }
      colnames(res) <- c("statistic", "pvalue", "R", "conf.low",
                        "conf.high", "mean.diff", "stderr")
    }
    if(!is.null(rx.H0) & !is.null(ry.H0)){
      if(parallel){
        data.comb.H0 <- cbind(data.x.H0, data.y.H0)
        res.H0 <- parRapply(cl = cl, x = data.comb.H0, FUN = row.boot.welch, 
                            nx = nx, ny = ny, conf.level = 1-sig.level, mu = mu, 
                            R = R)
      }else{
        res.H0 <- matrix(NA, nrow = iter, ncol = 7)
        for(i in 1:iter){
          tmp <- boot.t.test(x = data.x.H0[i,], y = data.y.H0[i,],
                             var.equal = FALSE, conf.level = 1-sig.level, 
                             conf.int = conf.int, mu = mu, 
                             alternative = alternative, R = R)
          res.H0[,i] <- c(tmp$statistic, tmp$boot.p.value, 
                          tmp$R, tmp$boot.conf.int[1], tmp$boot.conf.int[2],
                          tmp$boot.estimate, tmp$boot.stderr)
        }
        colnames(res.H0) <- c("statistic", "pvalue", "R", "conf.low",
                             "conf.high", "mean.diff", "stderr")
      }
      BOOT.WELCH <- list("H1" = res, "H0" = res.H0)
    }else{
      BOOT.WELCH <- list("H1" = res)  
    }
    RESULT <- c(RESULT, list(Boot.Welch = BOOT.WELCH))
  }
  if(parallel && PARA){
    if(STOP) stopCluster(cl)
    if(DETACH) detach("package:MKinfer")
  }
  SetUp <- c("nx" = nx, "rx" = rx, "rx.H0" = rx.H0, 
             "ny" = ny, "ry" = ry, "ry.H0" = ry.H0, 
             "sig.level" = sig.level, "mu" = mu,
             "alternative" = alternative, "iter" = iter)
  RESULT <- c(RESULT, list(SetUp = SetUp))
  class(RESULT) <- "sim.power.ttest"
  RESULT
}

print.sim.power.ttest <- function(x, digits = getOption("digits"), ...){
  cat("\n    Simulation Set-up\n")
  alpha <- x$SetUp$sig.level
  y0 <- x$SetUp
  iter <- x$SetUp$iter
  cat(paste(format(names(y0), width = 15L, justify = "right"), 
            format(y0, digits = digits), sep = " = "), sep = "\n")
  RES <- NULL
  if("Student" %in% names(x)){
    cat("\n    Student Two-sample t-Test\n")
    y1 <- c("emp.power" = sum(x$Student$H1$pvalue < alpha)/iter)
    if(!is.null(x$Student$H0)){
      y1 <- c(y1, "emp.type.I.error" = sum(x$Student$H0$pvalue < alpha)/iter)
    }
    cat(paste(format(names(y1), width = 15L, justify = "right"), 
              format(y1, digits = digits), sep = " = "), sep = "\n")
    RES <- c(RES, y1)
  }
  if("Welch" %in% names(x)){
    cat("\n    Welch Two-sample t-Test\n")
    y2 <- c("emp.power" = sum(x$Welch$H1$pvalue < alpha)/iter)
    if(!is.null(x$Welch$H0)){
      y2 <- c(y2, "emp.type.I.error" = sum(x$Welch$H0$pvalue < alpha)/iter)
    }
    cat(paste(format(names(y2), width = 15L, justify = "right"), 
              format(y2, digits = digits), sep = " = "), sep = "\n")
    RES <- c(RES, y2)
  }
  if("Hsu" %in% names(x)){
    cat("\n    Hsu Two-sample t-Test\n")
    y3 <- c("emp.power" = sum(x$Hsu$H1$pvalue < alpha)/iter)
    if(!is.null(x$Hsu$H0)){
      y3 <- c(y3, "emp.type.I.error" = sum(x$Hsu$H0$pvalue < alpha)/iter)
    }
    cat(paste(format(names(y3), width = 15L, justify = "right"), 
              format(y3, digits = digits), sep = " = "), sep = "\n")
    RES <- c(RES, y3)
  }
  if("Xiao" %in% names(x)){
    cat("\n    Xiao Two-sample t-Test\n")
    y4 <- c("emp.power" = sum(x$Xiao$H1$pvalue < alpha)/iter)
    if(!is.null(x$Xiao$H0)){
      y4 <- c(y4, "emp.type.I.error" = sum(x$Xiao$H0$pvalue < alpha)/iter)
    }
    cat(paste(format(names(y4), width = 15L, justify = "right"), 
              format(y4, digits = digits), sep = " = "), sep = "\n")
    RES <- c(RES, y4)
  }
  if("Perm.Student" %in% names(x)){
    cat("\n    Permutation Student Two-sample t-Test\n")
    y5 <- c("emp.power" = sum(x$Perm.Student$H1$pvalue < alpha)/iter)
    if(!is.null(x$Perm.Student$H0)){
      y5 <- c(y5, "emp.type.I.error" = sum(x$Perm.Student$H0$pvalue < alpha)/iter)
    }
    cat(paste(format(names(y5), width = 15L, justify = "right"), 
              format(y5, digits = digits), sep = " = "), sep = "\n")
    RES <- c(RES, y5)
  }
  if("Perm.Welch" %in% names(x)){
    cat("\n    Permutation Welch Two-sample t-Test\n")
    y6 <- c("emp.power" = sum(x$Perm.Welch$H1$pvalue < alpha)/iter)
    if(!is.null(x$Perm.Welch$H0)){
      y6 <- c(y6, "emp.type.I.error" = sum(x$Perm.Welch$H0$pvalue < alpha)/iter)
    }
    cat(paste(format(names(y6), width = 15L, justify = "right"), 
              format(y6, digits = digits), sep = " = "), sep = "\n")
    RES <- c(RES, y6)
  }
  if("Boot.Student" %in% names(x)){
    cat("\n    Bootstrap Student Two-sample t-Test\n")
    y7 <- c("emp.power" = sum(x$Boot.Student$H1$pvalue < alpha)/iter)
    if(!is.null(x$Boot.Student$H0)){
      y7 <- c(y7, "emp.type.I.error" = sum(x$Boot.Student$H0$pvalue < alpha)/iter)
    }
    cat(paste(format(names(y7), width = 15L, justify = "right"), 
              format(y7, digits = digits), sep = " = "), sep = "\n")
    RES <- c(RES, y7)
  }
  if("Boot.Welch" %in% names(x)){
    cat("\n    Bootstrap Welch Two-sample t-Test\n")
    y8 <- c("emp.power" = sum(x$Boot.Welch$H1$pvalue < alpha)/iter)
    if(!is.null(x$Boot.Welch$H0)){
      y8 <- c(y8, "emp.type.I.error" = sum(x$Boot.Welch$H0$pvalue < alpha)/iter)
    }
    cat(paste(format(names(y8), width = 15L, justify = "right"), 
              format(y8, digits = digits), sep = " = "), sep = "\n")
    RES <- c(RES, y8)
  }
  cat("\n")
  invisible(RES)
}
