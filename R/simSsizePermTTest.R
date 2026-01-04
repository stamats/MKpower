sim.ssize.perm.t.test <- function(rx, ry = NULL, mu = 0, sig.level = 0.05, power = 0.8, 
                                  type = c("two.sample", "one.sample", "paired"), 
                                  alternative = c("two.sided", "less", "greater"),
                                  var.equal = FALSE, R = 9999, symmetric = TRUE,
                                  useCombn = FALSE, n.min = 10, n.max = 200, 
                                  step.size = 10, iter = 10000, BREAK = TRUE,
                                  parallel = FALSE, cl = NULL){
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
  
  if(parallel){
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
  
  if(type == "two.sample"){
    stopifnot(!is.null(ry))
    stopifnot(is.function(ry))
    
    rowtest2 <- function(XY, n, alternative, mu, var.equal, conf.level, R, 
                        symmetric, useCombn){
      perm.t.test(XY[1:n], XY[(n+1):(2*n)], alternative = alternative,
                  mu = mu, var.equal = var.equal, conf.level = conf.level, 
                  R = R, symmetric = symmetric, useCombn = useCombn)$perm.p.value
    }
    
    ns <- seq.int(from = n.min, to = n.max, by = step.size)
    empPower <- numeric(length(ns))
    for(i in seq_len(length(ns))){
      data.x <- matrix(rx(ns[i]*iter), nrow = iter)
      data.y <- matrix(ry(ns[i]*iter), nrow = iter)
      data.xy <- cbind(data.x, data.y)
      if(parallel){
        res <- parRapply(cl = cl, x = data.xy, FUN = rowtest2, 
                         n = ns[i], alternative = alternative, 
                         mu = mu, var.equal = var.equal, conf.level = 1-sig.level, 
                         R = R, symmetric = symmetric, useCombn = useCombn)
      }else{
        res <- apply(data.xy, 1, rowtest2, n = ns[i], alternative = alternative, 
                     mu = mu, var.equal = var.equal, conf.level = 1-sig.level, 
                     R = R, symmetric = symmetric, useCombn = useCombn)
      }
      empPower[i] <- sum(res < sig.level)/iter
      if(empPower[i] > power) if(BREAK) break
    }
    if(var.equal){
      METHOD <- "Permutation Student t test"
    }else{
      METHOD <- "Permutation Welch t test"
    }
  }
  if(type == "one.sample"){
    
    rowtest1 <- function(X, alternative, mu, conf.level, R){
      perm.t.test(X, alternative = alternative, mu = mu, 
                  conf.level = conf.level, R = R)$perm.p.value
    }
    
    ns <- seq.int(from = n.min, to = n.max, by = step.size)
    empPower <- numeric(length(ns))
    for(i in seq_len(length(ns))){
      data.x <- matrix(rx(ns[i]*iter), nrow = iter)
      if(parallel){
        res <- parRapply(cl, data.x, FUN = rowtest1, alternative = alternative, 
                         mu = mu, conf.level = 1-sig.level, R = R)
        
      }else{
        res <- apply(data.x, 1, rowtest1, alternative = alternative, mu = mu, 
                     conf.level = 1-sig.level, R = R)
      }
      empPower[i] <- sum(res < sig.level)/iter
      if(empPower[i] > power) if(BREAK) break
    }
    METHOD <- "One-sample t test"
  }
  if(type == "paired"){
    ns <- seq.int(from = n.min, to = n.max, by = step.size)
    empPower <- numeric(length(ns))
    for(i in seq_len(length(ns))){
      data.xy <- matrix(rx(ns[i]*iter), nrow = iter)
      if(parallel){
        res <- parRapply(cl, data.xy, FUN = rowtest1, alternative = alternative, 
                         mu = mu, conf.level = 1-sig.level, R = R)
      }else{
        res <- apply(data.xy, 1, rowtest1, alternative = alternative, mu = mu, 
                     conf.level = 1-sig.level, R = R)
      }
      empPower[i] <- sum(res < sig.level)/iter
      if(empPower[i] > power) if(BREAK) break
    }
    METHOD <- "Paired t test"
  }
  empPower[empPower == 0] <- NA
  names(empPower) <- ns
  if(parallel){
    if(STOP) stopCluster(cl)
    if(DETACH) detach("package:MKinfer")
  }
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
