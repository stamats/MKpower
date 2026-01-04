hist.sim.power.ttest <- function(x, color.hline = "orange", ...){
  iter <- x$SetUp$iter
  alpha <- x$SetUp$sig.level
  
  DF1 <- data.frame(NULL)
  DF0 <- data.frame(NULL)
  
  if(!is.null(x$SetUp$rx) && !is.null(x$SetUp$ry)){
    if("Student" %in% names(x)){
      tmp <- data.frame(pvalue = x$Student$H1$pvalue,
                        test = rep("Student two-sample t-test", iter),
                        hypothesis = rep("H1", iter))
      DF1 <- rbind(DF1, tmp)
    }
    if("Welch" %in% names(x)){
      tmp <- data.frame(pvalue = x$Welch$H1$pvalue,
                        test = rep("Welch two-sample t-test", iter),
                        hypothesis = rep("H1", iter))
      DF1 <- rbind(DF1, tmp)
    }
    if("Hsu" %in% names(x)){
      tmp <- data.frame(pvalue = x$Hsu$H1$pvalue,
                        test = rep("Hsu two-sample t-test", iter),
                        hypothesis = rep("H1", iter))
      DF1 <- rbind(DF1, tmp)
    }
    if("Xiao" %in% names(x)){
      tmp <- data.frame(pvalue = x$Xiao$H1$pvalue,
                        test = rep("Xiao two-sample t-test", iter),
                        hypothesis = rep("H1", iter))
      DF1 <- rbind(DF1, tmp)
    }
    if("Perm.Student" %in% names(x)){
      tmp <- data.frame(pvalue = x$Perm.Student$H1$pvalue,
                        test = rep("Permutation Student two-sample t-test", iter),
                        hypothesis = rep("H1", iter))
      DF1 <- rbind(DF1, tmp)
    }
    if("Perm.Welch" %in% names(x)){
      tmp <- data.frame(pvalue = x$Perm.Welch$H1$pvalue,
                        test = rep("Permutation Welch two-sample t-test", iter),
                        hypothesis = rep("H1", iter))
      DF1 <- rbind(DF1, tmp)
    }
    if("Boot.Student" %in% names(x)){
      tmp <- data.frame(pvalue = x$Boot.Student$H1$pvalue,
                        test = rep("Bootstrap Student two-sample t-test", iter),
                        hypothesis = rep("H1", iter))
      DF1 <- rbind(DF1, tmp)
    }
    if("Boot.Welch" %in% names(x)){
      tmp <- data.frame(pvalue = x$Boot.Welch$H1$pvalue,
                        test = rep("Bootstrap Welch two-sample t-test", iter),
                        hypothesis = rep("H1", iter))
      DF1 <- rbind(DF1, tmp)
    }
    DF <- DF1
  }
  
  if(!is.null(x$SetUp$rx.H0) && !is.null(x$SetUp$ry.H0)){
    if("Student" %in% names(x)){
      tmp <- data.frame(pvalue = x$Student$H0$pvalue,
                        test = rep("Student two-sample t-test", iter),
                        hypothesis = rep("H0", iter))
      DF0 <- rbind(DF0, tmp)
    }
    if("Welch" %in% names(x)){
      tmp <- data.frame(pvalue = x$Welch$H0$pvalue,
                        test = rep("Welch two-sample t-test", iter),
                        hypothesis = rep("H0", iter))
      DF0 <- rbind(DF0, tmp)
    }
    if("Hsu" %in% names(x)){
      tmp <- data.frame(pvalue = x$Hsu$H0$pvalue,
                        test = rep("Hsu two-sample t-test", iter),
                        hypothesis = rep("H0", iter))
      DF0 <- rbind(DF0, tmp)
    }
    if("Xiao" %in% names(x)){
      tmp <- data.frame(pvalue = x$Xiao$H0$pvalue,
                        test = rep("Xiao two-sample t-test", iter),
                        hypothesis = rep("H0", iter))
      DF0 <- rbind(DF0, tmp)
    }
    if("Perm.Student" %in% names(x)){
      tmp <- data.frame(pvalue = x$Perm.Student$H0$pvalue,
                        test = rep("Permutation Student two-sample t-test", iter),
                        hypothesis = rep("H0", iter))
      DF1 <- rbind(DF1, tmp)
    }
    if("Perm.Welch" %in% names(x)){
      tmp <- data.frame(pvalue = x$Perm.Welch$H0$pvalue,
                        test = rep("Permutation Welch two-sample t-test", iter),
                        hypothesis = rep("H0", iter))
      DF1 <- rbind(DF1, tmp)
    }
    if("Boot.Student" %in% names(x)){
      tmp <- data.frame(pvalue = x$Boot.Student$H0$pvalue,
                        test = rep("Bootstrap Student two-sample t-test", iter),
                        hypothesis = rep("H0", iter))
      DF1 <- rbind(DF1, tmp)
    }
    if("Boot.Welch" %in% names(x)){
      tmp <- data.frame(pvalue = x$Boot.Welch$H0$pvalue,
                        test = rep("Bootstrap Welch two-sample t-test", iter),
                        hypothesis = rep("H0", iter))
      DF1 <- rbind(DF1, tmp)
    }
    DF <- rbind(DF1, DF0)
  }
  methods <- c("Student two-sample t-test", "Welch two-sample t-test",
               "Hsu two-sample t-test", "Xiao two-sample t-test",
               "Permutation Student two-sample t-test",
               "Permutation Welch two-sample t-test",
               "Bootstrap Student two-sample t-test",
               "Bootstrap Welch two-sample t-test")
  selected <- unique(DF$test)
  DF$test <- factor(DF$test, levels = methods[methods %in% selected])
  
  Lab1 <- Lab0 <- Test <- NULL
  if("Student" %in% names(x)){
    if(!is.null(x$Student$H1)){
      Lab1 <- c(Lab1, sum(x$Student$H1$pvalue < alpha)/nrow(x$Student$H1))
    }
    if(!is.null(x$Student$H0)){
      Lab0 <- c(Lab0, sum(x$Student$H0$pvalue < alpha)/nrow(x$Student$H0))
    }
    Test <- c(Test, "Student two-sample t-test")
  }
  if("Welch" %in% names(x)){
    if(!is.null(x$Welch$H1)){
      Lab1 <- c(Lab1, sum(x$Welch$H1$pvalue < alpha)/nrow(x$Welch$H1))
    }
    if(!is.null(x$Welch$H0)){
      Lab0 <- c(Lab0, sum(x$Welch$H0$pvalue < alpha)/nrow(x$Welch$H0))
    }
    Test <- c(Test, "Welch two-sample t-test")
  }
  if("Hsu" %in% names(x)){
    if(!is.null(x$Hsu$H1)){
      Lab1 <- c(Lab1, sum(x$Hsu$H1$pvalue < alpha)/nrow(x$Hsu$H1))
    }
    if(!is.null(x$Hsu$H0)){
      Lab0 <- c(Lab0, sum(x$Hsu$H0$pvalue < alpha)/nrow(x$Hsu$H0))
    }
    Test <- c(Test, "Hsu two-sample t-test")
  }
  if("Xiao" %in% names(x)){
    if(!is.null(x$Xiao$H1)){
      Lab1 <- c(Lab1, sum(x$Xiao$H1$pvalue < alpha)/nrow(x$Xiao$H1))
    }
    if(!is.null(x$Hsu$H0)){
      Lab0 <- c(Lab0, sum(x$Xiao$H0$pvalue < alpha)/nrow(x$Xiao$H0))
    }
    Test <- c(Test, "Xiao two-sample t-test")
  }
  if("Perm.Student" %in% names(x)){
    if(!is.null(x$Perm.Student$H1)){
      Lab1 <- c(Lab1, sum(x$Perm.Student$H1$pvalue < alpha)/nrow(x$Perm.Student$H1))
    }
    if(!is.null(x$Perm.Student$H0)){
      Lab0 <- c(Lab0, sum(x$Perm.Student$H0$pvalue < alpha)/nrow(x$Perm.Student$H0))
    }
    Test <- c(Test, "Permutation Student two-sample t-test")
  }
  if("Perm.Welch" %in% names(x)){
    if(!is.null(x$Perm.Welch$H1)){
      Lab1 <- c(Lab1, sum(x$Perm.Welch$H1$pvalue < alpha)/nrow(x$Perm.Welch$H1))
    }
    if(!is.null(x$Perm.Welch$H0)){
      Lab0 <- c(Lab0, sum(x$Perm.Welch$H0$pvalue < alpha)/nrow(x$Perm.Welch$H0))
    }
    Test <- c(Test, "Permutation Welch two-sample t-test")
  }
  if("Boot.Student" %in% names(x)){
    if(!is.null(x$Boot.Student$H1)){
      Lab1 <- c(Lab1, sum(x$Boot.Student$H1$pvalue < alpha)/nrow(x$Boot.Student$H1))
    }
    if(!is.null(x$Boot.Student$H0)){
      Lab0 <- c(Lab0, sum(x$Boot.Student$H0$pvalue < alpha)/nrow(x$Boot.Student$H0))
    }
    Test <- c(Test, "Bootstrap Student two-sample t-test")
  }
  if("Boot.Welch" %in% names(x)){
    if(!is.null(x$Boot.Welch$H1)){
      Lab1 <- c(Lab1, sum(x$Boot.Welch$H1$pvalue < alpha)/nrow(x$Boot.Welch$H1))
    }
    if(!is.null(x$Boot.Welch$H0)){
      Lab0 <- c(Lab0, sum(x$Boot.Welch$H0$pvalue < alpha)/nrow(x$Boot.Welch$H0))
    }
    Test <- c(Test, "Bootstrap Welch two-sample t-test")
  }
  if(!is.null(Lab1)) Lab1 <- paste0("emp. power: ", Lab1)
  if(!is.null(Lab0)) Lab0 <- paste0("emp. type-I-error: ", Lab0)
  
  if(!is.null(Lab0) && !is.null(Lab1)){
    DF.text <- data.frame(test = rep(Test, 2),
                          hypothesis = c(rep("H1", length(Lab1)),
                                         rep("H0", length(Lab0))),
                          label = c(Lab1, Lab0))
    DF.text$test <- factor(DF.text$test, levels = methods[methods %in% selected])
    gg <- ggplot(data = DF, aes(x = .data$pvalue)) + 
      geom_histogram(aes(y = after_stat(density)), binwidth = 0.01) + 
      geom_hline(yintercept = 1.0, color = color.hline) + xlab("p value") +
      geom_text(data = DF.text, aes(x = 0.5, y = Inf, label = .data$label), 
                vjust = 2, inherit.aes = FALSE) +
      facet_grid(hypothesis ~ test, scales = "free_y")
  }
  if(is.null(Lab0) && !is.null(Lab1)){
    DF.text <- data.frame(test = Test,
                          hypothesis = rep("H1", length(Lab1)),
                          label = Lab1)
    DF.text$test <- factor(DF.text$test, levels = methods[methods %in% selected])
    gg <- ggplot(data = DF, aes(x = .data$pvalue)) + 
      geom_histogram(aes(y = after_stat(density)), binwidth = 0.01) + 
      geom_text(data = DF.text, aes(x = 0.5, y = Inf, label = .data$label), 
                vjust = 2, inherit.aes = FALSE) +
      geom_hline(yintercept = 1.0, color = color.hline) + xlab("p value") +
      facet_grid(~ test)
  }
  if(!is.null(Lab0) && is.null(Lab1)){
    DF.text <- data.frame(test = Test,
                          hypothesis = rep("H0", length(Lab0)),
                          label = Lab0)
    DF.text$test <- factor(DF.text$test, levels = methods[methods %in% selected])
    gg <- ggplot(data = DF, aes(x = .data$pvalue)) + 
      geom_histogram(aes(y = after_stat(density)), binwidth = 0.01) + 
      geom_text(data = DF.text, aes(x = 0.5, y = Inf, label = .data$label), 
                vjust = 2, inherit.aes = FALSE) +
      geom_hline(yintercept = 1.0, color = color.hline) + xlab("p value") +
      facet_grid(~ test)
  }
  gg
}

hist.sim.power.wtest <- function(x, color.hline = "orange", ...){
  iter <- x$SetUp$iter
  sig.level <- x$SetUp$sig.level
  approximate <- x$SetUp$approximate
  if(approximate){
    DF <- data.frame(pvalue = c(x$Exact$H1$pvalue,
                                x$Asymptotic$H1$pvalue,
                                x$Approximate$H1$pvalue),
                     test = c(rep("Exact Wilcoxon-Mann-Whitney test", iter),
                              rep("Asymptotic Wilcoxon-Mann-Whitney test", iter),
                              rep("Approximate Wilcoxon-Mann-Whitney test", iter)),
                     hypothesis = rep("H1", 3*iter))
  }else{
    DF <- data.frame(pvalue = c(x$Exact$H1$pvalue,
                                x$Asymptotic$H1$pvalue),
                     test = c(rep("Exact Wilcoxon-Mann-Whitney test", iter),
                              rep("Asymptotic Wilcoxon-Mann-Whitney test", iter)),
                     hypothesis = rep("H1", 2*iter))
  }
  if(!is.null(x$Exact$H0)){
    if(approximate){
      DF1 <- data.frame(pvalue = c(x$Exact$H0$pvalue,
                                   x$Asymptotic$H0$pvalue,
                                   x$Approximate$H0$pvalue),
                        test = c(rep("Exact Wilcoxon-Mann-Whitney test", iter),
                                 rep("Asymptotic Wilcoxon-Mann-Whitney test", iter),
                                 rep("Approximate Wilcoxon-Mann-Whitney test", iter)),
                        hypothesis = rep("H0", 3*iter))
    }else{
      DF1 <- data.frame(pvalue = c(x$Exact$H0$pvalue,
                                   x$Asymptotic$H0$pvalue),
                        test = c(rep("Exact Wilcoxon-Mann-Whitney test", iter),
                                 rep("Asymptotic Wilcoxon-Mann-Whitney test", iter)),
                        hypothesis = rep("H0", 2*iter))
    }
    DF <- rbind(DF, DF1)
  }
  if(approximate){
    DF$test <- factor(DF$test, levels = c("Exact Wilcoxon-Mann-Whitney test",
                                          "Asymptotic Wilcoxon-Mann-Whitney test",
                                          "Approximate Wilcoxon-Mann-Whitney test"))
  }else{
    DF$test <- factor(DF$test, levels = c("Exact Wilcoxon-Mann-Whitney test",
                                          "Asymptotic Wilcoxon-Mann-Whitney test"))
  }
  if(!is.null(x$Exact$H0)){
    if(approximate){
      Lab <- round(c(sum(x$Exact$H1$pvalue < sig.level)/iter,
                     sum(x$Asymptotic$H1$pvalue < sig.level)/iter,
                     sum(x$Approximate$H1$pvalue < sig.level)/iter,
                     sum(x$Exact$H0$pvalue < sig.level)/iter,
                     sum(x$Asymptotic$H0$pvalue < sig.level)/iter,
                     sum(x$Approximate$H0$pvalue < sig.level)/iter), 4)
      Lab[1:3] <- paste("emp. power:", Lab[1:3])
      Lab[4:6] <- paste("emp. type-I-error:", Lab[4:6])
      DF.text <- data.frame(test = rep(c("Exact Wilcoxon-Mann-Whitney test",
                                         "Asymptotic Wilcoxon-Mann-Whitney test",
                                         "Approximate Wilcoxon-Mann-Whitney test"), 2),
                            hypothesis = c(rep("H1", 3), rep("H0", 3)),
                            label = Lab)
    }else{
      Lab <- round(c(sum(x$Exact$H1$pvalue < sig.level)/iter,
                     sum(x$Asymptotic$H1$pvalue < sig.level)/iter,
                     sum(x$Exact$H0$pvalue < sig.level)/iter,
                     sum(x$Asymptotic$H0$pvalue < sig.level)/iter), 4)
      Lab[1:2] <- paste("emp. power:", Lab[1:2])
      Lab[3:4] <- paste("emp. type-I-error:", Lab[3:4])
      DF.text <- data.frame(test = rep(c("Exact Wilcoxon-Mann-Whitney test",
                                         "Asymptotic Wilcoxon-Mann-Whitney test"), 2),
                            hypothesis = c(rep("H1", 2), rep("H0", 2)),
                            label = Lab)
    }
    gg <- ggplot(data = DF, aes(x = .data$pvalue)) + 
      geom_histogram(aes(y = after_stat(density)), binwidth = 0.01) + 
      geom_hline(yintercept = 1.0, color = color.hline) + xlab("p value") +
      geom_text(data = DF.text, aes(x = 0.5, y = Inf, label = .data$label), 
                vjust = 2, inherit.aes = FALSE) +
      facet_grid(hypothesis ~ test, scales = "free_y")
  }else{
    if(approximate){
      Lab <- round(c(sum(x$Exact$H1$pvalue < sig.level)/iter,
                     sum(x$Asymptotic$H1$pvalue < sig.level)/iter,
                     sum(x$Approxmiate$H1$pvalue < sig.level)/iter), 4)
      Lab <- paste("emp. power:", Lab)
      DF.text <- data.frame(test = c("Exact Wilcoxon-Mann-Whitney test",
                                     "Asymptotic Wilcoxon-Mann-Whitney test",
                                     "Approximate Wilcoxon-Mann-Whitney test"),
                            hypothesis = rep("H1", 3),
                            label = Lab)
    }else{
      Lab <- round(c(sum(x$Exact$H1$pvalue < sig.level)/iter,
                     sum(x$Asymptotic$H1$pvalue < sig.level)/iter), 4)
      Lab <- paste("emp. power:", Lab)
      DF.text <- data.frame(test = c("Exact Wilcoxon-Mann-Whitney test",
                                     "Asymptotic Wilcoxon-Mann-Whitney test"),
                            hypothesis = rep("H1", 2),
                            label = Lab)
    }
    gg <- ggplot(data = DF, aes(x = .data$pvalue)) + 
      geom_histogram(aes(y = after_stat(density)), binwidth = 0.01) + 
      geom_text(data = DF.text, aes(x = 0.5, y = Inf, label = .data$label), 
                vjust = 2, inherit.aes = FALSE) +
      geom_hline(yintercept = 1.0, color = color.hline) + xlab("p value") +
      facet_grid(~ test)
  }
  gg
}
