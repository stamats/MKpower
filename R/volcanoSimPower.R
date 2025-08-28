volcano.sim.power.ttest <- function(x, alpha = 1, shape = 19, 
                                    hex = FALSE, bins = 50, ...){
  iter <- x$SetUp$iter
  mu <- x$SetUp$mu
  sig.level <- x$SetUp$sig.level
  DF <- data.frame(pvalue = c(x$Classical$H1$pvalue,
                              x$Welch$H1$pvalue,
                              x$Hsu$H1$pvalue),
                   MD = c(x$Classical$H1$mean.diff,
                          x$Welch$H1$mean.diff,
                          x$Hsu$H1$mean.diff),
                   test = c(rep("Classical two-sample t-test", iter),
                            rep("Welch two-sample t-test", iter),
                            rep("Hsu two-sample t-test", iter)),
                   hypothesis = rep("H1", 3*iter))
  if(!is.null(x$Classical$H0)){
    DF1 <- data.frame(pvalue = c(x$Classical$H0$pvalue,
                                 x$Welch$H0$pvalue,
                                 x$Hsu$H0$pvalue),
                      MD = c(x$Classical$H0$mean.diff,
                             x$Welch$H0$mean.diff,
                             x$Hsu$H0$mean.diff),
                      test = c(rep("Classical two-sample t-test", iter),
                               rep("Welch two-sample t-test", iter),
                               rep("Hsu two-sample t-test", iter)),
                      hypothesis = rep("H0", 3*iter))
    DF <- rbind(DF, DF1)
  }
  DF$test <- factor(DF$test, levels = c("Classical two-sample t-test",
                                        "Welch two-sample t-test",
                                        "Hsu two-sample t-test"))
  if(!is.null(x$Classical$H0)){
    Lab <- round(c(sum(x$Classical$H1$pvalue < sig.level)/nrow(x$Classical$H1),
                   sum(x$Welch$H1$pvalue < sig.level)/nrow(x$Welch$H1),
                   sum(x$Hsu$H1$pvalue < sig.level)/nrow(x$Hsu$H1),
                   sum(x$Classical$H0$pvalue < sig.level)/nrow(x$Classical$H0),
                   sum(x$Welch$H0$pvalue < sig.level)/nrow(x$Welch$H0),
                   sum(x$Hsu$H0$pvalue < sig.level)/nrow(x$Hsu$H0)), 4)
    Lab[1:3] <- paste("emp. power:", Lab[1:3])
    Lab[4:6] <- paste("emp. type-I-error:", Lab[4:6])
    DF.text <- data.frame(test = rep(c("Classical two-sample t-test",
                                       "Welch two-sample t-test",
                                       "Hsu two-sample t-test"), 2),
                          hypothesis = c(rep("H1", 3), rep("H0", 3)),
                          label = Lab)
    if(hex){
      gg <- ggplot(DF, aes(x = .data$MD, y = .data$pvalue)) + 
        geom_hex(bins = bins) + scale_y_neglog10() + 
        geom_line(data = data.frame(x = c(mu, mu), y = c(min(DF$pvalue)/10, max(DF$pvalue))),
                  aes(x = .data$x, y = .data$y)) +
        xlab("mean difference") + ylab("-log10(p value)") +
        geom_text(data = DF.text, aes(x = mu, y = min(DF$pvalue)/100, label = .data$label), 
                  vjust = 2, inherit.aes = FALSE) +
        geom_hline(yintercept = sig.level) +
        facet_grid(hypothesis ~ test, scales = "free_y")
    }else{
      gg <- ggplot(DF, aes(x = .data$MD, y = .data$pvalue)) + 
        geom_point(alpha = alpha, shape = shape) + scale_y_neglog10() + 
        geom_line(data = data.frame(x = c(mu, mu), y = c(min(DF$pvalue)/10, max(DF$pvalue))),
                  aes(x = .data$x, y = .data$y)) +
        xlab("mean difference") + ylab("-log10(p value)") +
        geom_text(data = DF.text, aes(x = mu, y = min(DF$pvalue)/100, label = .data$label), 
                  vjust = 2, inherit.aes = FALSE) +
        geom_hline(yintercept = sig.level) +
        facet_grid(hypothesis ~ test, scales = "free_y")
    }
  }else{
    Lab <- round(c(sum(x$Classical$H1$pvalue < sig.level)/nrow(x$Classical$H1),
                   sum(x$Welch$H1$pvalue < sig.level)/nrow(x$Welch$H1),
                   sum(x$Hsu$H1$pvalue < sig.level)/nrow(x$Hsu$H1)), 4)
    Lab <- paste("emp. power:", Lab)
    DF.text <- data.frame(test = c("Classical two-sample t-test",
                                   "Welch two-sample t-test",
                                   "Hsu two-sample t-test"),
                          hypothesis = rep("H1", 3),
                          label = Lab)
    if(hex){
      gg <- ggplot(DF, aes(x = .data$MD, y = .data$pvalue)) + 
        geom_hex(bins = bins) + scale_y_neglog10()+
        geom_vline(xintercept = mu) + 
        xlab("mean difference") + ylab("-log10(p value)") +
        geom_text(data = DF.text, aes(x = mu, y = Inf, label = .data$label), 
                  vjust = 2, inherit.aes = FALSE) +
        geom_hline(yintercept = sig.level) +
        facet_grid(~ test)
    }else{
      gg <- ggplot(DF, aes(x = .data$MD, y = .data$pvalue)) + 
        geom_point(alpha = alpha, shape = shape) + scale_y_neglog10()+
        geom_vline(xintercept = mu) + 
        xlab("mean difference") + ylab("-log10(p value)") +
        geom_text(data = DF.text, aes(x = mu, y = Inf, label = .data$label), 
                  vjust = 2, inherit.aes = FALSE) +
        geom_hline(yintercept = sig.level) +
        facet_grid(~ test)
    }
  }
  gg
}

volcano.sim.power.wtest <- function(x, alpha = 1, shape = 19, 
                                    hex = FALSE, bins = 50, ...){
  if(!x$SetUp$conf.int)
    stop("You have to run sim.power.wilcox.test with 'conf.int = TRUE'")
  iter <- x$SetUp$iter
  mu <- x$SetUp$mu
  sig.level <- x$SetUp$sig.level
  approximate <- x$SetUp$approximate
  if(approximate){
    DF <- data.frame(pvalue = c(x$Exact$H1$pvalue,
                                x$Asymptotic$H1$pvalue,
                                x$Approximate$H1$pvalue),
                     MD = c(x$Exact$H1$loc.diff,
                            x$Asymptotic$H1$loc.diff,
                            x$Approximate$H1$loc.diff),
                     test = c(rep("Exact Wilcoxon-Mann-Whitney test", iter),
                              rep("Asymptotic Wilcoxon-Mann-Whitney test", iter),
                              rep("Approximate Wilcoxon-Mann-Whitney test", iter)),
                     hypothesis = rep("H1", 3*iter))
  }else{
    DF <- data.frame(pvalue = c(x$Exact$H1$pvalue,
                                x$Asymptotic$H1$pvalue),
                     MD = c(x$Exact$H1$loc.diff,
                            x$Asymptotic$H1$loc.diff),
                     test = c(rep("Exact Wilcoxon-Mann-Whitney test", iter),
                              rep("Asymptotic Wilcoxon-Mann-Whitney test", iter)),
                     hypothesis = rep("H1", 2*iter))
  }
  if(!is.null(x$Classical$H0)){
    if(approximate){
      DF1 <- data.frame(pvalue = c(x$Exact$H0$pvalue,
                                   x$Asymptotic$H0$pvalue,
                                   x$Approximate$H0$pvalue),
                        MD = c(x$Exact$H0$loc.diff,
                               x$Asymptotic$H0$loc.diff,
                               x$Approximate$H0$loc.diff),
                        test = c(rep("Exact Wilcoxon-Mann-Whitney test", iter),
                                 rep("Asymptotic Wilcoxon-Mann-Whitney test", iter),
                                 rep("Approximate Wilcoxon-Mann-Whitney test", iter)),
                        hypothesis = rep("H0", 3*iter))
    }else{
      DF1 <- data.frame(pvalue = c(x$Exact$H0$pvalue,
                                   x$Asymptotic$H0$pvalue),
                        MD = c(x$Exact$H0$loc.diff,
                               x$Asymptotic$H0$loc.diff),
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
  if(!is.null(x$Classical$H0)){
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
    if(hex){
      gg <- ggplot(DF, aes(x = .data$MD, y = .data$pvalue)) + 
        geom_hex(bins = bins) + scale_y_neglog10() + 
        geom_line(data = data.frame(x = c(mu, mu), y = c(min(DF$pvalue)/10, max(DF$pvalue))),
                  aes(x = .data$x, y = .data$y)) +
        xlab("mean difference") + ylab("-log10(p value)") +
        geom_text(data = DF.text, aes(x = mu, y = min(DF$pvalue)/100, label = .data$label), 
                  vjust = 2, inherit.aes = FALSE) +
        geom_hline(yintercept = sig.level) +
        facet_grid(hypothesis ~ test, scales = "free_y")
    }else{
      gg <- ggplot(DF, aes(x = .data$MD, y = .data$pvalue)) + 
        geom_point(alpha = alpha, shape = shape) + scale_y_neglog10() + 
        geom_line(data = data.frame(x = c(mu, mu), y = c(min(DF$pvalue)/10, max(DF$pvalue))),
                  aes(x = .data$x, y = .data$y)) +
        xlab("mean difference") + ylab("-log10(p value)") +
        geom_text(data = DF.text, aes(x = mu, y = min(DF$pvalue)/100, label = .data$label), 
                  vjust = 2, inherit.aes = FALSE) +
        geom_hline(yintercept = sig.level) +
        facet_grid(hypothesis ~ test, scales = "free_y")
    }
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
    if(hex){
      gg <- ggplot(DF, aes(x = .data$MD, y = .data$pvalue)) + 
        geom_hex(bins = bins) + scale_y_neglog10()+
        geom_vline(xintercept = mu) + 
        xlab("mean difference") + ylab("-log10(p value)") +
        geom_text(data = DF.text, aes(x = mu, y = Inf, label = .data$label), 
                  vjust = 2, inherit.aes = FALSE) +
        geom_hline(yintercept = sig.level) +
        facet_grid(~ test)
    }else{
      gg <- ggplot(DF, aes(x = .data$MD, y = .data$pvalue)) + 
        geom_point(alpha = alpha, shape = shape) + scale_y_neglog10()+
        geom_vline(xintercept = mu) + 
        xlab("mean difference") + ylab("-log10(p value)") +
        geom_text(data = DF.text, aes(x = mu, y = Inf, label = .data$label), 
                  vjust = 2, inherit.aes = FALSE) +
        geom_hline(yintercept = sig.level) +
        facet_grid(~ test)
    }
  }
  gg
}
