ggplot(marker.region.data) +
  geom_histogram(aes(rep.meth.fraction, y = ..density..), colour="black", fill="white", binwidth = .005) +
  stat_function(fun = function(x) dbeta(x, alpha0, beta0), color = "red",
                size = 1) +
  stat_function(fun = function(x) dbeta(x, alpha1, beta1), color = "blue",
                size = 1) +
  stat_function(fun = function(x) dbeta(x, alpha2, beta2), color = "purple",
                size = 1) +
  xlab("Methylation Fraction")

df <- data.frame(rep.meth.fraction = rep.meth.fraction)


ggplot(data = df, aes(x = rep.meth.fraction)) +
  geom_histogram() +
  stat_function(fun = function(x) dbeta(x, fit.emp[1], fit.emp[2]), color = "red") +
  stat_function(fun = function(x) dbeta(x, fit.optim[1], fit.optim[2]), color = "blue") +
  xlim(0, 1.1)


ks.test(x = rbeta(n = 1000, fit.emp[1], fit.emp[2]), y = rbeta(n = 1000, fit.optim[1], fit.optim[2]))