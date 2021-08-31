library(tidyverse)
library(ggplot2)
library(binom)

###################################################################################

# excel solutions - for troubleshooting purposes
out.c.lambda <-  0.00189705440760631
out.c.tau<-  3.1891416476047
out.v.lambda <-  0.002237103
out.v.tau<-  1.658323567

###################################################################################

# Reporting the observed data from excel file (manually inputted)

# Ovarian category
ov.cat <- c(0:7) 

# Number of infected 
t.viv.inf <- c(0, 21, 58, 96, 260, 258, 214, 124)    # T. vivax 
t.con.inf <- c(0, 17, 44, 52, 188, 212, 157, 113)    # T. congolense 

# Sample size
t.viv.n <- c(0, 5133, 5334, 4284, 5452, 4758, 3209, 1768)    # T. vivax 
t.con.n <- c(0, 5131, 5333, 4283, 5451, 4758, 3209, 1768)    # T. congolense 

# Prevalence
prev.viv <- t.viv.inf/t.viv.n    # T. vivax
prev.con <- t.con.inf/t.con.n    # T.congolense

# Separate dataframes for species
viv.data <- data.frame(ov.cat, prev.viv)
colnames(viv.data) <- c("ovarian_cat", "observed_prev")
con.data <- data.frame(ov.cat, prev.con)
colnames(con.data) <- c("ovarian_cat", "observed_prev")

# Including species column to help when merging the dataframes & for plotting
viv.data$Species <- "vivax"
con.data$Species <- "congolense"

# Combine dataframes => easier to plot
viv.con.df <- rbind(viv.data, con.data)

# Remove NaNs from prevalence 
viv.con.df[is.na(viv.con.df)] <- 0
viv.data[is.na(viv.data)] <- 0
con.data[is.na(con.data)] <- 0

# Plot prevalences of T. vivax and T.congolense
axis.ov.cat <- c("0","1","2","3","4+4n","5+4n","6+4n","7+4n")

ggplot(viv.con.df, aes(x = ovarian_cat, y = observed_prev, 
                       group = Species, colour = Species))+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = seq(0,7,1), labels = axis.ov.cat)+
  theme_bw()+
  labs(x = "Ovarian Category", y = "Observed Prevalence", 
       title = "Prevalence by age and trypanosome species for female tsetse flies in Zimbabwe")+
  theme(plot.title = element_text(hjust = 0.5, size = 10))


###################################################################################

### FUNCTION: takes in parameters, spits out predicted

# Catalytic function

prev.func <- function(a=0:47, loglambda, logtau){
  lambda = exp(loglambda)
  tau = exp(logtau)
  prev = 1-exp(-lambda*(a-tau))
  prev[prev<0] = 0                         # only want positive prevalence values
  return(prev)
}

prev_vec <- prev.func(0:47, log(0.0019), log(3.189)) # testing function
prev_vec

###################################################################################

### FUNCTION: gives pooled (predicted) prevalence relative to ovarian categories

# Prevalence for every bloodmeal -> averaged over the 15 meals -> averaged over ovarian category

pooled.prev <- function(a=0:47, loglambda, logtau){
  # create a dataframe of prevalence, meals and ovarian age
  prev.df <- tibble(newprev = (prev.func(a, loglambda, logtau)),
                    meal = c(0:47),
                    ovarianAge = c(rep(0:15, each = 3)))
  
  # create a target df will contain ovarian categories with their relative pooled prevalences
  target.df <- prev.df %>%
    group_by(ovarianAge) %>%
    # find the mean of prevalences for the grouped meals in threes (0-15)
    summarize(ave_meal = mean(newprev)) %>%
    # change the ovarian age to ovarian category by pooling the ages 4+4n, 5+4n, 6+4n, 7+4n
    mutate(ovarianCat = ovarianAge) %>%
    mutate(ovarianCat = ifelse(ovarianCat %in% c(4,8,12), 4, ovarianCat)) %>%
    mutate(ovarianCat = ifelse(ovarianCat %in% c(5,9,13), 5, ovarianCat)) %>%
    mutate(ovarianCat = ifelse(ovarianCat %in% c(6,10,14), 6, ovarianCat)) %>%
    mutate(ovarianCat = ifelse(ovarianCat %in% c(7,11,15), 7, ovarianCat)) %>%
    group_by(ovarianCat) %>%
    # find the mean of prevalences for the pooled ovarian categories => our target
    summarize(predicted_prev = mean(ave_meal))
  
  return(target.df)
}

pooled.prev(loglambda = log(0.0019), logtau = log(3.189)) # test the function
pooled.prev(loglambda = log(0.00196), logtau = log(4)) # test the function again

###################################################################################

## FUNCTION: takes in predicted and observed prevalences and returns error

# T. congolense
RSS.c <- function(a = 0:47, loglambda, logtau){
  congo <- cbind(con.data, pooled.prev(a, loglambda, logtau))
  congo1 <- congo %>%
    mutate(RSS.c = (observed_prev-predicted_prev)^2)
  return(sum(congo1$RSS.c))
}

#T.vivax
RSS.v <- function(a = 0:47, loglambda, logtau){
  vivax <- cbind(viv.data, pooled.prev(a, loglambda, logtau))
  vivax1 <- vivax %>%
    mutate(RSS.v = (observed_prev-predicted_prev)^2)
  return(sum(vivax1$RSS.v))
}

RSS.c(loglambda = log(0.004441465), logtau = log(19.066947005))
RSS.v(loglambda = log(0.002235216), logtau = log(1.660081055))

###################################################################################

# Different values or parameters lambda and tau using RSS

# CONGOLENSE 

# ** Looking at what RSS.c looks like for different values of loglambda
loglambda.vec <-  seq(log(0.00001), log(0.01), length.out = 100)
plot(exp(loglambda.vec)
     , log(sapply(X=loglambda.vec, FUN = RSS.c, a = 0:47
                  , logtau = log(3.2))), type = 'l'
     , main = "RSS function values of lambda for T. congolense"
     , ylab = "RSS.c(lambda)"
     , xlab = "Different values of lambda")

abline(v = out.c.lambda, col = 'red')
# minimum is at the right place

# Looking at the contour plot in an "inelegant" way
con.n <- 20
loglambda.vec <-  seq(log(0.001), log(0.003), length.out = con.n)
loglambda.n <- length(loglambda.vec)
logtau.vec <-  seq(log(2), log(4), length.out = con.n)
logtau.n <- length(logtau.vec)
rss.c <- matrix(NA, nrow = loglambda.n, ncol = logtau.n)
# lambda by row, tau by column
for (ii in 1:loglambda.n){
  ii.loglambda <- loglambda.vec[ii]
  for (jj in 1:logtau.n){
    jj.logtau <- logtau.vec[jj]
    rss.c[ii,jj] <- RSS.c(loglambda = ii.loglambda, logtau = jj.logtau)
  }
}
# look at contour plot
contour(x = exp(loglambda.vec),
        y = exp(logtau.vec),
        log(rss.c), nlevels = 40, col = terrain.colors(50),
        plot.title = title(main = "Contour plot for RSS values of lambda and tau for T.congolense", 
                           xlab='lambda',  ylab='tau'))
points(out.c.lambda, out.c.tau, pch = 20)

# VIVAX

# ** Looking at what RSS.v looks like for different values of loglambda
loglambda.vec <-  seq(log(0.00001), log(0.01), length.out = 100)
plot(exp(loglambda.vec)
     , log(sapply(X=loglambda.vec, FUN = RSS.v, a = 0:47
                  , logtau = log(1.7))), type = 'l'
     , main = "RSS function values of lambda for T. vivax"
     , ylab = "RSS.v(lambda)"
     , xlab = "Different values of lambda")

abline(v = out.v.lambda, col = 'red')
# minimum is at the right place

# Looking at the contour plot in an "inelegant" way
viv.n <- 20
loglambda.vec <-  seq(log(0.001), log(0.003), length.out = viv.n)
loglambda.n <- length(loglambda.vec)
logtau.vec <-  seq(log(1), log(3), length.out = viv.n)
logtau.n <- length(logtau.vec)
rss.v <- matrix(NA, nrow = loglambda.n, ncol = logtau.n)

# lambda by row, tau by column
for (ii in 1:loglambda.n){
  ii.loglambda <- loglambda.vec[ii]
  for (jj in 1:logtau.n){
    jj.logtau <- logtau.vec[jj]
    rss.v[ii,jj] <- RSS.v(loglambda = ii.loglambda, logtau = jj.logtau)
  }
}

# look at contour plot
contour(x = exp(loglambda.vec),
        y = exp(logtau.vec),
        log(rss.v), nlevels = 40, col = terrain.colors(50),
        plot.title = title(main = "Contour plot for RSS values of lambda and tau for T. vivax", 
                           xlab='lambda',  ylab='tau'))
points(out.v.lambda, out.v.tau, pch = 20)

###################################################################################

# FUNCTION: minimize RSS and plot
# ?optim

# CONGOLENSE

trace <- 3

params1 <- c(loglambda = log(0.0001), logtau = log(5))     # initial parameter guesses
params2 <- c(loglambda = log(0.002), logtau = log(3.2))    # initial parameter guesses
params3 <- c(loglambda = log(0.001), logtau = log(4))      # initial parameter guesses
params4 <- c(loglambda = log(0.002), logtau = log(2))      # initial parameter guesses

opt.con.logRSS <- optim(par = params1,
                        fn = function(x) 10*log(RSS.c(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.con.logRSS <- exp(opt.con.logRSS$par)
opt.con.logRSS1 <- opt.con.logRSS
opt.con.logRSS <- optim(par = params2,
                        fn = function(x) 10*log(RSS.c(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.con.logRSS <- exp(opt.con.logRSS$par)
opt.con.logRSS2 <- opt.con.logRSS
opt.con.logRSS <- optim(par = params3,
                        fn = function(x) 10*log(RSS.c(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.con.logRSS <- exp(opt.con.logRSS$par)
opt.con.logRSS3 <- opt.con.logRSS
opt.con.logRSS <- optim(par = params4,
                        fn = function(x) 10*log(RSS.c(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.con.logRSS <- exp(opt.con.logRSS$par)
opt.con.logRSS4 <- opt.con.logRSS
opt.con.logRSS1
opt.con.logRSS2
opt.con.logRSS3
opt.con.logRSS4

# Stable and reasonable, plot to see differece between predicted and observed
# plot with confidence intervals

ci.con.rss <- binom.confint(x = t.con.inf, n = t.con.n, conf.level = 0.95, methods = "exact")

con.obs.pred.rss <- cbind(con.data, pooled.prev(loglambda = log(opt.con.logRSS1[1])
                      , logtau = log(opt.con.logRSS1[2]))
                      , lower = ci.con.rss$lower, upper=ci.con.rss$upper)
con.obs.pred.rss$observed_prev[1] <- con.obs.pred.rss$lower[1] <- con.obs.pred.rss$upper[1]  <- 0



ggplot(con.obs.pred.rss, aes(x = ovarian_cat, y = observed_prev, colour = "blue"))+
  geom_ribbon(aes(ymin = lower, ymax=upper), alpha = 0.3, colour = NA)+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = seq(0,7,1), labels = axis.ov.cat)+
  geom_point(aes(x = ovarian_cat, y = predicted_prev, colour = "red"))+
  geom_line(aes(x = ovarian_cat, y = predicted_prev, colour = "red"))+
  theme_bw()+
  labs(x = "Ovarian Category", y = "Prevalence", 
       title = "Predicted versus observed prevalence of T. congolense infection using RSS")+
  theme(plot.title = element_text(hjust = 0.5, size = 10))+
  scale_color_identity(name = "Data",
                       breaks = c("blue", "red"),
                       labels = c("Observed", "Predicted"),
                       guide = "legend")

# VIVAX

opt.viv.logRSS <- optim(par = params1,
                        fn = function(x) 10*log(RSS.v(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.viv.logRSS <- exp(opt.viv.logRSS$par)
opt.viv.logRSS1 <- opt.viv.logRSS
opt.viv.logRSS <- optim(par = params2,
                        fn = function(x) 10*log(RSS.v(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.viv.logRSS <- exp(opt.viv.logRSS$par)
opt.viv.logRSS2 <- opt.viv.logRSS
opt.viv.logRSS <- optim(par = params3,
                        fn = function(x) 10*log(RSS.v(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.viv.logRSS <- exp(opt.viv.logRSS$par)
opt.viv.logRSS3 <- opt.viv.logRSS
opt.viv.logRSS <- optim(par = params4,
                        fn = function(x) 10*log(RSS.v(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.viv.logRSS <- exp(opt.viv.logRSS$par)
opt.viv.logRSS4 <- opt.viv.logRSS
opt.viv.logRSS1
opt.viv.logRSS2
opt.viv.logRSS3
opt.viv.logRSS4

# Stable and reasonable, plot to see differece between predicted and observed
# plot with confidence intervals 

ci.viv.rss <- binom.confint(x = t.viv.inf, n = t.viv.n, conf.level = 0.95, methods = "exact")

viv.obs.pred.rss <- cbind(viv.data, pooled.prev(loglambda = log(opt.viv.logRSS1[1])
                                                , logtau = log(opt.viv.logRSS1[2]))
                          , lower = ci.viv.rss$lower, upper=ci.viv.rss$upper)
viv.obs.pred.rss$observed_prev[1] <- viv.obs.pred.rss$lower[1] <- viv.obs.pred.rss$upper[1]  <- 0



ggplot(viv.obs.pred.rss, aes(x = ovarian_cat, y = observed_prev, colour = "blue"))+
  geom_ribbon(aes(ymin = lower, ymax=upper), alpha = 0.3, colour = NA)+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = seq(0,7,1), labels = axis.ov.cat)+
  geom_point(aes(x = ovarian_cat, y = predicted_prev, colour = "red"))+
  geom_line(aes(x = ovarian_cat, y = predicted_prev, colour = "red"))+
  theme_bw()+
  labs(x = "Ovarian Category", y = "Prevalence", 
       title = "Predicted versus observed prevalence of T. vivax infection using RSS")+
  theme(plot.title = element_text(hjust = 0.5, size = 10))+
  scale_color_identity(name = "Data",
                       breaks = c("blue", "red"),
                       labels = c("Observed", "Predicted"),
                       guide = "legend")


########################################################################################

# FUNCTION: negative log likelihood for MLE

# congolense
nllikelihood.con <- function(a=0:47, loglambda, logtau) {
  pooled = pooled.prev(a, loglambda, logtau)
  probability = (pooled$predicted_prev)
  nlls <- -dbinom(t.con.inf, t.con.n, probability, log = T)
  return(sum(nlls))
}

nllikelihood.con(loglambda = log(0.1), logtau = log(3.0))

# vivax
nllikelihood.viv <- function(a=0:47, loglambda, logtau) {
  pooled = pooled.prev(a, loglambda, logtau)
  probability = (pooled$predicted_prev)
  nlls <- -dbinom(t.viv.inf, t.viv.n, probability, log = T)
  return(sum(nlls))
}

nllikelihood.viv(loglambda = log(0.1), logtau = log(3.0))

########################################################################################

# Different values or parameters lambda and tau using MLE

# CONGOLENSE

# Testing using contour plots

# MLE for different values of loglambda
loglambda.mle <-  seq(log(0.00001), log(0.01), length.out = 100)
plot(exp(loglambda.mle)
     , log(sapply(X=loglambda.mle, FUN = nllikelihood.con, a = 0:47
                  , logtau = log(3.2))), type = 'l'
     , main = "MLE function values of lambda for T. congolense"
     , ylab = "nllikelihood(lambda)"
     , xlab = "Different values of lambda")
abline(v = out.c.lambda, col = 'red')
# looks fine, minimum is at the right place

# Looking at the contour plot in an "inelegant" way
con.n <- 20
loglambda.mle <-  seq(log(0.001), log(0.003), length.out = con.n)
loglambda.n.mle <- length(loglambda.mle)
logtau.mle <-  seq(log(2), log(4), length.out = con.n)
logtau.n.mle <- length(logtau.mle)
mle.c <- matrix(NA, nrow = loglambda.n.mle, ncol = logtau.n.mle)

# lambda by row, tau by column
for (ii in 1:loglambda.n.mle){
  ii.loglambda <- loglambda.mle[ii]
  for (jj in 1:logtau.n.mle){
    jj.logtau <- logtau.mle[jj]
    mle.c[ii,jj] <- nllikelihood.con(loglambda = ii.loglambda, logtau = jj.logtau)
  }
}

# look at contour plot
contour(x = exp(loglambda.mle),
        y = exp(logtau.mle),
        log(mle.c), nlevels = 40, col = terrain.colors(50),
        plot.title = title(main = "Contour plot for MLE values of lambda and tau for T. congolense", 
                           xlab='lambda',  ylab='tau'))
points(out.c.lambda, out.c.tau, pch = 20)

# VIVAX

# Testing using contour plots

# MLE for different values of loglambda
loglambda.vec <-  seq(log(0.00001), log(0.01), length.out = 100)
plot(exp(loglambda.vec)
     , log(sapply(X=loglambda.vec, FUN = nllikelihood.viv, a = 0:47
                  , logtau = log(1.7))), type = 'l'
     , main = "MLE function values of lambda for T. vivax"
     , ylab = "RSS.v(lambda)"
     , xlab = "Different values of lambda")

abline(v = out.v.lambda, col = 'red')
# minimum is at the right place

# Looking at the contour plot in an "inelegant" way
viv.n <- 20
loglambda.vec <-  seq(log(0.001), log(0.003), length.out = viv.n)
loglambda.n <- length(loglambda.vec)
logtau.vec <-  seq(log(1), log(3), length.out = viv.n)
logtau.n <- length(logtau.vec)
mle.v <- matrix(NA, nrow = loglambda.n, ncol = logtau.n)

# lambda by row, tau by column
for (ii in 1:loglambda.n){
  ii.loglambda <- loglambda.vec[ii]
  for (jj in 1:logtau.n){
    jj.logtau <- logtau.vec[jj]
    mle.v[ii,jj] <- nllikelihood.viv(loglambda = ii.loglambda, logtau = jj.logtau)
  }
}

# look at contour plot
contour(x = exp(loglambda.vec),
        y = exp(logtau.vec),
        log(mle.v), nlevels = 40, col = terrain.colors(50),
        plot.title = title(main = "Contour plot for MLE values of lambda and tau for T. vivax", 
                           xlab='lambda',  ylab='tau'))
points(out.v.lambda, out.v.tau, pch = 20)

########################################################################################

# Maximum Likelihood estimation

# CONGOLENSE

trace <- 3

params1 <- c(loglambda = log(0.001), logtau = log(4.8))      # initial parameter guesses
params2 <- c(loglambda = log(0.002), logtau = log(3.2))      # initial parameter guesses
params3 <- c(loglambda = log(0.001), logtau = log(4))      # initial parameter guesses
params4 <- c(loglambda = log(0.002), logtau = log(2))      # initial parameter guesses

# Getting error for the first run of the optimizer. 

opt.con.logMLE <- optim(par = params1,
                        fn = function(x) 10*log(nllikelihood.con(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.con.logMLE <- exp(opt.con.logMLE$par)

opt.con.logMLE1 <- opt.con.logMLE
opt.con.logMLE <- optim(par = params2,
                        fn = function(x) 10*log(nllikelihood.con(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.con.logMLE <- exp(opt.con.logMLE$par)
opt.con.logMLE2 <- opt.con.logMLE
opt.con.logMLE <- optim(par = params3,
                        fn = function(x) 10*log(nllikelihood.con(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.con.logMLE <- exp(opt.con.logMLE$par)
opt.con.logMLE3 <- opt.con.logMLE
opt.con.logMLE <- optim(par = params4,
                        fn = function(x) 10*log(nllikelihood.con(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.con.logMLE <- exp(opt.con.logMLE$par)
opt.con.logMLE4 <- opt.con.logMLE
opt.con.logMLE1
opt.con.logMLE2
opt.con.logMLE3
opt.con.logMLE4

# Plotting prevalences with confidence intervals 

ci.con.mle <- binom.confint(x = t.con.inf, n = t.con.n, conf.level = 0.95, methods = "exact")

con.obs.pred.mle <- cbind(con.data, pooled.prev(loglambda = log(opt.con.logMLE1[1])
                                                , logtau = log(opt.con.logMLE1[2]))
                          , lower = ci.con.mle$lower, upper=ci.con.mle$upper)
con.obs.pred.mle$observed_prev[1] <- con.obs.pred.mle$lower[1] <- con.obs.pred.mle$upper[1]  <- 0

ggplot(con.obs.pred.mle, aes(x = ovarian_cat, y = observed_prev, colour = "blue"))+
  geom_ribbon(aes(ymin = lower, ymax=upper), alpha = 0.3, colour = NA)+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = seq(0,7,1), labels = axis.ov.cat)+
  geom_point(aes(x = ovarian_cat, y = predicted_prev, colour = "red"))+
  geom_line(aes(x = ovarian_cat, y = predicted_prev, colour = "red"))+
  theme_bw()+
  labs(x = "Ovarian Category", y = "Prevalence", 
       title = "Predicted versus observed prevalence of T. congolense infection using MLE")+
  theme(plot.title = element_text(hjust = 0.5, size = 10))+
  scale_color_identity(name = "Data",
                       breaks = c("blue", "red"),
                       labels = c("Observed", "Predicted"),
                       guide = "legend")

# VIVAX

opt.viv.logMLE <- optim(par = params1,
                        fn = function(x) 10*log(nllikelihood.viv(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.viv.logMLE <- exp(opt.viv.logMLE$par)
opt.viv.logMLE1 <- opt.viv.logMLE
opt.viv.logMLE <- optim(par = params2,
                        fn = function(x) 10*log(nllikelihood.viv(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.viv.logMLE <- exp(opt.viv.logMLE$par)
opt.viv.logMLE2 <- opt.viv.logMLE
opt.viv.logMLE <- optim(par = params3,
                        fn = function(x) 10*log(nllikelihood.viv(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.viv.logMLE <- exp(opt.viv.logMLE$par)
opt.viv.logMLE3 <- opt.viv.logMLE
opt.viv.logMLE <- optim(par = params4,
                        fn = function(x) 10*log(nllikelihood.viv(a=0:47, x[1], x[2])),
                        control = list(trace = trace, maxit = 800))
opt.viv.logMLE <- exp(opt.viv.logMLE$par)
opt.viv.logMLE4 <- opt.viv.logMLE
opt.viv.logMLE1
opt.viv.logMLE2
opt.viv.logMLE3
opt.viv.logMLE4

# Stable and reasonable, plot to see differece between predicted and observed

ci.viv.mle <- binom.confint(x = t.viv.inf, n = t.viv.n, conf.level = 0.95, methods = "exact")

viv.obs.pred.mle <- cbind(viv.data, pooled.prev(loglambda = log(opt.viv.logMLE1[1])
                                                , logtau = log(opt.viv.logMLE1[2]))
                          , lower = ci.viv.mle$lower, upper=ci.viv.mle$upper)
viv.obs.pred.mle$observed_prev[1] <- viv.obs.pred.mle$lower[1] <- viv.obs.pred.mle$upper[1]  <- 0

ggplot(viv.obs.pred.mle, aes(x = ovarian_cat, y = observed_prev, colour = "blue"))+
  geom_ribbon(aes(ymin = lower, ymax=upper), alpha = 0.3, colour = NA)+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = seq(0,7,1), labels = axis.ov.cat)+
  geom_point(aes(x = ovarian_cat, y = predicted_prev, colour = "red"))+
  geom_line(aes(x = ovarian_cat, y = predicted_prev, colour = "red"))+
  theme_bw()+
  labs(x = "Ovarian Category", y = "Prevalence", 
       title = "Predicted versus observed prevalence of T. vivax infection using MLE")+
  theme(plot.title = element_text(hjust = 0.5, size = 10))+
  scale_color_identity(name = "Data",
                       breaks = c("blue", "red"),
                       labels = c("Observed", "Predicted"),
                       guide = "legend")