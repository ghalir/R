#--------------------------------------------------
# © Aliakbar Rasekhi, 2019
# Codes in R for the paper
# Incomplete gamma distribution: a new two parameter lifetime distribution with survival regression model
# A Rasekhi, M Rasekhi, G Hamedani
# UPB Scientific Bulletin, Series A: Applied Mathematics and Physics 81(2), 177-184
#--------------------------------------------------
rm(list = ls())
library(expint)
library(flexsurv)
# density
ding <- function(x, c=0, b=1, log = FALSE) {
  a = exp(c)
  if (log) {
    -log((c + exp(a)*expint(a))) - log(b) + log(log(x/b + a)) - x/b
  } else {
    1/(c + exp(a)*expint(a))/b*log(x/b + a)*exp(-x/b)
  }
}
ding <- Vectorize(ding)
# distribution function
ping <- function(q, c=0, b=1, lower.tail = TRUE, log.p = FALSE) {
  a = exp(c)
  if (log.p) {
  if (lower.tail) { 
      log(1 - 1/(c + exp(a)*expint(a))*(exp(a)*expint(q/b + a) + log(q/b + a)*exp(-q/b)))
    } else {
      log(1/(c + exp(a)*expint(a))*(exp(a)*expint(q/b + a) + log(q/b + a)*exp(-q/b)))
    }
  } else {
  if (lower.tail) { 
    1 - 1/(c + exp(a)*expint(a))*(exp(a)*expint(q/b + a) + log(q/b + a)*exp(-q/b))
    } else {
    1/(c + exp(a)*expint(a))*(exp(a)*expint(q/b + a) + log(q/b + a)*exp(-q/b))  
    }
  }
}
ping <- Vectorize(ping)
# quantile function
qing <- function(p, c=0, b=1) {   # lower.tail = TRUE, log.p = FALSE
  eq <- function(q) { ping(q,c,b) - p }
  sol <- uniroot(eq, lower = 0, upper = b*100)
  return(sol$root)
}
#
# InG distribution for use in flexsurvreg()
InG <- list(name = "ing", pars = c("c","b"), location = "b", 
            transforms = c(log, log), inv.transforms = c(exp, exp),
            inits = function(t){ c(1, median(t)) })
#--------------------------------------------------
# Fit to a dataset; for help run: ?bc
data(bc)
head(bc)
#--------------------------------------------------
fitw <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group,
                    data = bc, dist = "weibull")
#fitw
fiti <- flexsurvreg(formula = Surv(recyrs, censrec) ~ group,
                    data = bc, dist = InG)
#fiti
fiti$opt$convergence
#--------------------------------------------------
# weibull table
tablew = fitw$res.t
phw = coef(fitw)
sew = sqrt(diag(vcov(fitw)))
pvw = 2*(1 - pnorm(abs(phw/sew)))
cbind(tablew,p.value = pvw)
fitw$AIC
# InG table
tablei = fiti$res.t
phi = coef(fiti)
sei = sqrt(diag(vcov(fiti)))
pvi = 2*(1 - pnorm(abs(phi/sei)))
cbind(tablei,p.value = pvi)
fiti$AIC
#--------------------------------------------------
plot(fitw, xlab = "Time (years)", ylab = "Survival",
     xlim = c(0,8), lwd.ci = 1, lty.ci = 1)
lines(fiti, col = 3, lwd.ci = 1, lty.ci = 1)
legend("bottomleft",c("Kaplan-Meier","Weibull","ING"),
       col = 1:3, lty = 1, lwd = 2, cex = .9)
text(7.5,.7,"Good", cex = .8)
text(7.2,.4,"Medium", cex = .8)
text(7.2,.1,"Poor", cex = .8)
title("Figure 2. BC data")
#--------------------------------------------------
