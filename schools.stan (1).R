set.seed(123)
library("rstan") # observe startup messages

# see: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE) 

# we used a text editor to copy and paste the model
# see schools.stan.txt for the model

# We let theta be a transformation of mu, eta, and tau instead of
# declaring theta in the parameters block
# this allows the sampler to run more efficiently
# see: https://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html

schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

#find /usr/local/lib/R/site-library/ /usr/lib/R/library/ /usr/lib/R/site-library/ ~/.local/lib/ find -iname '*rds' -a -size 0 

fit1 <- stan(file = 'schools.stan', 
            data = schools_dat,          # named list of data
            chains = 4,                  # number of Markov chains
            warmup = 1000,               # number of warmup iterations per chain
            iter = 2000,                 # total number of iterations per chain
            cores = 2,                   # number of cores (using 2 just for the vignet)
            refresh = 1000)              # show progress every 'refresh' iterations

print(fit1, pars=c("theta", "mu", "tau", "lp__"), probs=c(.1,.5,.9))

# print out should appear

traceplot(fit1, pars = c("mu", "tau"), inc_warmup = TRUE, nrow = 2)

print(fit1, pars = c("mu", "tau"))

plot(fit1)

library(shinystan)

sampler_params <- get_sampler_params(fit1, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)

pairs(fit1, pars=c("mu", "tau", "lp__"), las = 1)

# lp__ means log of unnormalized joint posterior density

#maybe acceptance rate?
sum(do.call(rbind, sampler_params)[,5] == 0) / length( do.call(rbind, sampler_params)[,5] )

mean(do.call(rbind, sampler_params)[,1])
