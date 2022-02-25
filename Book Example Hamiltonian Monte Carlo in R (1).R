#programming Hamiltonian Monte Carlo in R
#BDA page 603

library(rstan)
library(arm)

schools = data.frame(                     # data for schools example
  school = LETTERS[1:8],
  estimate = c(28,  8, -3,  7, -1,  1, 18, 12),
  sd = c(15, 10, 16, 11,  9, 11, 10, 18)
)

J <- nrow(schools)
y <- schools$estimate
sigma <- schools$sd

# Log Posterior Density below:

log_p_th <- function(th, y, sigma){
  J <- length(th) - 2                # number of schools
  theta <- th[1:J]                   # theta from "th" vector (J total)
  mu <- th[J+1]                      # mu from "th" vector
  tau <- th[J+2]                     # tau from "th" vector

# Tau is restricted to be positive. Set condition accordingly.   
  
  if (is.nan(tau) | tau<=0)          # tau restricted to be positive
    return(-Inf)                     # posterior density = zero if tau<=0

# We define our Priors for the final Posterior Density 
  
   else{
    log_hyperprior <- 1              
    log_prior <- sum(dnorm(theta, mu, tau, log=TRUE))        # p(theta | mu, tau)
    log_likelihood <- sum(dnorm(y, theta, sigma, log=TRUE))  # p(y | theta, sig)
    return(log_hyperprior + log_prior + log_likelihood)      # log(posterior)
  }
}

# Gradient of the Log Posterior Density below

gradient_th <- function(th, y, sigma){ #th is all params for model
  J <- length(th) - 2
  theta <- th[1:J]
  mu <- th[J+1]
  tau <- th[J+2]
  
# We let gradient = 0 when tau<=0
  if (tau<=0)                        # if tau <=0, gradient=0 for all parameters
    return(c(0,0,0))
  
# We define the gradients wrt their parameters  
  else {
    d_theta <- - (theta-y)/sigma^2 - (theta-mu)/tau^2    # gradient wrt theta
    d_mu <- -sum(mu-theta)/tau^2                         # gradient wrt mu
    d_tau <- -J/tau + sum((mu-theta)^2)/tau^3            # gradient wrt tau
    return(c(d_theta, d_mu, d_tau))                    # returns each gradient
  }
}

# For a more stable gradient we can use the function below:

#
#gradient_th_numerical <- function(th, y, sigma){
#  d <- length(th)                    # d=10 in this example
#  e <- .0001                         # adds small positive amount to stabilize
#  diff <- rep(NA, d)
#  for (k in 1:d){
#    th_hi <- th
#    th_lo <- th
#    th_hi[k] <- th[k] + e
#    th_lo[k] <- th[k] - e 

# Approximate gradient with first differences
#    diff[k]<-(log_p_th(th_hi,y,sigma)-log_p_th(th_lo,y,sigma))/(2*e)
#  }
#  return(diff)
#}


# Single HMC iteration below:

hmc_iteration <- function(th, y, sigma, epsilon, L, M) {

# Arguments for Function:  
# Epsilon is step size
# L is number of leapfrog steps per iteration
# M is diagonal mass matrix
  
  # we will need this for leapfrog step:
  M_inv <- 1/M
  
  # d = 10 in this example:
  d <- length(th)
  
# Initialize all parameters:
  
  phi <- rnorm(d, 0, sqrt(M))
  
  # initial value for all parameters
  th_old <- th
  
# Define the log posterior density of theta_old and phi_old :
  log_p_old <- log_p_th(th,y,sigma) - 0.5*sum(M_inv*phi^2)
  
# Step 2a. complete half step of phi:
  phi <- phi + 0.5*epsilon*gradient_th(th, y, sigma)
  
# Step 2b. complete full step of theta:
  for (l in 1:L){
   
    th <- th + epsilon*M_inv*phi
    
# Step 2c. complete full step of phi:
    # In first/last Leapfrog step, the phi must be found before and after theta
    phi <- phi + (if (l==L) 0.5 else 1)*epsilon*gradient_th(th,y,sigma)
  }
  
  phi <- -phi #preserves property of Hamiltonian
  
# Define the log posterior density of theta* and phi*:
  
  log_p_star <- log_p_th(th,y,sigma) - 0.5*sum(M_inv*phi^2)
  
# Accept / reject step via the ratio of log posterior densities:
  r <- exp(log_p_star - log_p_old)
  
  # control for instability
  if (is.nan(r)) r <- 0
  
  # probability of jumping to new value:
  p_jump <- min(r,1)
  
  # defining new theta value:
  th_new <- if (runif(1) < p_jump) th else th_old
  
  # HMC iteration return:
  return(list(th=th_new, p_jump=p_jump))
}


# HMC algorithm with multiple iterations:

hmc_run <- function(starting_values, iter, epsilon_0, L_0, M) {
  # epsilon_0 is the baseline step size
  # L_0 is the number of steps
  
# Starting_values is a  (m x d) matrix
# m is the number of chains
# d is the number of parameters in the vector, 10 in this example
  chains <- nrow(starting_values)
  d <- ncol(starting_values)
  
  # empty arrays to store results:
  sims <- array(NA, c(iter, chains, d),
                dimnames=list(NULL, NULL, colnames(starting_values)))
  
# We expect about half of our iterations to be part of the warm-up period
  warmup <- 0.5*iter
  p_jump <- array(NA, c(iter, chains))
  
  for (j in 1:chains){
    th <- starting_values[j,]
    for (t in 1:iter){
      epsilon <- runif(1, 0, 2*epsilon_0)     # randomly drawn at each iteration
      L <- ceiling(2*L_0*runif(1))            # randomly drawn to explore joint dist.
      
      # this completes an HMC iteration and creates temp parameter vector:
      temp <- hmc_iteration(th, y, sigma, epsilon, L, M)
      
      # we identify the p_jump var in temp:
      p_jump[t,j] <- temp$p_jump
      
      # we identify the parameters in temp, called "sims"
      # sims is a 3 dimensional array defined earlier
      sims[t,j,] <- temp$th
      
      # we use this th for next iteration
      th <- temp$th
    }
  }
  
  # monitors convergence
  # excludes warmup iterations when computing the summaries
  monitor(sims, warmup)
  
  # average acceptance, after the warm-up period, below:
  
  cat("Avg acceptance probs:",
      fround(colMeans(p_jump[(warmup+1):iter,]),2),"\n")
  
# Track and print the parameter simulations and the average acceptance rate:
  return(list(sims=sims, p_jump=p_jump))
}


# we define a vector with the names of the parameters:
parameter_names <- c(paste("theta[",1:8,"]",sep=""), "mu", "tau")

# define d. it's 10 in this example
d <- length(parameter_names)

# set the number of chains. book said to do 4:
chains <- 4

# We define the diagonal mass matrix to be proportional to 
# the inverse variance matrix of the posterior distribution
# Textbook recommends 15 for this example
mass_vector <- rep(1/15^2, d)

# initialize array for starting points: 
starts <- array(NA,c(chains,d),dimnames=list(NULL,parameter_names))

# Then we assign all 10 starting points from Normal(0,15^2)
for (j in 1:chains){
  starts[j,] <- rnorm(d,0,15)
  starts[j,10] <- runif(1,0,15) # starting point for tau must be positive
}

# FIRST RUN:
# we run only for 20 iterations to make sure our program doesn't crash
M1 <- hmc_run (starting_values=starts, iter=20,
               epsilon_0=.1, L_0=10, M=mass_vector)

# SECOND RUN:
# increase iterations to 100 after first run is successful
M2 <- hmc_run(starting_values=starts, iter=100,
              epsilon_0=.1, L_0=10, M=mass_vector)

# THIRD RUN:
# Decrease step size, epsilon, to improve acceptance rate
# Increase base number of steps also to improve acceptance rate
M3 <- hmc_run(starting_values=starts, iter=100,
              epsilon_0=.05, L_0=20, M=mass_vector)

# run for thousands of iterations until we obtain stable inferences

# So much work! Let's see Stan....




