## Higher than log 2


#------------------------------------------------------------------------------------
## Clear the environment in R ##
#------------------------------------------------------------------------------------
rm(list = ls())


#------------------------------------------------------------------------------------
## Install the required packages into ##
## R if it's not already installed.   ##
#------------------------------------------------------------------------------------
if(!require(latex2exp)) install.packages("latex2exp")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(parallel)) install.packages("parallel")
if(!require(doParallel)) install.packages("doParallel")


#------------------------------------------------------------------------------------
## Load the packages into the R environment ##
#------------------------------------------------------------------------------------
library(latex2exp)
library(ggplot2)
library(parallel)
library(doParallel)


#------------------------------------------------------------------------------------
seed <- 123

H <- function(x, beta, p){
  return(-(beta*x^p - 1/2 *((1+x)*log(1+x)+(1-x)*log(1-x))))
}

second_der_H <- function(x, beta, p){
  return(beta*p*(p-1)*x^(p-2) - 1/(1-x^2))
}

m_star <- function(beta, p){
  if(beta > 0){
    return(optim(par=0.50, fn = H, beta = beta, p = p, method="L-BFGS-B", lower = 0.8,
                 upper = 1-exp(-20))$par)
  }
  else{
    return("Choose beta > 0")
  }
}

normal_var <- function(m_star, beta, p){
  return(-second_der_H(m_star, beta, p)/(p^2 * m_star^(2*p-2)))
}

beta_star <- function(p){
  if(p == 2){
    return(0.5)
  }
  else if(p == 3){
    return(0.672)
  }
  else if(p == 4){
    return(0.6888)
  }
}


#------------------------------------------------------------------------------------
## Theorem 1  ##
#------------------------------------------------------------------------------------
min_N_ML <- function(delta, beta, beta0, p){
  denom1 <- -H(m_star(beta, p), beta0, p)
  denom2 <- -H(m_star(beta0, p), beta0, p)
  return(log(delta)/(denom1 - denom2))
}

# Theorem 1
min_N_MPL <- function(delta, beta, beta0, p){
  if(beta > beta0 && beta0 > beta_star(p)){
    if(p == 2){
      denom1 <- -H(m_star(beta, p), beta0, p)
    }
    else if(p > 2){
      denom1 <- max(-H(m_star(beta, p), beta0, p),0)
    }
    denom2 <- -H(m_star(beta0, p), beta0, p)
    return(log(delta)/(denom1 - denom2))
  }
  else{
    return("Try again.")
  }
}


#------------------------------------------------------------------------------------
## Some Calculations ##
## Calculate the PMF of the CW model at n, beta, p ##
#------------------------------------------------------------------------------------
pmf <- function(n, beta, p){
  m <- seq(-1, 1, by = 2/n)
  boltzmann <- function(x){
    return(choose(n, n*(1+x)/2) * exp(n * (beta * x^p - log(2))))
  }
  y <- sapply(m, boltzmann)
  return(y/sum(y))
}

# Calculate the pth moment of magnetization (for MLE)
expec_magp <- function(n, beta, p){
  mf <- pmf(n, beta, p)
  m <- seq(-1, 1, by = 2/n)
  mp <- m^p
  return((mf %*% mp)[1])
}

# Generate (number = reps) random average magnetization from the CW model at n, beta, p
pmf_rng <- function(reps, n, beta, p){
  mf <- pmf(n, beta, p)
  m <- seq(-1, 1, by = 2/n)
  return(round(sample(m, reps, replace = TRUE, prob = mf), 6))
}

# Calculate the MPLE given the average magnetization
mple <- function(m, p){
  if(m != 0){
    return(m^(1-p)* atanh(m)/p)
  }
  else{
    return(0)
  }
}

# Calculate the MLE given the average magnetization
mle <- function(n, m, p, list_exp, beta_seq){
  loc <- which.min(abs(list_exp - m^p))
  return(beta_seq[loc])
}


#------------------------------------------------------------------------------------
## Find N in Theorem 1 using exact distribution ##
## Find MPLE ##
#------------------------------------------------------------------------------------
find_N_mpl <- function(n, interval, delta, reps, beta, beta0, p, seed = 123){
  set.seed(seed)
  n_seq <- seq(n-interval, n+interval, by = 1)
  l <- integer(length(n_seq))
  mple_p <- function(m){
    return(mple(m, p))
  }
  for(i in n_seq){
    set.seed(seed)
    m <- seq(-1, 1, by = 2/i)
    mple_seq <- sapply(m, mple_p)
    mf <- pmf(i, beta0, p)
    df <- function(x){
      w <- which(mple_seq <= x)
      return(sum(mf[w]))
    }
    get_p_value <- function(mple){
      return(1-df(mple))
    }
    samples <- pmf_rng(reps, i, beta, p)
    mples <- sapply(samples, mple_p)
    all_p_values <- sapply(mples, get_p_value)
    l[(i-(n-interval - 1))] <- mean(all_p_values)
  }
  return(l)
}


# Find MLE
find_N_ml <- function(n, interval, delta, reps, beta, beta0, p, seed = 123){
  set.seed(seed)
  n_seq <- seq(n-interval, n+interval, by = 1)
  l <- integer(length(n_seq))
  beta_seq <- seq(0, 2.5, by = 0.001)
  for(i in n_seq){
    set.seed(seed)
    m <- seq(-1, 1, by = 2/i)
    mf <- pmf(i, beta0, p)
    expec <- function(beta){
      return(expec_magp(i, beta, p))
    }
    list_exp <- sapply(beta_seq, expec)
    mle_p <- function(m){
      return(mle(i, m, p, list_exp, beta_seq))
    }
    mle_seq <- sapply(m, mle_p)
    df <- function(x){
      w <- which(mle_seq <= x)
      return(sum(mf[w]))
    }
    get_p_value <- function(mle){
      return(1-df(mle))
    }
    samples <- pmf_rng(reps, i, beta, p)
    mles <- sapply(samples, mle_p)
    all_p_values <- sapply(mles, get_p_value)
    l[(i-(n-interval - 1))] <- mean(all_p_values)
  }
  return(l)
}


#------------------------------------------------------------------------------------
## Parallel Computing ##
#------------------------------------------------------------------------------------
par_find_N_mpl <- function(n, interval, delta, reps, beta, beta0, p, seed = 123){
  set.seed(seed)
  n_seq <- seq(n-interval, n+interval, by = 1)
  mple_p <- function(m){
    return(mple(m, p))
  }
  f <- function(i){
    set.seed(seed)
    m <- seq(-1, 1, by = 2/i)
    mple_seq <- sapply(m, mple_p)
    mf <- pmf(i, beta0, p)
    df <- function(x){
      w <- which(mple_seq <= x)
      return(sum(mf[w]))
    }
    get_p_value <- function(mple){
      return(1-df(mple))
    }
    samples <- pmf_rng(reps, i, beta, p)
    mples <- sapply(samples, mple_p)
    all_p_values <- sapply(mples, get_p_value)
    return(mean(all_p_values))
  }
  l <- mcmapply(f, n_seq)
  return(l)
}

par_find_N_ml <- function(n, interval, delta, reps, beta, beta0, p, seed = 123){
  set.seed(seed)
  n_seq <- seq(n-interval, n+interval, by = 1)
  beta_seq <- seq(0, 2.5, by = 0.001)
  f <- function(i){
    set.seed(seed)
    m <- seq(-1, 1, by = 2/i)
    mf <- pmf(i, beta0, p)
    expec <- function(beta){
      return(expec_magp(i, beta, p))
    }
    list_exp <- sapply(beta_seq, expec)
    mle_p <- function(m){
      return(mle(i, m, p, list_exp, beta_seq))
    }
    mle_seq <- sapply(m, mle_p)
    df <- function(x){
      w <- which(mle_seq <= x)
      return(sum(mf[w]))
    }
    get_p_value <- function(mle){
      return(1-df(mle))
    }
    samples <- pmf_rng(reps, i, beta, p)
    mles <- sapply(samples, mle_p)
    all_p_values <- sapply(mles, get_p_value)
    return(mean(all_p_values))
  }
  l <- mcmapply(f, n_seq)
  return(l)
}


#------------------------------------------------------------------------------------
# Plot 1
delta <- 0.01
beta <- 0.8
beta0 <- 0.7
p <- 2
interval <- 100
reps <- 10000

min_ML <- min_N_ML(delta, beta, beta0, p)
min_MPL <- min_N_MPL(delta, beta, beta0, p)
n <- round(min_ML)
n_seq <- seq(n-interval, n+interval, by = 1)

set.seed(seed)
p_values <- par_find_N_mpl(n, interval, delta, reps, beta, beta0, p)
red_seq <- n_seq[which(p_values < delta)]
red_p_values <- p_values[which(p_values < delta)]
black_seq <- n_seq[!(n_seq %in% red_seq)]
black_p_values <- p_values[which(p_values >= delta)]

dfblack = data.frame("black_seq"=black_seq,"black_p_vals"=black_p_values)
dfred = data.frame("red_seq"=red_seq,"black_p_vals"=red_p_values)

pdf("mpleplot80_p2_null_above_log2_point70_delta_1_multiple_seed_123.pdf")
ggplot() +
  geom_point(data=dfblack,aes(black_seq, black_p_values),col="black",shape=1) +
  geom_point(data=dfred,aes(red_seq, red_p_values),col="red",shape=2) +
  labs(x = 'Sample Size', y = 'p-values', xlim = c(min(n_seq), max(n_seq)),
       ylim = c(min(p_values), max(p_values)) ) +
  geom_vline(xintercept = max(black_seq), col = 'blue')  
dev.off()

print(paste('min_ML = ', round(min_ML, 3), ', min_MPL = ', round(min_MPL, 3), ', delta = ', delta,
            ', beta0 = ', beta0, ', beta = ', beta, ' p = ', p, ', blue line = ', max(black_seq)))

# MLE
set.seed(seed)
p_values <- par_find_N_ml(n, interval, delta, reps, beta, beta0, p)
red_seq <- n_seq[which(p_values < delta)]
red_p_values <- p_values[which(p_values < delta)]
black_seq <- n_seq[!(n_seq %in% red_seq)]
black_p_values <- p_values[which(p_values >= delta)]

dfblack = data.frame("black_seq"=black_seq,"black_p_vals"=black_p_values)
dfred = data.frame("red_seq"=red_seq,"black_p_vals"=red_p_values)

pdf("mleplot80_p2_null_above_log2_point70_delta_1_multiple_seed_123.pdf")
ggplot() +
  geom_point(data=dfblack,aes(black_seq, black_p_values),col="black",shape=1) +
  geom_point(data=dfred,aes(red_seq, red_p_values),col="red",shape=2) +
  labs(x = 'Sample Size', y = 'p-values', xlim = c(min(n_seq), max(n_seq)),
       ylim = c(min(p_values), max(p_values)) ) +
  geom_vline(xintercept = max(black_seq), col = 'blue') 
dev.off()

print(paste('min_ML = ', round(min_ML, 3), ', min_MPL = ', round(min_MPL, 3), ', delta = ', delta,
            ', beta0 = ', beta0, ', beta = ', beta, ' p = ', p, ', blue line = ', max(black_seq)))


#------------------------------------------------------------------------------------
# Plot 2
delta <- 0.01
beta <- 0.85
beta0 <- 0.7
p <- 2
interval <- 100
reps <- 10000


min_ML <- min_N_ML(delta, beta, beta0, p)
min_MPL <- min_N_MPL(delta, beta, beta0, p)
n <- round(min_ML)
n_seq <- seq(n-interval, n+interval, by = 1)

set.seed(seed)
p_values <- par_find_N_mpl(n, interval, delta, reps, beta, beta0, p)
red_seq <- n_seq[which(p_values < delta)]
red_p_values <- p_values[which(p_values < delta)]
black_seq <- n_seq[!(n_seq %in% red_seq)]
black_p_values <- p_values[which(p_values >= delta)]

dfblack = data.frame("black_seq"=black_seq,"black_p_vals"=black_p_values)
dfred = data.frame("red_seq"=red_seq,"black_p_vals"=red_p_values)

pdf("mpleplot85_p2_null_above_log2_point70_delta_1_multiple_seed_123.pdf")
ggplot() +
  geom_point(data=dfblack,aes(black_seq, black_p_values),col="black",shape=1) +
  geom_point(data=dfred,aes(red_seq, red_p_values),col="red",shape=2) +
  labs(x = 'Sample Size', y = 'p-values', xlim = c(min(n_seq), max(n_seq)),
       ylim = c(min(p_values), max(p_values)) ) +
  geom_vline(xintercept = max(black_seq), col = 'blue')
dev.off()

print(paste('min_ML = ', round(min_ML, 3), ', min_MPL = ', round(min_MPL, 3), ', delta = ', delta,
            ', beta0 = ', beta0, ', beta = ', beta, ' p = ', p, ', blue line = ', max(black_seq)))

# MLE
p_values <- par_find_N_ml(n, interval, delta, reps, beta, beta0, p)
red_seq <- n_seq[which(p_values < delta)]
red_p_values <- p_values[which(p_values < delta)]
black_seq <- n_seq[!(n_seq %in% red_seq)]
black_p_values <- p_values[which(p_values >= delta)]

dfblack = data.frame("black_seq"=black_seq,"black_p_vals"=black_p_values)
dfred = data.frame("red_seq"=red_seq,"black_p_vals"=red_p_values)

pdf("mleplot85_p2_null_above_log2_point70_delta_1_multiple_seed_123.pdf")
ggplot() +
  geom_point(data=dfblack,aes(black_seq, black_p_values),col="black",shape=1) +
  geom_point(data=dfred,aes(red_seq, red_p_values),col="red",shape=2) +
  labs(x = 'Sample Size', y = 'p-values', xlim = c(min(n_seq), max(n_seq)),
       ylim = c(min(p_values), max(p_values)) ) +
  geom_vline(xintercept = max(black_seq), col = 'blue')  
dev.off()

print(paste('min_ML = ', round(min_ML, 3), ', min_MPL = ', round(min_MPL, 3), ', delta = ', delta,
            ', beta0 = ', beta0, ', beta = ', beta, ' p = ', p, ', blue line = ', max(black_seq)))


#------------------------------------------------------------------------------------
# Plot 3
delta <- 0.01
beta <- 0.9
beta0 <- 0.68
p <- 3
interval <- 100
reps <- 10000


min_ML <- min_N_ML(delta, beta, beta0, p)
min_MPL <- min_N_MPL(delta, beta, beta0, p)
n <- round(min_MPL)
n_seq <- seq(n-interval, n+interval, by = 1)

set.seed(seed)
p_values <- par_find_N_mpl(n, interval, delta, reps, beta, beta0, p)
red_seq <- n_seq[which(p_values < delta)]
red_p_values <- p_values[which(p_values < delta)]
black_seq <- n_seq[!(n_seq %in% red_seq)]
black_p_values <- p_values[which(p_values >= delta)]


dfblack = data.frame("black_seq"=black_seq,"black_p_vals"=black_p_values)
dfred = data.frame("red_seq"=red_seq,"black_p_vals"=red_p_values)

pdf("mpleplot90_p3_null_below_log2_point68_delta_1_multiple_seed_123.pdf")
ggplot() +
  geom_point(data=dfblack,aes(black_seq, black_p_values),col="black",shape=1) +
  geom_point(data=dfred,aes(red_seq, red_p_values),col="red",shape=2) +
  labs(x = 'Sample Size', y = 'p-values', xlim = c(min(n_seq), max(n_seq)),
       ylim = c(min(p_values), max(p_values)) ) +
  geom_vline(xintercept = max(black_seq), col = 'blue')  
dev.off()

print(paste('min_ML = ', round(min_ML, 3), ', min_MPL = ', round(min_MPL, 3), ', delta = ', delta,
            ', beta0 = ', beta0, ', beta = ', beta, ' p = ', p, ', blue line = ', max(black_seq)))


# MLE
n <- round(min_ML)
n <- n - interval
n_seq <- seq(n-interval, n+interval, by = 1)
set.seed(seed)
p_values <- par_find_N_ml(n, interval, delta, reps, beta, beta0, p)
red_seq <- n_seq[which(p_values < delta)]
red_p_values <- p_values[which(p_values < delta)]
black_seq <- n_seq[!(n_seq %in% red_seq)]
black_p_values <- p_values[which(p_values >= delta)]


dfblack = data.frame("black_seq"=black_seq,"black_p_vals"=black_p_values)
dfred = data.frame("red_seq"=red_seq,"black_p_vals"=red_p_values)

pdf("mleplot90_p3_null_below_log2_point68_delta_1_multiple_seed_123.pdf")
ggplot() +
  geom_point(data=dfblack,aes(black_seq, black_p_values),col="black",shape=1) +
  geom_point(data=dfred,aes(red_seq, red_p_values),col="red",shape=2) +
  labs(x = 'Sample Size', y = 'p-values', xlim = c(min(n_seq), max(n_seq)),
       ylim = c(min(p_values), max(p_values)) ) +
  geom_vline(xintercept = max(black_seq), col = 'blue')  
dev.off()

print(paste('min_ML = ', round(min_ML, 3), ', min_MPL = ', round(min_MPL, 3), ', delta = ', delta,
            ', beta0 = ', beta0, ', beta = ', beta, ' p = ', p, ', blue line = ', max(black_seq)))


#------------------------------------------------------------------------------------
# Plot 4
delta <- 0.01
beta <- 0.90
beta0 <- 0.7
p <- 2
interval <- 100
reps <- 10000


min_ML <- min_N_ML(delta, beta, beta0, p)
min_MPL <- min_N_MPL(delta, beta, beta0, p)
n <- round(min_ML)
n_seq <- seq(n-interval, n+interval, by = 1)

set.seed(seed)
p_values <- par_find_N_mpl(n, interval, delta, reps, beta, beta0, p)
red_seq <- n_seq[which(p_values < delta)]
red_p_values <- p_values[which(p_values < delta)]
black_seq <- n_seq[!(n_seq %in% red_seq)]
black_p_values <- p_values[which(p_values >= delta)]


dfblack = data.frame("black_seq"=black_seq,"black_p_vals"=black_p_values)
dfred = data.frame("red_seq"=red_seq,"black_p_vals"=red_p_values)

pdf("mpleplot90_p2_null_above_log2_point7_delta_1_multiple_seed_123.pdf")
ggplot() +
  geom_point(data=dfblack,aes(black_seq, black_p_values),col="black",shape=1) +
  geom_point(data=dfred,aes(red_seq, red_p_values),col="red",shape=2) +
  labs(x = 'Sample Size', y = 'p-values', xlim = c(min(n_seq), max(n_seq)),
       ylim = c(min(p_values), max(p_values)) ) +
  geom_vline(xintercept = max(black_seq), col = 'blue')  
dev.off()


print(paste('min_ML = ', round(min_ML, 3), ', min_MPL = ', round(min_MPL, 3), ', delta = ', delta,
            ', beta0 = ', beta0, ', beta = ', beta, ' p = ', p, ', blue line = ', max(black_seq)))

# MLE
set.seed(seed)
p_values <- par_find_N_ml(n, interval, delta, reps, beta, beta0, p)
red_seq <- n_seq[which(p_values < delta)]
red_p_values <- p_values[which(p_values < delta)]
black_seq <- n_seq[!(n_seq %in% red_seq)]
black_p_values <- p_values[which(p_values >= delta)]


dfblack = data.frame("black_seq"=black_seq,"black_p_vals"=black_p_values)
dfred = data.frame("red_seq"=red_seq,"black_p_vals"=red_p_values)

pdf("mleplot90_p2_null_above_log2_point7_delta_1_multiple_seed_123.pdf")
ggplot() +
  geom_point(data=dfblack,aes(black_seq, black_p_values),col="black",shape=1) +
  geom_point(data=dfred,aes(red_seq, red_p_values),col="red",shape=2) +
  labs(x = 'Sample Size', y = 'p-values', xlim = c(min(n_seq), max(n_seq)), 
       ylim = c(min(p_values), max(p_values)) ) +
  geom_vline(xintercept = max(black_seq), col = 'blue') 
dev.off()


print(paste('min_ML = ', round(min_ML, 3), ', min_MPL = ', round(min_MPL, 3), ', delta = ', delta,
            ', beta0 = ', beta0, ', beta = ', beta, ' p = ', p, ', blue line = ', max(black_seq)))


