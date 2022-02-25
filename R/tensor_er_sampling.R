#------------------------------------------------------------------------------------
## Install the required packages into ##
## R if it's not already installed.   ##
#------------------------------------------------------------------------------------
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(parallel)) install.packages("parallel")
if(!require(doParallel)) install.packages("doParallel")
if(!require(doSNOW)) install.packages("doSNOW")
if(!require(doRNG)) install.packages("doRNG")


#------------------------------------------------------------------------------------
## Load the required packages into the R environment ##
#------------------------------------------------------------------------------------
library(ggplot2)
library(parallel)
library(doParallel)
library(doSNOW)
library(doRNG)


#------------------------------------------------------------------------------------
# Number of cores
cores <- detectCores()


#### Energy ####
energy_matrix <- function(x, A, alpha, N){
  return(2 * (alpha * N)^{-1} * (t(x) %*% A %*% x))
}

energy_tensor <- function(x, A, p = 3){
  N <- length(x)
  X <- x %*% t(x)
  X_tensor <- array(NA, dim = rep(N, p))
  for(i in 1:N){
    X_tensor[i, , ] <- X * x[i]
  }
  return(2 * sum(A * X_tensor))
}


#### Generate the random ER tensor ####
# alpha is edge probability
generate_matrix_A <- function(N, alpha){
  A <- matrix(rbinom(N^2, 1, alpha), N, N)
  A[lower.tri(A, diag = TRUE)] <- 0
  return(A)
}

generate_tensor_A <- function(N, alpha, p = 3){
  A_tensor <- array(NA, dim = rep(N, p))
  for(i in 1:N){
    A <- matrix(rbinom(N^2, 1, alpha), N, N)
    A[lower.tri(A, diag = TRUE)] <- 0
    A_tensor[i, , ] <- A
  }
  return(A_tensor)
}


#### Metropolis Hastings ####
# Flip Probability for MH algorithm
flip_probability <- function(energy1, energy2, beta){
  return(exp((-energy1 + energy2) * beta))
}

MH_matrix <- function(iter, A, beta0, alpha, N){
  set.seed(123)
  spins = rep(1, N)
  current_energy <- energy_matrix(spins, A, alpha, N)
  for(i in 1:iter){
    spins_new <- spins
    i <- sample(1:N, 1)
    spins_new[i] <- spins_new[i] * (-1)
    new_energy <- energy_matrix(spins_new, A, alpha, N)
    if(new_energy > current_energy){
      current_energy <- new_energy; spins <- spins_new
    }
    else{
      if(flip_probability(current_energy, new_energy, beta0) > runif(1)){
        current_energy <- new_energy; spins <- spins_new
      }
    }
  }
  return(spins)
}

MH_tensor <- function(iter, A, beta0){
  set.seed(123)
  spins = rep(1, N)
  current_energy <- energy_tensor(spins, A)
  for(i in 1:iter){
    spins_new <- spins
    i <- sample(1:N, 1)
    spins_new[i] <- spins_new[i] * (-1)
    new_energy <- energy_tensor(spins_new, A)
    if(new_energy > current_energy){
      current_energy <- new_energy; spins <- spins_new
    }
    else{
      if(flip_probability(current_energy, new_energy, beta0) > runif(1)){
        current_energy <- new_energy; spins <- spins_new
      }
    }
  }
  return(spins)
}


#### Generate Samples ####
set.seed(123)
beta0 <- 0.7
alpha <- 0.25
iter <- 40000
n_samples <- 10000
N_seq <- seq(150, 250, by = 10)
for(j in N_seq){
  print(j)
  N <- j
  generate_spin <- function(i){
    A <- generate_matrix_A(N, alpha)
    return(list(MH_matrix(iter, A, beta0, alpha, N), A))
  }
  cl <- makeCluster(cores - 2)
  registerDoSNOW(cl)
  samples <- foreach(i = 1:n_samples) %dorng% {
    generate_spin(i)
  }
  saveRDS(samples, paste('beta_0.7_alpha_', as.character(alpha), '_N_', as.character(N), sep = ''))
  stopCluster(cl)
}


set.seed(123)
beta0 <- 0.9
iter <- 40000
n_samples <- 10000
N_seq <- seq(210, 250, by = 10)
alpha <- 0.25
for(j in N_seq){
  print(j)
  N <- j
  generate_spin <- function(i){
    A <- generate_matrix_A(N, alpha)
    return(list(MH_matrix(iter, A, beta0, alpha, N), A))
  }
  cl <- makeCluster(cores - 2)
  registerDoSNOW(cl)
  samples <- foreach(i = 1:n_samples) %dorng% {
    generate_spin(i)
  }
  saveRDS(samples, paste('beta_0.9_alpha_', as.character(alpha), '_N_', as.character(N), sep = ''))
  stopCluster(cl)
}


#### Find Optimal Sample Size from Theorem 1 for MPLE  ####
# Calculate the MPLE given the average magnetization
mple <- function(spins, A, alpha, N, p){
  if(p == 2){
    left <- energy_matrix(spins, A, alpha, N)
    B <- A + t(A)
    vector_m <- 1/(alpha * N) * B %*% spins # Define vector containing m_i(X)
    right <- function(beta){
      return(sum(vector_m * tanh(2 * beta * vector_m)) - left)
    }
    return(uniroot(right, lower = -10, upper = 10, tol = 1e-7)$root)
  }
}


##### Parallel Computing ####
par_find_N_mpl <- function(n_seq, reps, beta, beta0, p, alpha){
  f <- function(i){
    nullsamples <- readRDS(paste('beta_', as.character(beta0), '_alpha_0.25_N_', as.character(i), sep = ''))
    mples <- rep(0, n_samples)
    for(j in 1:n_samples){
      mples[j] <- mple(nullsamples[[j]][[1]], nullsamples[[j]][[2]], alpha, i, p)
    }
    emp_cdf <- ecdf(mples)
    altsamples <- readRDS(paste('beta_', as.character(beta), '_alpha_0.25_N_', as.character(i), sep = ''))
    alt_mples <- rep(0, n_samples)
    for(j in 1:n_samples){
      alt_mples[j] <- mple(altsamples[[j]][[1]], altsamples[[j]][[2]], alpha, i, p)
    }
    df <- sapply(alt_mples, emp_cdf)
    p_values <- 1 - df
    return(mean(p_values))
  }
  l <- mcmapply(f, n_seq)
  return(l)
}


#### Plots ####
# Plot 1
delta <- 0.01
beta <- 0.9
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


