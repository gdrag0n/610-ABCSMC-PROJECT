# Load required libraries
library(tidyverse)

# Initialize variables: tolerances and population sizes
e <- c(3000, 1400, 600, 140, 40)  # Tolerance schedule
T <- length(e)  # Number of populations
N <- 1000       # Number of particles

# Priors for models and parameters
models <- 1:2
parameter_prior <- list(
  model1 = function() runif(1, 0, 100),  # Uniform prior for Model 1
  model2 = function() runif(1, 0, 100)  # Uniform prior for Model 2
)

# Transition probabilities
transition_probabilities <- matrix(c(
  0.8, 0.2,  # Model 1
  0.2, 0.8   # Model 2
), nrow = 2, byrow = TRUE)

# Observed data
observed_data <- data.frame(
  time = c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045,
           0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095),
  Y = c(3, 8, 15, 18, 23, 27, 32, 35, 37, 37, 38, 39, 40, 41, 41, 42, 42, 42, 42, 42)
)

# Distance function: Mean-squared error
distance <- function(D0, D_s) {
  min_len <- min(length(D0), length(D_s))
  return(mean((D0[1:min_len] - D_s[1:min_len])^2))
}

# Model refinement
refined <- function(m, models, transition_probabilities) {
  return(sample(models, 1, prob = transition_probabilities[m, ]))
}

# Parameter perturbation kernel
perturb_parameters <- function(theta, sd) {
  return(rnorm(1, mean = theta, sd = sd))
}

# Function to normalize weights safely
normalize_weights <- function(weights) {
  if (all(!is.na(weights)) && sum(weights) > 0) {
    return(weights / sum(weights))  # Normalize weights
  } else {
    warning("Weights contain NA or sum to zero. Resetting to uniform weights.")
    return(rep(1 / length(weights), length(weights)))  # Reset to uniform weights
  }
}

# Function to calculate ESS
calculate_ESS <- function(weights) {
  if (all(!is.na(weights)) && sum(weights) > 0) {
    return(1 / sum(weights^2))
  } else {
    warning("Invalid weights detected. Setting ESS to 0.")
    return(0)
  }
}


# Gillespie Model 1: X + Y -> 2Y
gillespie_model1 <- function(X0, Y0, k1, t_max) {
  t <- 0
  X <- X0
  Y <- Y0
  times <- c(t)
  X_values <- c(X)
  Y_values <- c(Y)
  
  while (t < t_max) {
    a1 <- k1 * X * Y
    if (a1 <= 0) break  # Stop if propensity is zero
    
    t <- t + rexp(1, rate = a1)  # Time until next reaction
    X <- X - 1
    Y <- Y + 1
    
    # Save values
    times <- c(times, t)
    X_values <- c(X_values, X)
    Y_values <- c(Y_values, Y)
  }
  
  return(data.frame(time = times, X_values = X_values, Y_values = Y_values))
}

# Gillespie Model 2: X -> Y
gillespie_model2 <- function(X0, Y0, k2, t_max) {
  t <- 0
  X <- X0
  Y <- Y0
  times <- c(t)
  X_values <- c(X)
  Y_values <- c(Y)
  
  while (t < t_max) {
    a2 <- k2 * X
    if (a2 <= 0) break  # Stop if propensity is zero
    
    t <- t + rexp(1, rate = a2)  # Time until next reaction
    X <- X - 1
    Y <- Y + 1
    
    # Save values
    times <- c(times, t)
    X_values <- c(X_values, X)
    Y_values <- c(Y_values, Y)
  }
  
  return(data.frame(time = times, X_values = X_values, Y_values = Y_values))
}


# Simulate data function
simulate_data <- function(model, parameter, observed_data) {
  X0 <- 40
  Y0 <- 3
  t_max <- max(observed_data$time)
  
  if (model == 1) {
    return(gillespie_model1(X0, Y0, k1 = parameter, t_max = t_max))
  } else if (model == 2) {
    return(gillespie_model2(X0, Y0, k2 = parameter, t_max = t_max))
  } else {
    stop("Invalid model selected")
  }
}



#Garrett Plot
# Initial conditions and parameters
X0 <- 40
Y0 <- 3
t_max <- 0.1
k1 <- 2.1
k2 <- 30

# Data from Table 1 (Supplementary Material)
table1_time <- c(0.000, 0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045,
                 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095)
table1_Y <- c(3, 8, 15, 18, 23, 27, 32, 35, 37, 37, 38, 39, 40, 41, 41, 42, 42, 42, 42, 42)

# Run Gillespie simulations
results1 <- gillespie_model1(X0, Y0, k1, t_max)
results2 <- gillespie_model2(X0, Y0, k2, t_max)

#plotting
ggplot() +
  geom_line(data = results1, aes(x = time, y = X_values), color = "red", linetype = "dashed", size = 1) +
  geom_line(data = results1, aes(x = time, y = Y_values), color = "blue", linetype = "dashed", size = 1) +
  geom_line(data = results2, aes(x = time, y = X_values), color = "red", linewidth = 1) +
  geom_line(data = results2, aes(x = time, y = Y_values), color = "blue", linewidth = 1) +
  geom_point(aes(x = table1_time, y = table1_Y), color = "blue", size = 2) +
  labs(title = "Stochastic trajectories of species X (red) and Y (blue)",
       x = "Time", y = "Number of molecules") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))




# Weight computation
compute_S <- function(model_star, parameter_star, particles_prev, weights_prev, models, sd_param) {
  prior_prob <- 1 / length(models)
  likelihood_prior <- prior_prob * dnorm(parameter_star, mean = 0, sd = 100)
  
  sum_weights <- sum(sapply(1:length(particles_prev), function(j) {
    trans_prob <- transition_probabilities[particles_prev[[j]]$model, model_star]
    perturb_density <- dnorm(parameter_star, mean = particles_prev[[j]]$parameter, sd = sd_param)
    weights_prev[j] * trans_prob * perturb_density
  }))
  
  return(likelihood_prior / sum_weights)
}

# Initialize particles and weights
particles <- vector("list", T)
weights <- vector("list", T)
particles[[1]] <- lapply(1:N, function(i) {
  model <- sample(models, 1)
  parameter <- parameter_prior[[paste0("model", model)]]()
  list(model = model, parameter = parameter)
})
weights[[1]] <- rep(1 / N, N)

# ABC SMC Algorithm
for (t in 2:T) {
  particles[[t]] <- vector("list", N)
  weights[[t]] <- numeric(N)
  
  for (i in 1:N) {
    repeat {
      idx <- sample(1:N, 1, prob = weights[[t - 1]])
      model_star <- refined(particles[[t - 1]][[idx]]$model, models, transition_probabilities)
      parameter_star <- perturb_parameters(particles[[t - 1]][[idx]]$parameter, 5)
      
      D_star <- simulate_data(model_star, parameter_star, observed_data)
      if (distance(observed_data$Y, D_star$Y_values) <= e[t]) {
        particles[[t]][[i]] <- list(model = model_star, parameter = parameter_star)
        break
      }
    }
    
    weights[[t]][i] <- compute_S(model_star, parameter_star, particles[[t - 1]], weights[[t - 1]], models, 5)
  }
  
  weights[[t]] <- normalize_weights(weights[[t]])
  ESS <- calculate_ESS(weights[[t]])
  
  if (ESS < N / 2) {
    indices <- sample(1:N, N, replace = TRUE, prob = weights[[t]])
    particles[[t]] <- particles[[t]][indices]
    weights[[t]] <- rep(1 / N, N)
  }
}


# Posterior Probabilities
posterior <- sapply(models, function(m) {
  sum(sapply(1:N, function(i) if (particles[[T]][[i]]$model == m) weights[[T]][i] else 0))
})
posterior <- posterior / sum(posterior)




# Results
post_df <- data.frame(probability = posterior, 
                      model_type = c("Model 1", "Model 2"))
ggplot(data = post_df, aes(x=model_type, y=probability)) + 
  geom_bar(color = "black", fill = "tan", stat = "identity") +
  labs(title = "Figure 1b: Posterior Probabilities",
       x = "Model Type",
       y = "Probability") +
  theme_bw()
ggsave("fig-1b.PNG")
             
print(posterior)
plot(1:T, e, type = "b", main = "Tolerance Schedule", xlab = "Population", ylab = "Epsilon")
hist(unlist(lapply(particles[[T]], function(p) p$parameter)), breaks = 20, main = "Posterior Distribution", xlab = "Parameter Values")
