# SEIRS_model.R
# Modelling SEIRS malaria dynamics using synthetic cohort data

# Load libraries
library(deSolve)
library(ggplot2)

# Time steps
times <- seq(0, 365, by = 1)  # 1 year, daily steps

# Initial population and compartments
N <- 1000  # Total population
init <- c(S = 900, E = 20, Ia = 50, Is = 20, R = 10)

# Parameters
params <- c(
  beta = 0.25,     # Transmission rate
  sigma = 1/10,    # Incubation rate (1/latent period)
  gamma_a = 1/30,  # Recovery rate for asymptomatic
  gamma_s = 1/10,  # Recovery rate for symptomatic
  p = 0.6,         # Probability of becoming asymptomatic
  omega = 1/180    # Loss of immunity (R to S)
)

# SEIRS model function
seirs_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dS  <- -beta * S * (Ia + Is) / N + omega * R
    dE  <- beta * S * (Ia + Is) / N - sigma * E
    dIa <- sigma * p * E - gamma_a * Ia
    dIs <- sigma * (1 - p) * E - gamma_s * Is
    dR  <- gamma_a * Ia + gamma_s * Is - omega * R
    
    return(list(c(dS, dE, dIa, dIs, dR)))
  })
}

# Run simulation
output <- ode(y = init, times = times, func = seirs_model, parms = params)
output <- as.data.frame(output)

# Plot results
output_long <- reshape2::melt(output, id = "time")

ggplot(output_long, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +
  labs(title = "SEIRS Malaria Dynamics (Synthetic Data)",
       x = "Time (days)", y = "Population",
       color = "Compartment") +
  theme_minimal()
Addition of the basic SEIRS model
