########################### Modèle SIRS avec saisonnalité et perte d’immunité ##########################
library(deSolve)

# Paramètres
beta0 <- 0.3       # Transmission moyenne
gamma <- 0.1       # Guérison
delta <- 0.01      # Perte d’immunité
A <- 0.4           # Amplitude saisonnière
T <- 365           # Période (1 an)
N <- 1000

# Conditions initiales
init <- c(S = 990, I = 10, R = 0)

# Modèle SIRS avec saisonnalité
sirs_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta_t <- beta0 * (1 + A * sin(2 * pi * time / T))
    
    dS <- -beta_t * S * I / N + delta * R
    dI <- beta_t * S * I / N - gamma * I
    dR <- gamma * I - delta * R
    
    return(list(c(dS, dI, dR)))
  })
}

# Temps
times <- seq(0, 3*T, by = 1)

# Résolution
out <- ode(y = init, times = times, func = sirs_model,
           parms = c(beta0 = beta0, gamma = gamma, delta = delta, A = A, T = T))

out_df <- as.data.frame(out)

# Tracer
plot(out_df$time, out_df$I, type = "l", col = "red", lwd = 2,
     xlab = "Temps (jours)", ylab = "Infectés", main = "Modèle SIRS avec saisonnalité")
lines(out_df$time, out_df$S, col = "blue", lwd = 2)
lines(out_df$time, out_df$R, col = "green", lwd = 2)
legend("right", legend = c("Susceptibles", "Infectés", "Rétablis"),
       col = c("blue", "red", "green"), lwd = 2)


#### Ajout des naissances et décès naturels dans le modèle

library(deSolve)

# Paramètres
beta0 <- 0.3
gamma <- 0.1
delta <- 0.01
A <- 0.4
T <- 365
mu <- 1 / (60*365)  # espérance de vie de 60 ans
N <- 1000

# Conditions initiales
init <- c(S = 990, I = 10, R = 0)

# Modèle SIRS avec natalité et mortalité
sirs_birth_death <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta_t <- beta0 * (1 + A * sin(2 * pi * time / T))
    
    dS <- mu*N - beta_t * S * I / N + delta * R - mu * S
    dI <- beta_t * S * I / N - gamma * I - mu * I
    dR <- gamma * I - delta * R - mu * R
    
    return(list(c(dS, dI, dR)))
  })
}

# Simulation
times <- seq(0, 3*T, by = 1)

out <- ode(y = init, times = times, func = sirs_birth_death,
           parms = c(beta0 = beta0, gamma = gamma, delta = delta, A = A, T = T, mu = mu, N = N))

out_df <- as.data.frame(out)

# Tracer
plot(out_df$time, out_df$I, type = "l", col = "red", lwd = 2,
     xlab = "Temps (jours)", ylab = "Infectés",
     main = "SIRS avec natalité et mortalité")
lines(out_df$time, out_df$S, col = "blue", lwd = 2)
lines(out_df$time, out_df$R, col = "green", lwd = 2)
legend("right", legend = c("Susceptibles", "Infectés", "Rétablis"),
       col = c("blue", "red", "green"), lwd = 2)

####### Ajout des naissances et décès naturels dans le modèle

library(deSolve)

# Paramètres
beta0 <- 0.3
gamma <- 0.1
delta <- 0.01
A <- 0.4
T <- 365
mu <- 1 / (60*365)  # espérance de vie de 60 ans
N <- 1000

# Conditions initiales
init <- c(S = 990, I = 10, R = 0)

# Modèle SIRS avec natalité et mortalité
sirs_birth_death <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta_t <- beta0 * (1 + A * sin(2 * pi * time / T))
    
    dS <- mu*N - beta_t * S * I / N + delta * R - mu * S
    dI <- beta_t * S * I / N - gamma * I - mu * I
    dR <- gamma * I - delta * R - mu * R
    
    return(list(c(dS, dI, dR)))
  })
}

# Simulation
times <- seq(0, 3*T, by = 1)

out <- ode(y = init, times = times, func = sirs_birth_death,
           parms = c(beta0 = beta0, gamma = gamma, delta = delta, A = A, T = T, mu = mu, N = N))

out_df <- as.data.frame(out)

# Tracer
plot(out_df$time, out_df$I, type = "l", col = "red", lwd = 2,
     xlab = "Temps (jours)", ylab = "Infectés",
     main = "SIRS avec natalité et mortalité")
lines(out_df$time, out_df$S, col = "blue", lwd = 2)
lines(out_df$time, out_df$R, col = "green", lwd = 2)
legend("right", legend = c("Susceptibles", "Infectés", "Rétablis"),
       col = c("blue", "red", "green"), lwd = 2)

######### simulation d’interventions (ITNs et vaccination) dans le modèle SIRS
######### Ajout d’interventions (ITNs / traitement / vaccination)
#### Code R complet du modèle avec ITNs + vaccination

library(deSolve)

# Paramètres
beta0 <- 0.3
gamma <- 0.1
delta <- 0.01
mu <- 1 / (60 * 365)  # mortalité naturelle (espérance de vie de 60 ans)
A <- 0.4              # amplitude de la saisonnalité
T <- 365              # période saisonnière
v <- 0.2              # 20% des nouveau-nés sont vaccinés
N <- 1000

# Durée de protection ITNs (en jours)
itn_start <- 180      # jour 180 = début de saison des pluies
itn_duration <- 90    # 90 jours de réduction de beta
itn_reduction <- 0.5  # réduction de 50%

# Conditions initiales
init <- c(S = 990, I = 10, R = 0)

# Modèle avec interventions
sirs_interventions <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Réduction de beta pendant la période ITN
    beta_eff <- beta0 * (1 + A * sin(2 * pi * time / T))
    if (time %% T >= itn_start && time %% T <= (itn_start + itn_duration)) {
      beta_eff <- beta_eff * (1 - itn_reduction)
    }
    
    dS <- (1 - v) * mu * N - beta_eff * S * I / N + delta * R - mu * S
    dI <- beta_eff * S * I / N - gamma * I - mu * I
    dR <- v * mu * N + gamma * I - delta * R - mu * R
    
    return(list(c(dS, dI, dR)))
  })
}

# Simulation
times <- seq(0, 3*T, by = 1)

out <- ode(y = init, times = times, func = sirs_interventions,
           parms = c(beta0 = beta0, gamma = gamma, delta = delta,
                     mu = mu, A = A, T = T, v = v, N = N,
                     itn_start = itn_start, itn_duration = itn_duration,
                     itn_reduction = itn_reduction))

out_df <- as.data.frame(out)

# Tracer les résultats
plot(out_df$time, out_df$I, type = "l", col = "red", lwd = 2,
     xlab = "Temps (jours)", ylab = "Infectés",
     main = "Modèle SIRS avec ITNs et vaccination")
lines(out_df$time, out_df$S, col = "blue", lwd = 2)
lines(out_df$time, out_df$R, col = "green", lwd = 2)
legend("right", legend = c("Susceptibles", "Infectés", "Rétablis"),
       col = c("blue", "red", "green"), lwd = 2)









