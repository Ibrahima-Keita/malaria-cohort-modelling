##################### Modèle SIR avec saisonnalité sur le paramètre β(t) ##########################
# Charger la bibliothèque
library(deSolve)

# Paramètres
beta0 <- 0.3     # taux de transmission moyen
gamma <- 0.1     # taux de guérison
A <- 0.5         # amplitude saisonnière (entre 0 et 1)
T <- 365         # période saisonnière (1 an)
N <- 1000

# Conditions initiales
init <- c(S = 999, I = 1, R = 0)

# Modèle SIR avec saisonnalité
sir_seasonal <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta_t <- beta0 * (1 + A * sin(2 * pi * time / T))
    
    dS <- -beta_t * S * I / N
    dI <- beta_t * S * I / N - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# Vecteur temps (3 ans)
times <- seq(0, 3*T, by = 1)

# Résolution du modèle
out <- ode(y = init, times = times, func = sir_seasonal,
           parms = c(beta0 = beta0, gamma = gamma, A = A, T = T))

# Conversion en data.frame
out_df <- as.data.frame(out)

# Tracer le graphique
plot(out_df$time, out_df$I, type = "l", col = "red", lwd = 2,
     xlab = "Temps (jours)", ylab = "Nombre d'infectés",
     main = "Modèle SIR avec saisonnalité (I(t))")

lines(out_df$time, out_df$S, col = "blue", lwd = 2)
lines(out_df$time, out_df$R, col = "green", lwd = 2)

legend("right", legend = c("Susceptibles", "Infectés", "Rétablis"),
       col = c("blue", "red", "green"), lwd = 2)
