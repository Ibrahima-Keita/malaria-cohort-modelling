################################## Modélisation compartimentale (SIR) #################################
####################################### Modèle SIR de base ##################################

####Exercice pratique en R — Modèle SIR de base
#Objectif : Simuler un modèle SIR sur 160 jours pour observer la dynamique d'une épidémie.

# Charger les bibliothèques nécessaires
library(deSolve)  # pour résoudre les équations différentielles

# Paramètres du modèle
beta <- 0.3       # taux de transmission
gamma <- 0.1      # taux de guérison
N <- 1000         # population totale

# Conditions initiales
S0 <- 999         # Nombre initial de susceptibles
I0 <- 1           # Nombre initial d'infectés
R0 <- 0           # Nombre initial de rétablis

# Mettre tout dans un vecteur
init <- c(S = S0, I = I0, R = R0)

# Créer une fonction pour le système d'équations SIR
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# Créer le vecteur temps
times <- seq(0, 160, by = 1)

# Résoudre le système
out <- ode(y = init, times = times, func = sir_model,
           parms = c(beta = beta, gamma = gamma))

# Convertir le résultat en data.frame
out_df <- as.data.frame(out)

# Tracer les résultats
plot(out_df$time, out_df$S, type = "l", col = "blue", ylim = c(0, 1000),
     xlab = "Temps (jours)", ylab = "Population", lwd = 2,
     main = "Modèle SIR - Évolution d'une épidémie")
lines(out_df$time, out_df$I, col = "red", lwd = 2)
lines(out_df$time, out_df$R, col = "green", lwd = 2)
legend("right", legend = c("Susceptibles", "Infectés", "Rétablis"),
       col = c("blue", "red", "green"), lwd = 2)
