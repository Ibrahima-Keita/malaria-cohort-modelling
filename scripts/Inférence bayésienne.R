############## estimation descriptive (fréquentiste) des paramètres du modèle SEIRS ##########
# estimate_SEIRS_parameters.R

# estimate_SEIRS_parameters.R
# Estimation descriptive des paramètres p, gamma_a et gamma_s

library(tidyverse)

# Charger les données
data <- read_csv("data/synthetic_malaria_cohort.csv")

# Classifier les états d'infection
data <- data %>%
  mutate(InfectionStatus = case_when(
    RDT_Result == "Positive" & Fever == 1 ~ "Is",  # Symptomatique
    RDT_Result == "Positive" & Fever == 0 ~ "Ia",  # Asymptomatique
    TRUE ~ "S"  # Négatif
  ))

#  Estimation de p (proportion d’asymptomatiques parmi les infectés)
infected <- data %>% filter(InfectionStatus %in% c("Ia", "Is"))
n_total_infected <- nrow(infected)
n_asymptomatic <- sum(infected$InfectionStatus == "Ia")
p_hat <- n_asymptomatic / n_total_infected

#  Estimation de la durée moyenne d’infection (pour gamma)
data <- data %>% arrange(ChildID, VisitMonth)

# Créer un identifiant de groupe pour chaque épisode d'infection
data <- data %>%
  group_by(ChildID) %>%
  mutate(InfectionGroup = cumsum(InfectionStatus != lag(InfectionStatus, default = first(InfectionStatus)))) %>%
  ungroup()

# Durée des épisodes d'infection (Ia et Is)
durations <- data %>%
  filter(InfectionStatus %in% c("Ia", "Is")) %>%
  group_by(ChildID, InfectionGroup, InfectionStatus) %>%
  summarise(Duration = n(), .groups = "drop")

# Moyennes des durées par type
mean_durations <- durations %>%
  group_by(InfectionStatus) %>%
  summarise(MeanDuration = mean(Duration), .groups = "drop")

# Estimation des taux gamma
gamma_estimates <- mean_durations %>%
  mutate(Gamma = round(1 / MeanDuration, 3))

# Résumé final
cat("-------- Résumé des estimations descriptives --------\n")
cat(paste0("Proportion asymptomatique (p) : ", round(p_hat, 3), "\n\n"))
print(gamma_estimates)


##################################### Inférence bayésienne ###################################

#### Étape 1 : Installation des packages dans R

#install.packages("brms")
#install.packages("tidyverse")

#### Étape 2 : Script complet avecbrms

# estimate_p_brms.R
# Bayesian estimation of asymptomatic proportion using brms

library(brms)
library(tidyverse)

# Filter only infected children (RDT positive)
infected <- data %>%
  filter(RDT_Result == "Positive") %>%
  mutate(
    Asymptomatic = ifelse(Fever == 0, 1, 0)
  )

# Summarize as a binomial outcome
n_total <- nrow(infected)
n_asymptomatic <- sum(infected$Asymptomatic)

asympt_data <- tibble(
  asympt_success = n_asymptomatic,
  asympt_trials = n_total
)

# Bayesian model (binomial with logit link)
fit <- brm(
  asympt_success | trials(asympt_trials) ~ 1,
  data = asympt_data,
  family = binomial(),
  prior = prior(beta(1, 1), class = "Intercept"),
  iter = 2000,
  chains = 4,
  seed = 123
)

# Summary of posterior
print(fit)

# Plot posterior distribution
plot(fit)

# Extract posterior samples and credible interval
posterior_samples <- posterior_summary(fit, probs = c(0.025, 0.975))
posterior_samples

#### Ajouter une visualisation avecggplot2

# Load posterior samples and transform logit scale to probability
posterior_df <- as_draws_df(fit)
posterior_df$p <- plogis(posterior_df$b_Intercept)  # transform log-odds to probability

# Plot posterior distribution
library(ggplot2)

ggplot(posterior_df, aes(x = p)) +
  geom_density(fill = "skyblue", alpha = 0.6, color = "darkblue") +
  geom_vline(xintercept = quantile(posterior_df$p, c(0.025, 0.975)), 
             linetype = "dashed", color = "red") +
  geom_vline(xintercept = mean(posterior_df$p), 
             linetype = "solid", color = "blue") +
  labs(
    title = "Posterior Distribution of Asymptomatic Proportion (p)",
    subtitle = "Bayesian estimation using brms",
    x = "Proportion of asymptomatic cases (p)",
    y = "Density"
  ) +
  theme_minimal()
