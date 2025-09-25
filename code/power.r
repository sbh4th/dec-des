library(fixest)    # fast estimation with fixed effects
library(tidyverse)

set.seed(12345)

# Parameters
n_ads     <- 1000      # number of job ads
cv_per_ad <- 4         # 4 CVs per ad
n_sims    <- 500       # number of simulations for power

# True effects
beta_female   <- -0.05   # effect of male vs female
beta_married  <-  0.03   # effect of married vs unmarried
beta_interact <-  0.00   # interaction effect
base_p        <-  0.50   # baseline callback rate

# --- Helper to simulate one data set under matched design ---
simulate_matched <- function() {
  # Balanced assignment
  combos <- expand.grid(gender=c("female","male"), 
    marital=c("unmarried","married"))
  
  df <- do.call(rbind, lapply(1:n_ads, function(ad){
    combos$ad_id <- ad
    combos
  }))
  
  # Outcome probability
  df <- df %>%
    mutate(
      female = as.integer(gender=="female"),
      married = as.integer(marital=="married"),
      prob = base_p + beta_female*female + 
        beta_married*married + 
        beta_interact*female*married,
      callback = rbinom(n(),1,pmin(pmax(prob,0.01),0.99))
    )
  df
}

# --- Helper to simulate one data set under unmatched design ---
simulate_unmatched <- function() {
  df <- expand.grid(ad_id=1:(n_ads * 4))
  
  # random assignment
  df$gender <- sample(c("female","male"), 
                      nrow(df), replace=TRUE)
  df$marital <- sample(c("unmarried","married"),
                       nrow(df), replace=TRUE)
  
  df <- df %>%
    mutate(
      female = as.integer(gender=="female"),
      married = as.integer(marital=="married"),
      prob = base_p + beta_female * female + 
        beta_married * married + 
        beta_interact * female * married,
      callback = rbinom(n(),1,pmin(pmax(prob,0.01),0.99))
    )
  df
}

# --- Estimation ---
estimate_matched <- function(df){
  fit <- feols(callback ~ female * married | ad_id, 
               data=df, cluster="ad_id")
  coef_female <- coef(fit)["female"]
  se_female <- se(fit)["female"]
  p_female <- abs(coef_female/se_female) > 1.96
  coef_married <- coef(fit)["married"]
  se_married <- se(fit)["married"]
  p_married <- abs(coef_married/se_married) > 1.96
  coef_int <- coef(fit)["female:married"]
  se_int <- se(fit)["female:married"]
  p_int <- abs(coef_int/se_int) > 1.96
  c(est_f = coef_female, se_f = se_female,
    p_f = p_female, est_m = coef_married, 
    se_m = se_married, p_m = p_married,
    est_i = coef_int, se_i = se_int, 
    p_i = p_int)
}

estimate_unmatched <- function(df){
  fit <- feols(callback ~ female * married, 
               data=df, vcov="hetero")
  coef_female <- coef(fit)["female"]
  se_female <- se(fit)["female"]
  p_female <- abs(coef_female/se_female) > 1.96
  coef_married <- coef(fit)["married"]
  se_married <- se(fit)["married"]
  p_married <- abs(coef_married/se_married) > 1.96
  coef_int <- coef(fit)["female:married"]
  se_int <- se(fit)["female:married"]
  p_int <- abs(coef_int/se_int) > 1.96
  c(est_f = coef_female, se_f = se_female,
    p_f = p_female, est_m = coef_married, 
    se_m = se_married, p_m = p_married,
    est_i = coef_int, se_i = se_int, 
    p_i = p_int)
}

# --- Simulation loop ---
results <- replicate(n_sims, {
  df_m <- simulate_matched()
  df_u <- simulate_unmatched()
  
  out_m <- estimate_matched(df_m)
  out_u <- estimate_unmatched(df_u)
  
  c(estf_m = out_m["est_f.female"], 
    sef_m = out_m["se_f.female"],
    pf_m = out_m["p_f.female"],
    estm_m = out_m["est_m.married"], 
    sem_m = out_m["se_m.married"],
    pm_m = out_m["p_m.married"],
    esti_m = out_m["est_i.female:married"],
    sei_m = out_m["se_i.female:married"],
    pi_m = out_m["p_i.female:married"],
    estf_u = out_u["est_f.female"], 
    sef_u = out_u["se_f.female"],
    pf_u = out_u["p_f.female"],
    estm_u = out_u["est_m.married"], 
    sem_u = out_u["se_m.married"],
    pm_u = out_u["p_m.married"],
    esti_u = out_u["est_i.female:married"],
    sei_u = out_u["se_i.female:married"],
    pi_u = out_u["p_i.female:married"])
})

# --- Convert to data frame with names ---
results <- as.data.frame(t(results))
colnames(results) <- c("estf_m","sef_m",
  "pf_m", "estm_m", "sem_m", "pm_m", 
  "esti_m", "sei_m", "pi_m", "estf_u", 
  "sef_u", "pf_u", "estm_u", "sem_u", 
  "pm_u", "esti_u", "sei_u", "pi_u")

# --- Compute power ---
alpha <- 0.05
z_crit <- qnorm(1 - alpha/2)

results <- results %>%
  mutate(
    rejectf_m = abs(estf_m / sef_m) > z_crit,
    rejectm_m = abs(estm_m / sem_m) > z_crit,
    rejecti_m = abs(esti_m / sei_m) > z_crit,
    rejectf_u = abs(estf_u / sef_u) > z_crit,
    rejectm_u = abs(estm_u / sem_u) > z_crit,
    rejecti_u = abs(esti_u / sei_u) > z_crit
  )

power_m <- mean(results$reject_m)
power_u <- mean(results$reject_u)

sd_m <- sd(results$est_m)
sd_u <- sd(results$est_u)

cat("Power (matched): ", power_m, "\n")
cat("Power (unmatched): ", power_u, "\n")
cat("Relative efficiency (var ratio unmatched/matched):", (sd_u^2)/(sd_m^2), "\n")

# One simulated dataset for matched design
df_m_example <- simulate_matched()

# One simulated dataset for unmatched design
df_u_example <- simulate_unmatched()


# make a table
tab <-results %>% 
  summarize(across(everything(), mean)) %>% 
  pivot_longer(cols = everything(),
    cols_vary = "slowest",
    names_to = c(".value", "design"),
    names_pattern = "(.*)_(m|u)") %>%
  mutate(design = ifelse(design=="m", 
    "Matched", "Unmatched"))

colnames(tab) <- c("Design", "Est", "SE", "Power", 
 "Est", "SE", "Power", "Est", "SE",
 "Power")

tt(tab, digits = 2) %>%
  group_tt(j = list(
    "Female" = 2:4,
    "Married" = 5:7,
    "Interaction" = 8:10))

