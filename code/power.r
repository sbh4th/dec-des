library(fixest)    # fast estimation with fixed effects
library(dplyr)

set.seed(12345)

# Parameters
n_ads     <- 1000      # number of job ads
cv_per_ad <- 4         # 4 CVs per ad
n_sims    <- 1000       # number of simulations for power

# True effects
beta_female   <- -0.03   # effect of male vs female
beta_married  <-  0.02   # effect of married vs unmarried
beta_interact <-  0.00   # interaction effect
base_p        <-  0.15   # baseline callback rate

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
  df <- expand.grid(ad_id=1:n_ads)
  
  # random assignment
  df$gender <- sample(c("female","male"), nrow(df), replace=TRUE)
  df$marital <- sample(c("unmarried","married"), nrow(df), replace=TRUE)
  
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
  c(est = coef_female, se = se_female)
}

estimate_unmatched <- function(df){
  fit <- feols(callback ~ female * married, 
               data=df, vcov="hetero")
  coef_female <- coef(fit)["female"]
  se_female <- se(fit)["female"]
  c(est = coef_female, se = se_female)
}

# --- Simulation loop ---
results <- replicate(n_sims, {
  df_m <- simulate_matched()
  df_u <- simulate_unmatched()
  
  out_m <- estimate_matched(df_m)
  out_u <- estimate_unmatched(df_u)
  
  c(est_m = out_m["est.female"], se_m = out_m["se.female"],
    est_u = out_u["est.female"], se_u = out_u["se.female"])
})

# --- Convert to data frame with names ---
results <- as.data.frame(t(results))
colnames(results) <- c("est_m","se_m","est_u","se_u")

# --- Compute power ---
alpha <- 0.05
z_crit <- qnorm(1 - alpha/2)

results <- results %>%
  mutate(
    t_m = est_m / se_m,
    t_u = est_u / se_u,
    reject_m = abs(t_m) > z_crit,
    reject_u = abs(t_u) > z_crit
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
results %>% summarize(across(everything(), mean)) %>% pivot_longer(
  cols = everything(),
  names_to = c(".value", "design"),
  names_pattern = "(.*)_(m|u)"
)

