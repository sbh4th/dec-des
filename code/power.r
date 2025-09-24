library(fixest)    # fast estimation with fixed effects
library(dplyr)

set.seed(12345)

# Parameters
n_ads <- 1000        # number of job ads
cv_per_ad <- 4       # 4 CVs per ad
n_sims <- 100       # number of simulations for power

# True effects
beta_gender  <- -0.03   # effect of male vs female
beta_married <- -0.02   # effect of married vs unmarried
beta_interact <- 0.00   # interaction effect
base_p <- 0.15          # baseline callback rate

# --- Helper to simulate one data set under matched design ---
simulate_matched <- function() {
  # Balanced assignment
  combos <- expand.grid(gender=c("female","male"), married=c("unmarried","married"))
  
  df <- do.call(rbind, lapply(1:n_ads, function(ad){
    combos$ad_id <- ad
    combos
  }))
  
  # Outcome probability
  df <- df %>%
    mutate(
      male = as.integer(gender=="male"),
      married01 = as.integer(married=="married"),
      prob = base_p + beta_gender*male + beta_married*married01 + beta_interact*male*married01,
      callback = rbinom(n(),1,pmin(pmax(prob,0.01),0.99))
    )
  df
}

# --- Helper to simulate one data set under unmatched design ---
simulate_unmatched <- function() {
  df <- expand.grid(cv=1:cv_per_ad, ad_id=1:n_ads)
  
  # random assignment
  df$gender <- sample(c("female","male"), nrow(df), replace=TRUE)
  df$married <- sample(c("unmarried","married"), nrow(df), replace=TRUE)
  
  df <- df %>%
    mutate(
      male = as.integer(gender=="male"),
      married01 = as.integer(married=="married"),
      prob = base_p + beta_gender*male + beta_married*married01 + beta_interact*male*married01,
      callback = rbinom(n(),1,pmin(pmax(prob,0.01),0.99))
    )
  df
}

# --- Estimation ---
estimate_matched <- function(df){
  fit <- feols(callback ~ male + married01 | ad_id, data=df, cluster="ad_id")
  coef(fit)["male"]
}

estimate_unmatched <- function(df){
  fit <- feols(callback ~ male + married01, data=df, cluster="ad_id")
  coef(fit)["male"]
}

# --- Simulation loop ---
results <- replicate(n_sims, {
  df_m <- simulate_matched()
  df_u <- simulate_unmatched()
  est_m <- estimate_matched(df_m)
  est_u <- estimate_unmatched(df_u)
  
  # also SEs for power
  se_m <- se(feols(callback ~ male + married01 | ad_id, data=df_m, cluster="ad_id"))["male"]
  se_u <- se(feols(callback ~ male + married01, data=df_u, cluster="ad_id"))["male"]
  
  c(est_m=est_m,se_m=se_m, est_u=est_u,se_u=se_u)
})

# results <- t(results)
# results <- as.data.frame(results)

# --- Convert to data frame with names ---
results <- as.data.frame(t(results))
colnames(results) <- c("est_m","se_m","est_u","se_u")

# --- Compute power ---
alpha <- 0.05
z_crit <- qnorm(1-alpha/2)

results <- results %>%
  mutate(
    t_m = est_m/se_m,
    t_u = est_u/se_u,
    reject_m = abs(t_m) > z_crit,
    reject_u = abs(t_u) > z_crit
  )

power_m <- mean(results$reject_m, na.rm=TRUE)
power_u <- mean(results$reject_u, na.rm=TRUE)

sd_m <- sd(results$est_m, na.rm=TRUE)
sd_u <- sd(results$est_u, na.rm=TRUE)

cat("Power (matched):  ", round(power_m,3), "\n")
cat("Power (unmatched):", round(power_u,3), "\n")
cat("SD of est. gender effect (matched):  ", round(sd_m,4), "\n")
cat("SD of est. gender effect (unmatched):", round(sd_u,4), "\n")
cat("Relative efficiency (var ratio unmatched/matched):", round((sd_u^2)/(sd_m^2),2), "\n")

# One simulated dataset for matched design
df_m_example <- simulate_matched()

# One simulated dataset for unmatched design
df_u_example <- simulate_unmatched()
