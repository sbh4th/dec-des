library(fixest)    # fast estimation with fixed effects
library(tidyverse)
library(purrr)
library(tidyr)

set.seed(12345)

# --- Helper to simulate one data set under matched design ---
simulate_matched <- function(n_ads = 1000,
                             cv_per_ad = 4,
                             beta_female = -0.05,
                             beta_married = 0.03,
                             beta_interact = 0.03,
                             base_p = 0.15) {
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

# Estimation
estimate_matched <- function(df) {
  m <- feols(callback ~ female + married + female:married | ad_id,
             data = df, cluster = "ad_id")
  est <- coef(m)
  se <- sqrt(diag(vcov(m)))
  tibble(term = names(est), est = est, se = se)
}

# Parameter grid
param_grid <- expand.grid(
  beta_female  = c(-0.03, -0.05),
  beta_married = c(-0.02, -0.04),
  beta_interact= c(0.00, 0.025),
  base_p       = c(0.075, 0.15),
  n_ads        = c(500, 1000)
)

n_sims <- 50
alpha <- 0.05
zcrit <- qnorm(1 - alpha/2)

simulate_power_matched_all <- function(params) {
  with(as.list(params), {
    results <- replicate(n_sims, {
      df <- simulate_matched(
        n_ads = n_ads,
        beta_female = beta_female,
        beta_married = beta_married,
        beta_interact = beta_interact
      )
      est_df <- estimate_matched(df)
      
      # Compute t-stats and reject for each coefficient
      est_df <- est_df %>%
        mutate(tstat = est/se,
               reject = abs(tstat) > zcrit)
      est_df
    }, simplify = FALSE)
    
    # Combine all simulation results
    results_df <- bind_rows(results, .id = "sim")
    
    # Summarize power for each term
    summary_df <- results_df %>%
      group_by(term) %>%
      summarise(
        mean_est = mean(est),
        mc_se    = sd(est),
        # mean_se  = mean(se),
        power    = mean(reject),
        .groups = "drop"
      ) %>%
      mutate(beta_female = beta_female,
             beta_married = beta_married,
             beta_interact = beta_interact,
             n_ads = n_ads)
    summary_df
  })
}

# Apply to all parameter combinations
power_results_matched_all <- map_dfr(split(param_grid, 
  seq(nrow(param_grid))), simulate_power_matched_all)

simulate_power_matched_all(
  params = list(n_ads = 1000,
                beta_female = -0.05,
                beta_married = 0.03,
                beta_interact = 0.02)
)







# Concordance
# Step 1: keep only female and male rows
sim_data <- simulate_matched()
sim_fm <- sim1 %>%
  filter(gender %in% c("female", "male"))

# Step 2: pivot to wide for each ad_id × marital status
df_wide <- sim_fm %>%
  select(ad_id, married, gender, callback) %>%
  pivot_wider(
    names_from = gender,
    values_from = callback
  )

# Now df_wide has: ad_id, married (0/1), female, male

# Step 3: summarise counts per marital status
totals_by_married <- df_wide %>%
  group_by(married) %>%
  summarise(
    n_00 = sum(female == 0 & male == 0),
    n_10 = sum(female == 0 & male == 1),
    n_01 = sum(female == 1 & male == 0),
    n_11 = sum(female == 1 & male == 1),
    .groups = "drop"
  )

totals_by_married