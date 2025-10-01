library(fixest)    # fast estimation with fixed effects
library(tidyverse)
library(purrr)
library(tidyr)

set.seed(12345)

# --- Helper to simulate one data set under matched design ---
simulate_matched <- function(
    n_ads = 1000,
    beta_female = -0.05,
    beta_married =  0.03,
    beta_interact = 0.03,
    base_p = 0.15,
    sectors = c("Consulting","Tech","Engineering",
      "Retail", "Management","Marketing",
      "Finance","Healthcare"),
    sector_sd = 0.05,   # sd of sector random intercepts
    ad_sd = 0.02        # sd of ad random intercepts
) {
  # balanced assignment: one row for each ad × (female/male) × (married/unmarried)
  combos <- expand.grid(gender = c("female","male"),
                        marital = c("unmarried","married"),
                        KEEP.OUT.ATTRS = FALSE,
                        stringsAsFactors = FALSE)

  # sector draws (named vector)
  sector_effects <- rnorm(length(sectors), 
    mean = 0, sd = sector_sd)
  names(sector_effects) <- sectors

  # realistic distribution across sectors
  random_values <- runif(length(sectors))
  sector_probs <- random_values / sum(random_values)
  ad_sectors <- sample(sectors, n_ads, 
    replace = TRUE, prob = sector_probs)

  # ad-level random intercepts (length n_ads)
  ad_effects <- rnorm(n_ads, mean = 0, sd = ad_sd)

  # assemble dataframe
  df <- do.call(rbind, lapply(1:n_ads, function(ad) {
    tmp <- combos
    tmp$ad_id <- ad
    tmp$sector <- ad_sectors[ad]
    tmp
  }))

  # compute prob, use correct indexing
  df <- df %>%
    mutate(
      female = as.integer(gender == "female"),
      married = as.integer(marital == "married"),
      sector_eff = sector_effects[sector],        # indexed by the sector column
      ad_eff = ad_effects[ad_id],
      prob = base_p + sector_eff + ad_eff +
             beta_female * female +
             beta_married * married +
             beta_interact * female * married,
      prob = pmin(pmax(prob, 0.001), 0.999),
      callback = rbinom(n(), 1, prob)
    )

  df
}

# Estimation
estimate_matched <- function(df) {
  m <- feols(callback ~ female * married,
    data = df, cluster = "ad_id")
  est <- coef(m)
  se <- sqrt(diag(vcov(m)))
  tibble(term = names(est), est = est, se = se)
}

# Parameter grid
param_grid <- expand.grid(
  beta_female  = c(-0.03, -0.05),
  beta_married = c(-0.03, -0.05),
  beta_interact= c(-0.05, -0.075),
  base_p       = 0.10,
  n_ads        = 2500
)

n_sims <- 1000
alpha <- 0.05
zcrit <- qnorm(1 - alpha/2)

simulate_power_matched_all <- function(params) {
  with(as.list(params), {
    results <- replicate(n_sims, {
      df <- simulate_matched(
        n_ads = n_ads,
        base_p = base_p,
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
      mutate(n_ads = n_ads, base_p = base_p,
             beta_female = beta_female,
             beta_married = beta_married,
             beta_interact = beta_interact,
             )
    summary_df
  })
}

# Apply to all parameter combinations
power_results_matched_all <- map_dfr(split(param_grid, 
  seq(nrow(param_grid))), simulate_power_matched_all)

one <- simulate_power_matched_all(
  params = list(n_ads = 2500,
                base_p = 0.10,
                beta_female = -0.04,
                beta_married = -0.03,
                beta_interact = -0.05,
                n_sims = 1000)
)



ggplot(power_results_matched_all,
       aes(x = beta_married, y = power,
           color = term,
           group = term)) +
  geom_line() +
  geom_point() +
  facet_grid(vars(beta_interact)) +
  labs(y = "Power", x = "Beta (Married)") +
  theme_minimal()


power_results_matched_all %>% 
  dplyr::select(n_ads, base_p, beta_female, beta_married, beta_interact, term, power) %>% 
  tt(, digits = 2) %>%
  style_tt(
    i = seq(1, by = 3, length.out = 32), j = 1:5, 
    rowspan = 3, alignv = "t")




# Concordance
# Step 1: keep only female and male rows
sim_data <- simulate_matched()
sim_fm <- df_small %>%
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

# function to run McNemar per marital status
run_mcnemar <- function(df) {
  # counts
  n_00 <- sum(df$female == 0 & df$male == 0)
  n_10 <- sum(df$female == 0 & df$male == 1)
  n_01 <- sum(df$female == 1 & df$male == 0)
  n_11 <- sum(df$female == 1 & df$male == 1)

  # build 2x2 matrix
  tab <- matrix(c(n_11, n_10, n_01, n_00),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(
                  Female = c("1","0"),   # female callback
                  Male   = c("1","0")    # male callback
                ))

  # run McNemar
  out <- mcnemar.test(tab)
  list(table = tab, mcnemar = out)
}

# Run separately by marital status
res_unmarried <- run_mcnemar(filter(df_wide, married == 0))
res_married   <- run_mcnemar(filter(df_wide, married == 1))

res_unmarried$mcnemar
res_married$mcnemar





