library(DeclareDesign)
library(fabricatr)
library(estimatr)
library(dplyr)
library(ggplot2)

set.seed(12345)

# Define parameters
N_ads <- 1000  # Number of job ads
N_cvs_per_ad <- 4  # CVs per ad (2x2 factorial)

# Population model - job ads with heterogeneous discrimination
population <- declare_model(
  ads = add_level(N = N_ads,
                  # Ad-level characteristics
                  sector = sample(c("Finance", "Tech", "Healthcare", "Retail"), N, replace = TRUE),
                  firm_size = sample(c("Small", "Medium", "Large"), N, replace = TRUE, 
                                   prob = c(0.4, 0.4, 0.2)),
                  # Baseline callback rate varies by ad
                  baseline_callback = rnorm(N, mean = 0.15, sd = 0.05),
                  baseline_callback = pmax(0.05, pmin(0.4, baseline_callback)),
                  # Discrimination parameters - vary by ad
                  gender_effect = rnorm(N, mean = -0.03, sd = 0.02),  # Male advantage
                  marriage_effect = rnorm(N, mean = -0.02, sd = 0.015), # Unmarried advantage
                  interaction_effect = rnorm(N, mean = 0.01, sd = 0.01)  # Small interaction
  )
)

# DESIGN 1: BLOCKED/MATCHED DESIGN
# Each ad gets exactly one CV of each type (2x2 factorial)

# Inquiry for blocked design
inquiry_blocked <- declare_inquiry(
  ATE_gender = mean(gender_effect),
  ATE_marriage = mean(marriage_effect), 
  ATE_interaction = mean(interaction_effect)
)

# Data strategy for blocked design - complete factorial within each ad
data_strategy_blocked <- declare_assignment(
  treatment_combo = conduct_ra(N = N, 
                              blocks = ads, 
                              conditions = c("unmarried_female", "unmarried_male", "married_female", "married_male"), 
                              block_m_each = c(1, 1, 1, 1))  # Exactly one of each type per ad
) + 
declare_step(
  handler = function(data) {
    data %>%
      mutate(
        female = case_when(
          treatment_combo %in% c("unmarried_female", "married_female") ~ 1,
          TRUE ~ 0
        ),
        married = case_when(
          treatment_combo %in% c("married_female", "married_male") ~ 1,
          TRUE ~ 0
        )
      )
  },
  label = "create_indicators"
)

# Potential outcomes for blocked design
potential_outcomes_blocked <- declare_potential_outcomes(
  Y ~ baseline_callback + 
      gender_effect * female + 
      marriage_effect * married + 
      interaction_effect * female * married +
      rnorm(n(), 0, 0.02),  # Individual noise
  conditions = list(
    female = c(0, 1),
    married = c(0, 1)
  )
)

# Reveal outcomes
reveal_blocked <- declare_reveal(outcome_variables = Y, assignment_variables = c(female, married))

# Estimator for blocked design - can control for ad fixed effects
estimator_blocked <- declare_estimator(
  Y ~ female * married + factor(ads),
  model = lm_robust,
  clusters = ads,
  inquiry = c("ATE_gender", "ATE_marriage", "ATE_interaction"),
  term = c("female", "married", "female:married"),
  label = "blocked_estimator"
)

# Complete blocked design
blocked_design <- population + inquiry_blocked + data_strategy_blocked + 
  potential_outcomes_blocked + reveal_blocked + estimator_blocked

# DESIGN 2: UNMATCHED DESIGN
# Random assignment of gender, marital status as covariate
# Same total sample size (4000 CVs)

population_unmatched <- declare_model(
  cvs = add_level(N = N_ads * N_cvs_per_ad,
                  ads = rep(1:N_ads, each = N_cvs_per_ad),
                  # CV-level baseline (varies more than in blocked design)
                  cv_noise = rnorm(N, 0, 0.01)
  ),
  ads = add_level(N = N_ads,
                  sector = sample(c("Finance", "Tech", "Healthcare", "Retail"), N, replace = TRUE),
                  firm_size = sample(c("Small", "Medium", "Large"), N, replace = TRUE, 
                                   prob = c(0.4, 0.4, 0.2)),
                  baseline_callback = rnorm(N, mean = 0.15, sd = 0.05),
                  baseline_callback = pmax(0.05, pmin(0.4, baseline_callback)),
                  gender_effect = rnorm(N, mean = -0.03, sd = 0.02),
                  marriage_effect = rnorm(N, mean = -0.02, sd = 0.015),
                  interaction_effect = rnorm(N, mean = 0.01, sd = 0.01)
  ),
  nest = FALSE
)

# Inquiry for unmatched design (same estimands)
inquiry_unmatched <- declare_inquiry(
  ATE_gender = mean(gender_effect),
  ATE_marriage = mean(marriage_effect),
  ATE_interaction = mean(interaction_effect)
)

# Random assignment of gender, marital status as covariate
data_strategy_unmatched <- declare_assignment(
  female = conduct_ra(N = N, prob = 0.5)
) + 
declare_step(
  handler = function(data) {
    data %>%
      mutate(married = rbinom(n(), 1, 0.5))  # Random marital status
  },
  label = "assign_marriage"
)

# Potential outcomes for unmatched design
potential_outcomes_unmatched <- declare_potential_outcomes(
  Y_Z_0 = baseline_callback[ads] + 
          marriage_effect[ads] * married + 
          cv_noise + rnorm(n(), 0, 0.02),
  Y_Z_1 = baseline_callback[ads] + 
          gender_effect[ads] + 
          marriage_effect[ads] * married + 
          interaction_effect[ads] * married +
          cv_noise + rnorm(n(), 0, 0.02)
)

reveal_unmatched <- declare_reveal(Y, female)

# Estimator for unmatched design - includes marital status as covariate
estimator_unmatched <- declare_estimator(
  Y ~ female * married,
  model = lm_robust,
  clusters = ads,
  inquiry = c("ATE_gender", "ATE_marriage", "ATE_interaction"),
  term = c("female", "married", "female:married"),
  label = "unmatched_estimator"
)

# Complete unmatched design
unmatched_design <- population_unmatched + inquiry_unmatched + data_strategy_unmatched + 
                   potential_outcomes_unmatched + reveal_unmatched + estimator_unmatched

# DESIGN DIAGNOSIS - Compare efficiency
print("Diagnosing Blocked Design...")
diagnosis_blocked <- diagnose_design(blocked_design, sims = 30)

print("Diagnosing Unmatched Design...")
diagnosis_unmatched <- diagnose_design(unmatched_design, sims = 500)

# Extract key results
blocked_results <- diagnosis_blocked$diagnosands_df %>%
  select(inquiry, estimator, bias, rmse, power, coverage) %>%
  mutate(design = "Blocked")

unmatched_results <- diagnosis_unmatched$diagnosands_df %>%
  select(inquiry, estimator, bias, rmse, power, coverage) %>%
  mutate(design = "Unmatched")

# Combine results
comparison_results <- bind_rows(blocked_results, unmatched_results)

print("=== DESIGN COMPARISON ===")
print(comparison_results)

# Calculate efficiency gains
efficiency_comparison <- comparison_results %>%
  select(inquiry, design, rmse, power) %>%
  pivot_wider(names_from = design, values_from = c(rmse, power), names_sep = "_") %>%
  mutate(
    rmse_ratio = rmse_Unmatched / rmse_Blocked,
    power_gain = power_Blocked - power_Unmatched,
    efficiency_gain = (rmse_Unmatched^2 - rmse_Blocked^2) / rmse_Unmatched^2
  ) %>%
  select(inquiry, rmse_ratio, power_gain, efficiency_gain)

print("=== EFFICIENCY GAINS FROM BLOCKING ===")
print("RMSE Ratio (Unmatched/Blocked): >1 indicates blocked design is more efficient")
print("Power Gain: Positive values indicate blocked design has higher power")
print("Efficiency Gain: Proportion reduction in MSE from blocking")
print(efficiency_comparison)

# Visualization
comparison_plot <- comparison_results %>%
  ggplot(aes(x = inquiry, y = rmse, fill = design)) +
  geom_col(position = "dodge") +
  labs(title = "RMSE Comparison: Blocked vs Unmatched Design",
       subtitle = "Lower RMSE indicates more efficient estimation",
       x = "Parameter", y = "RMSE", fill = "Design") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(comparison_plot)

# Summary statistics
cat("\n=== DESIGN SUMMARY ===\n")
cat("Blocked Design:\n")
cat("- Sample size: ", N_ads * N_cvs_per_ad, " CVs\n")
cat("- Structure: ", N_cvs_per_ad, " CVs per ad (complete 2x2 factorial)\n")
cat("- Controls for ad fixed effects\n\n")

cat("Unmatched Design:\n")
cat("- Sample size: ", N_ads * N_cvs_per_ad, " CVs\n") 
cat("- Structure: Random assignment of gender, marital status as covariate\n")
cat("- Does not control for ad fixed effects in this specification\n\n")

# Optional: Run single simulation to see data structure
cat("=== SAMPLE DATA STRUCTURE ===\n")
sample_blocked <- draw_data(blocked_design)
sample_unmatched <- draw_data(unmatched_design)

cat("Blocked design data preview:\n")
print(head(sample_blocked[, c("ads", "treatment_combo", "female", "married", "Y")], 8))

cat("\nUnmatched design data preview:\n") 
print(head(sample_unmatched[, c("ads", "female", "married", "Y")], 8))