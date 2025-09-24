library(DeclareDesign)
library(fabricatr)
library(estimatr)
library(dplyr)
library(ggplot2)
library(tidyr)

set.seed(12345)

# Define parameters
N_ads <- 1000  # Number of job ads
N_cvs_per_ad <- 4  # CVs per ad (2x2 factorial)

# True effect parameters (constant across all ads for simplicity)
true_gender_effect <- -0.03   # Male advantage (negative = discrimination against females)
true_marriage_effect <- -0.02  # Unmarried advantage  
true_interaction <- 0.01       # Small interaction

# DESIGN 1: BLOCKED/MATCHED DESIGN
# Each ad gets exactly one CV of each type

# Population: Start with CV-level data
population_blocked <- declare_model(
  N = N_ads * N_cvs_per_ad,
  ads = rep(1:N_ads, each = N_cvs_per_ad),
  # Each group of 4 CVs (per ad) gets the full factorial
  cv_within_ad = rep(1:4, times = N_ads),
  female = rep(c(1, 0, 1, 0), times = N_ads),
  married = rep(c(0, 0, 1, 1), times = N_ads),
  # Add ad-level heterogeneity in baseline callback rates
  baseline_rate = rep(rnorm(N_ads, 0.15, 0.05), each = N_cvs_per_ad),
  baseline_rate = pmax(0.05, pmin(0.4, baseline_rate))
)

# Inquiry: What we want to estimate
inquiry_blocked <- declare_inquiry(
  ATE_gender = true_gender_effect,
  ATE_marriage = true_marriage_effect,
  ATE_interaction = true_interaction
)

# No assignment step needed - treatments are built into population

# Potential outcomes
potential_outcomes_blocked <- declare_potential_outcomes(
  Y ~ baseline_rate + 
      true_gender_effect * female + 
      true_marriage_effect * married + 
      true_interaction * female * married + 
      rnorm(N, 0, 0.02)
)

# Reveal (deterministic - no assignment to reveal)
reveal_blocked <- declare_step(
  handler = function(data) {
    # The Y column should already exist from potential outcomes
    return(data)
  },
  label = "reveal_Y"
)

# Estimator with ad fixed effects
estimator_blocked <- declare_estimator(
  Y ~ female + married + female:married + factor(ads),
  .method = lm_robust,
  clusters = ads,
  inquiry = c("ATE_gender", "ATE_marriage", "ATE_interaction"),
  term = c("female", "married", "female:married"),
  label = "blocked_with_FE"
)

# Complete blocked design
blocked_design <- population_blocked + inquiry_blocked + 
                 potential_outcomes_blocked + reveal_blocked + estimator_blocked

# DESIGN 2: UNMATCHED DESIGN
# Same total sample size, but random assignment

# Population for unmatched design
population_unmatched <- declare_model(
  N = N_ads * N_cvs_per_ad,
  ads = rep(1:N_ads, each = N_cvs_per_ad),
  # Marital status is observed covariate (not randomized)
  married = rbinom(N, 1, 0.5),
  # Add ad-level heterogeneity
  baseline_rate = rep(rnorm(N_ads, 0.15, 0.05), each = N_cvs_per_ad),
  baseline_rate = pmax(0.05, pmin(0.4, baseline_rate))
)

# Same inquiry
inquiry_unmatched <- declare_inquiry(
  ATE_gender = true_gender_effect,
  ATE_marriage = true_marriage_effect,
  ATE_interaction = true_interaction
)

# Random assignment of gender
assignment_unmatched <- declare_assignment(
  female = conduct_ra(N = N, prob = 0.5)
)

# Potential outcomes
potential_outcomes_unmatched <- declare_potential_outcomes(
  Y_female_0 = baseline_rate + 
               true_marriage_effect * married + 
               rnorm(N, 0, 0.02),
  Y_female_1 = baseline_rate + 
               true_gender_effect + 
               true_marriage_effect * married + 
               true_interaction * married + 
               rnorm(N, 0, 0.02)
)

# Reveal
reveal_unmatched <- declare_reveal(Y, female)

# Estimator (without ad fixed effects to show difference)
estimator_unmatched <- declare_estimator(
  Y ~ female + married + female:married,
  .method = lm_robust,
  clusters = ads,
  inquiry = c("ATE_gender", "ATE_marriage", "ATE_interaction"),
  term = c("female", "married", "female:married"),
  label = "unmatched_no_FE"
)

# Complete unmatched design
unmatched_design <- population_unmatched + inquiry_unmatched + assignment_unmatched +
                   potential_outcomes_unmatched + reveal_unmatched + estimator_unmatched

# DIAGNOSE DESIGNS
cat("Diagnosing Blocked Design...\n")
diagnosis_blocked <- diagnose_design(blocked_design, sims = 30)

cat("Diagnosing Unmatched Design...\n")
diagnosis_unmatched <- diagnose_design(unmatched_design, sims = 500)

# Extract and compare results
blocked_results <- diagnosis_blocked$diagnosands_df %>%
  select(inquiry, bias, rmse, power, coverage) %>%
  mutate(design = "Blocked")

unmatched_results <- diagnosis_unmatched$diagnosands_df %>%
  select(inquiry, bias, rmse, power, coverage) %>%
  mutate(design = "Unmatched")

comparison_results <- bind_rows(blocked_results, unmatched_results)

cat("\n=== DESIGN COMPARISON ===\n")
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

cat("\n=== EFFICIENCY GAINS FROM BLOCKING ===\n")
cat("RMSE Ratio (Unmatched/Blocked): >1 indicates blocked design is more efficient\n")
cat("Power Gain: Positive values indicate blocked design has higher power\n")
cat("Efficiency Gain: Proportion reduction in MSE from blocking\n")
print(efficiency_comparison)

# Quick data preview
cat("\n=== SAMPLE DATA PREVIEW ===\n")
sample_blocked <- draw_data(blocked_design)
cat("Blocked design (first 12 rows):\n")
print(head(sample_blocked[, c("ads", "female", "married", "Y")], 12))

cat("\nUnmatched design (first 8 rows):\n")
sample_unmatched <- draw_data(unmatched_design)
print(head(sample_unmatched[, c("ads", "female", "married", "Y")], 8))

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