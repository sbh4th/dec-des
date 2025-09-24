library(DeclareDesign)
library(estimatr)

set.seed(12345)

N_ads <- 1000
N_cvs_per_ad <- 4
true_gender_effect <- -0.03
true_marriage_effect <- -0.02
true_interaction <- 0.01

population_blocked <- declare_model(
  N = N_ads * N_cvs_per_ad,
  ads = rep(1:N_ads, each = N_cvs_per_ad),
  female = rep(c(1, 0, 1, 0), times = N_ads),
  married = rep(c(0, 0, 1, 1), times = N_ads),
  baseline_rate = rep(rnorm(N_ads, 0.15, 0.05), each = N_cvs_per_ad),
  baseline_rate = pmax(0.05, pmin(0.4, baseline_rate)),
  Y = baseline_rate + 
      true_gender_effect * female + 
      true_marriage_effect * married + 
      true_interaction * female * married + 
      rnorm(N, 0, 0.02)
)

inquiry_blocked <- declare_inquiry(
  ATE_gender = true_gender_effect,
  ATE_marriage = true_marriage_effect,
  ATE_interaction = true_interaction
)

estimator_blocked <- declare_estimator(
  Y ~ female + married + female:married + factor(ads),
  .method = lm_robust,
  clusters = ads,
  inquiry = c("ATE_gender", "ATE_marriage", "ATE_interaction"),
  term = c("female", "married", "female:married"),
  label = "blocked_with_FE"
)

blocked_design <- population_blocked + inquiry_blocked + estimator_blocked

diagnosis_blocked <- diagnose_design(blocked_design, sims = 30)
diagnosis_blocked
