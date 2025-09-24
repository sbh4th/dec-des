# Install packages if necessary:
# install.packages(c("DeclareDesign","randomizr","estimatr"))

library(DeclareDesign)
library(dplyr)
library(fixest)

set.seed(2025)

# ---- parameters ----
n_ads <- 1000            # number of ads (blocks)
cv_per_ad <- 4           # 4 CVs per ad (complete set of 2x2 treatments)
beta_gender <- -0.30     # log-odds effect of being male
beta_married <- 0.10     # log-odds effect of being married
beta_interaction <- 0.00 # interaction effect
sigma_ad <- 1.0          # SD of ad-level baseline log-odds
base_logodds <- -1.0     # baseline intercept on logit scale

# -----------------------------------------------
# DeclarePopulation
# -----------------------------------------------
# ---- population ----
population <- declare_population(
  ad_id = add_level(
    N = n_ads,
    ad_baseline = rnorm(n = N, mean = base_logodds, sd = sigma_ad)
  ),
  cv_id = add_level(
    N = cv_per_ad   # 4 CVs per ad
  )
)

# -----------------------------------------------
# Potential Outcomes
# -----------------------------------------------
# We use a custom function to generate potential outcomes
# for each (gender, married) combination:
potential_outcomes <- declare_potential_outcomes(
  callback ~ rbinom(
    n = N,
    size = 1,
    prob = plogis(
      ad_baseline +
        beta_gender * (Z_gender == "male") +
        beta_married * (Z_married == "married") +
        beta_interaction * (Z_gender == "male") * (Z_married == "married")
    )
  ),
  conditions = list(
    Z_gender = c("female","male"),
    Z_married = c("unmarried","married")
  )
)

# -----------------------------------------------
# Assignment: matched vs unmatched
# -----------------------------------------------

# Matched assignment: exactly one of each combo per ad
assign_matched <- declare_assignment(handler = function(data) {
  combos <- expand.grid(
    Z_gender = c("female","male"),
    Z_married = c("unmarried","married"),
    stringsAsFactors = FALSE
  )
  # replicate for each ad
  assn <- data %>%
    group_by(ad_id) %>%
    mutate(row_in_ad = row_number()) %>%
    group_modify(~{
      c <- combos[sample(1:4, 4),]  # random order per ad
      tibble(Z_gender = c$Z_gender, Z_married = c$Z_married)
    }) %>% ungroup()
  bind_cols(data, assn %>% select(Z_gender,Z_married))
})

# Unmatched assignment: random across all CVs
assign_unmatched <- declare_assignment(handler = function(data) {
  data$Z_gender <- sample(c("female","male"), nrow(data), replace = TRUE)
  data$Z_married <- sample(c("unmarried","married"), nrow(data), replace = TRUE)
  data
})

# -----------------------------------------------
# Reveal outcomes
# -----------------------------------------------
reveal <- declare_reveal(
  outcome = callback,
  assignment_variables = c("Z_gender","Z_married")
)

# -----------------------------------------------
# Estimators
# -----------------------------------------------
# Matched: include ad fixed effects (factor(ad_id))
estimator_matched <- declare_estimator(
  handler = label_estimator(function(data) {
    data$male <- as.integer(data$Z_gender == "male")
    data$married <- as.integer(data$Z_married == "married")
    # Linear probability model with ad fixed effects
    fit <- feols(callback ~ male + married | ad_id, data = data, cluster = "ad_id")
    coef(fit)["male"]
  }),
  inquiry = declare_inquiry(beta_gender)
)

estimator_unmatched <- declare_estimator(
  handler = label_estimator(function(data) {
    data$male <- as.integer(data$Z_gender == "male")
    data$married <- as.integer(data$Z_married == "married")
    # No ad fixed effects here
    fit <- feols(callback ~ male + married, data = data, cluster = "ad_id")
    coef(fit)["male"]
  }),
  inquiry = declare_inquiry(beta_gender)
)

# -----------------------------------------------
# Build two designs
# -----------------------------------------------
design_matched <- population +
  potential_outcomes +
  assign_matched +
  reveal +
  estimator_matched

design_unmatched <- population +
  potential_outcomes +
  assign_unmatched +
  reveal +
  estimator_unmatched

# -----------------------------------------------
# Diagnose both designs
# -----------------------------------------------
diagnosis <- diagnose_design(
  design_matched,
  design_unmatched,
  sims = 30  # increase to 1000+ for more precision
)

print(diagnosis)

# Look at the table
diagnosis$diagnosands_df
