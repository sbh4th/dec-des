library(DeclareDesign)
library(fixest)
library(dplyr)

set.seed(2025)

# parameters
n_ads <- 100
cv_per_ad <- 4
beta_gender <- -0.10
beta_married <- 0.10
sigma_ad <- 1.0
base_logodds <- -1.0  # a bit higher to avoid rare-event separation

# population
population <- declare_population(
  ad_id = add_level(
    N = n_ads,
    ad_baseline = rnorm(n = N, mean = base_logodds, sd = sigma_ad)
  ),
  cv_id = add_level(N = cv_per_ad)
)

# potential outcomes
potential_outcomes <- declare_potential_outcomes(
  callback ~ rbinom(
    n = N,
    size = 1,
    prob = plogis(
      ad_baseline +
        beta_gender * (Z_gender == "male") +
        beta_married * (Z_married == "married")
    )
  ),
  conditions = list(
    Z_gender = c("female","male"),
    Z_married = c("unmarried","married")
  )
)

# assignments
assign_matched <- declare_assignment(handler = function(data) {
  combos <- expand.grid(
    Z_gender=c("female","male"), 
    Z_married=c("unmarried","married"))
  assn <- data %>% group_by(ad_id) %>% 
    mutate(row_in_ad=row_number()) %>%
    group_modify(~{
      c <- combos[sample(1:4,4),]
      tibble(Z_gender=c$Z_gender, Z_married=c$Z_married)
    }) %>% 
    ungroup()
  bind_cols(data, assn %>% select(Z_gender,Z_married))
})

assign_unmatched <- declare_assignment(handler=function(data){
  data$Z_gender <- sample(c("female","male"),nrow(data),replace=TRUE)
  data$Z_married <- sample(c("unmarried","married"),nrow(data),replace=TRUE)
  data
})

# reveal
reveal <- declare_reveal(callback, assignment_variables=c("Z_gender","Z_married"))

# inquiry (true value of gender effect)
inquiry <- declare_inquiry(value = beta_gender)

# estimators using fixest
estimator_matched <- declare_estimator(
  handler = label_estimator(function(data){
    data$male <- as.integer(data$Z_gender=="male")
    data$married <- as.integer(data$Z_married=="married")
    fit <- feols(callback ~ male + married | ad_id, data=data, cluster="ad_id")
    est <- coef(fit)["male"]
    if(is.null(est)) est <- NA_real_
    est
  }),
  inquiry = "value"
)

estimator_unmatched <- declare_estimator(
  handler = label_estimator(function(data){
    data$male <- as.integer(data$Z_gender=="male")
    data$married <- as.integer(data$Z_married=="married")
    fit <- feols(callback ~ male + married, data=data, cluster="ad_id")
    est <- coef(fit)["male"]
    if(is.null(est)) est <- NA_real_
    est
  }),
  inquiry = "value"
)

# build designs
design_matched <- population + potential_outcomes + assign_matched + reveal + inquiry + estimator_matched
design_unmatched <- population + potential_outcomes + assign_unmatched + reveal + inquiry + estimator_unmatched

relative_efficiency <- declare_diagnosands(
  # default diagnosands first
  # bias = mean(estimate - estimand),
  # rmse = sqrt(mean((estimate - estimand)^2)),
  # sd_estimate = sd(estimate),
  power = mean(p.value < 0.05),
  coverage = mean(estimand >= conf.low & estimand <= conf.high),
  # our custom diagnosand: variance per design
  # var_estimate = var(estimate),
  .method = "bootstrap"
)

# run diagnosis
diagnosis <- diagnose_design(design_matched, design_unmatched, sims=50)
print(diagnosis)
