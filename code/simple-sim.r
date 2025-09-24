library(here)
library(tidyverse)

set.seed(12345)

n_ads <- 1000
cvs_per_ad <- 4
sigma_ad <- 1.0
base_logodds <- -2.0
female_effect <- -0.3
married_effect <- -0.2
true_interaction <- 0.0

data <- data.frame(
  ads = rep(1:n_ads, each = cvs_per_ad),
  female = rep(c(1, 0, 1, 0), times = n_ads),
  married = rep(c(0, 0, 1, 1), times = n_ads)) %>%
  mutate(
    ad_baseline = rnorm(n = n(), 
      mean = base_logodds, sd = sigma_ad),
    callback = rbinom(n = n(), size = 1,
      prob = plogis(ad_baseline +
        female_effect * female +
        married_effect * married +
        true_interaction * female * married))) 

glm <- glm(callback ~ female * married, data = data, family = 'binomial')
fe <- feglm(callback ~ female * married | ads, data = data, family = 'binomial')
re <- glmer(callback ~ female * married + (1 | ads), data = data, family = 'binomial')

modelsummary(list(glm,fe,re))
 