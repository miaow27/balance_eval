library(tidyverse)

# 1. Generate Data  -------------------------------------------------------

N <- 1100 

# age 
age <- rnorm(N, 25, 5) %>% ifelse(. > 100, 100, .) %>% ifelse(. <= 1, 1, .) 

# female 
gender_f <- rbernoulli(N, 0.5) %>% as.integer()

# health risk score 1 - 10 (10 being worst)
health_risk <- rmultinom(N, 1, seq(10, 1)**2)
health_risk <- t(health_risk) %*% as.matrix(seq(1, 10), nrow = 10) 

# prior experience with customer service
prior_cust_service <- rbernoulli(N, 0.05) %>% as.integer()

# prior medical utilization (claim counts)
prior_medical_utilization <- rexp(N, 1/1000)  

# if did regular check up last year
prior_checkup <- rbernoulli(N, 0.15) %>% as.integer()

# combine 
tbl <- tibble(
  member_id = 1:N,
  age = age,
  gender_f = gender_f,
  health_risk = health_risk,
  prior_cust_service = prior_cust_service,
  prior_medical_utilization = prior_medical_utilization,
  prior_checkup = prior_checkup
)


# 2. Split into different group -------------------------------------------

p_holdout <- 1/11
p_outreach <- 0.8 
p_reach <- 0.625
p_engage <- 0.2

holdout <- numeric(N)
holdout[sample(1:N, p_holdout*N)] <- 1
sum(holdout)

target <- 1 - holdout
sum(target)

outreach <- numeric(N)
outreach[sample(which(target == 1), p_outreach * sum(target))] <- 1
sum(outreach)

reach <- numeric(N)
reach[sample(which(outreach == 1), p_reach * sum(outreach))] <- 1
sum(reach)

engage <- numeric(N)
engage[sample(which(reach == 1), p_engage * sum(reach))] <- 1
sum(engage)

tbl$target <- target
tbl$outreach <- outreach
tbl$reach <- reach
tbl$engage <- engage

tbl %>% count(target, outreach, reach, engage)


# 3. Evaluate bias between target and holdout -----------------------------

source('balance_eval.R')

# if we have an outcome of interest
# we can use it to get a composite score assessing balance
# usually better than the ps_score
tbl$next_year_spend = rnorm(n = N, mean = engage * (-100) + age * 10 + gender_f * 50 + health_risk * 200 + prior_checkup * (-100) + (1 - prior_checkup) * engage * (-100), sd = 50)


balance_eval(df = tbl, 
             grp = 'target', 
             cols = c('age', 'gender_f', 'health_risk', 'prior_cust_service', 
                      'prior_medical_utilization', 'prior_checkup'),
             prg = 'next_year_spend')
