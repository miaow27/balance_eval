alance_eval <- function(df, grp, cols, spark = FALSE, 
 m_threshold = 0.1, v_threshold = 2, 
 prg = ''){
 
 #' @title evaluate balance between two groups
 #'
 #' @description this function also works for sparklyr dataframe 
 #' 
 #' @param df: data frame 
 #' @param grp: string name for the group variable (has to have value of 1 and 0)
 #' @param cols: a list of all covarites to compare (has to be numeric)
 #' @param spark: if the df is a sparklyr df or not. default (False)
 #' @param m_threshold: how much of SME is considered imbalance (default = 0.1)
 #' @param v_threshold: how much of variation ratio is considered imbalance (default = 2)
 #' @param prg: string name for the prognostic variable. a model is fit on prg ~ cols (among control). then predicted valud is calculated for both group to determine the distance.
 #' 
 #' @return size: count of 
 #' @return balance: summary of normalized difference for overall and each individual covariate
 #' @return plot: the view of the balance
 #' 
 #' @references
 #' https://cran.r-project.org/web/packages/MatchIt/vignettes/assessing-balance.html#:~:text=Assessing%20balance%20involves%20assessing%20whether,joint%20distributional%20balance%20as%20well.
 #' https://pubmed.ncbi.nlm.nih.gov/30788363/
 #' 
 #' @examples
 #' ###############################
 #' ## small data
 #' ###############################
 #' library(tidyverse)
 #' 
 #' iris <- iris %>% mutate(sectosa_flag = if_else(Species == 'setosa', 1, 0))
 #' 
 #' balance_eval(df = iris, grp = 'sectosa_flag', cols = c("Sepal.Length", "Sepal.Width", "Petal.Length","Petal.Width"))
 #' 
 #' ## using Sepa.Length as the prognostic score for balance evaluation
 #' balance_eval(df = iris, grp = 'sectosa_flag', cols = c("Sepal.Width", "Petal.Length","Petal.Width"),
 #' prg = 'Sepal.Length')
 #' 
 #' 
 #' ##############################
 #' ## large data (in spark)
 #' ##############################
 #' library(sparklyr)
 #' library(tidyverse)
 #' sc <- spark_connect()
 #' 
 #' ## subset 
 #' tbl_1 <- dplyr::tbl(sc, 'mydb.mytbl') 
 #' 
 #' ## transform all covariets to numeric 
 #' ## please capped extreme value since we will perform t-statistics on each individual value
 #' tbl_1 <- tbl_1 %>% transmute(
 #' member_id, 
 #' recieve_new_drug_flag,
 #' medical_claim_count, 
 #' risk_of_cancer, 
 #' age, 
 #' emergency_claim_count, 
 #' inpatient_claim_count, 
 #' 
 #' ## need to turn all category to numeric
 #' gender_f = if_else(gender == "F", 1, 0),
 #' 
 #' ## capped the spending
 #' total_medical_spend_capp = if_else(total_medical_spend >= 1e4, 1e4, total_medical_spend),
 #' 
 #' ## capped the pcp cnt
 #' pcp_office_visit_capped = if_else(pcp_office_visit >= 20, 20, pcp_office_visit)
 #' )
 #' 
 #' ## check the transform is successful
 #' ## make sure all confounders are numeric 
 #' tbl_1 %>% head(5)
 #' 
 #' ## check NA
 #' tbl_1 %>% 
 #' mutate_all(is.na) %>%
 #' mutate_all(as.numeric) %>% 
 #' summarise_all(sum) %>% 
 #' as_tibble() %>% 
 #' gather() %>% 
 #' print(n = Inf)
 #' 
 #' ## fill NA
 #' tbl_1 <- tbl_1 %>% replace_na(., 0)
 #' 
 #' # covariate we care about
 #' confounders <- c("gender_f", "age", "risk_of_cancer", "emergency_claim_count", "inpatient_claim_count", "pcp_office_visit_capped")
 #' 
 #' ## check balance between member with vs without bh flag
 #' balance_eval(df = tbl_1, grp = 'recieve_new_drug_flag', cols = confounders, spark = TRUE)
 #' 
 #' ## using total_medical_spend as prognostic score to check overall imbalance
 #' ## ideally, the progostic variable is something we are interested as the outcome
 #' balance_eval(df = tbl_1, grp = 'recieve_new_drug_flag', prg = 'total_medical_spend', cols = confounders, spark = TRUE)



cat('Evaluating balance …')

 ## group size
 .size <- df %>% group_by(!!sym(grp)) %>% count() %>% as_tibble()
 n1 <- .size[.size[grp] == 1, 'n'] %>% pull()
 n0 <- .size[.size[grp] == 0, 'n'] %>% pull()
 
 #################################
 ## calculate propensity score 
 #################################
 f <- paste0(grp, ' ~ ', paste0(cols, collapse = ' + '))
 
 if (spark) {
 lr_model <- ml_logistic_regression(df, formula(f))
df_2 <- ml_predict(lr_model, df) %>% 
 mutate(ps_score = log(probability_1 / (1 - probability_1))) %>% 
 select(-features, -label, -rawPrediction, -probability, -prediction, 
 -probability_0, -probability_1)
} else {
 lr_model <- glm(formula(f), data = df, family = "binomial")
 
 df_2 <- df %>% 
 mutate(ps_score = lr_model$fitted.values,
 ps_score = log(ps_score/(1-ps_score)))
 }
 
 ###################################
 ## calculate prognostic score 
 ###################################
 if (prg != ''){
 f <- paste0(prg, ' ~ ', paste0(cols, collapse = ' + '))
 
 ## determine if prg is binary or continous
 max_val <- df_2 %>% summarise(max(!!sym(prg))) %>% as_tibble() %>% pull()
 min_val <- df_2 %>% summarise(min(!!sym(prg))) %>% as_tibble() %>% pull()
 if ((as.integer(max_val) == 1) & (as.integer(min_val) == 0)){ ## binomial
 
 if (spark) {
 lr_model <- ml_logistic_regression(df_2 %>% filter(!!sym(grp) == 0), formula(f))
 
 df_2 <- ml_predict(lr_model, df_2) %>% 
 mutate(prg_score = log(probability_1 / (1 - probability_1))) %>% 
 select(-features, -label, -rawPrediction, -probability, -prediction, 
 -probability_0, -probability_1)
 
 } else {
 lr_model <- glm(formula(f), data = df_2 %>% filter(!!sym(grp) == 0), family = "binomial")
 
 df_2 <- df_2 %>% 
 mutate(prg_score = predict(lr_model, df_2),
 prg_score = log(prg_score/(1-prg_score)))
 }
 
 }else{ ## gaussian
 
 if (spark) {
 l_model <- ml_linear_regression(df_2 %>% filter(!!sym(grp) == 0), formula(f))
 
 df_2 <- ml_predict(l_model, df_2) %>% 
 mutate(prg_score = prediction) %>% 
 select(-prediction)
 
 } else {
 l_model <- lm(formula(f), data = df_2 %>% filter(!!sym(grp) == 0))
 
 df_2 <- df_2 %>% 
 mutate(prg_score = predict(l_model, df_2))
 }
 
 }
 
 cols <- c('prg_score', cols)
 }
 
 cols <- c('ps_score', cols)
 
 ###################################
 ## compare normalized difference 
 ##################################
 .mean <- df_2 %>% 
 group_by(!!sym(grp)) %>% 
 select(all_of(cols), !!sym(grp)) %>%
 summarise_all(mean) %>% 
 as_tibble() %>% 
 pivot_longer(cols = cols, names_to = 'characteristic', values_to = 'mean') 
 
 .sd <- df_2 %>% 
 group_by(!!sym(grp)) %>% 
 select(all_of(cols), !!sym(grp)) %>%
 summarise_all(sd) %>% 
 as_tibble() %>% 
 pivot_longer(cols = cols, names_to = 'characteristic', values_to = 'sd')
 
 .summary <- .mean %>% 
 left_join(.sd, by = c('characteristic', grp)) %>% 
 pivot_wider(id_cols = 'characteristic', names_from = !!sym(grp), values_from = c('mean', 'sd')) %>% 
 mutate(
 raw_diff = abs(mean_1 - mean_0),
 per_diff = abs(raw_diff / mean_0),
 norm_diff = abs(mean_1 - mean_0)/sqrt((sd_0² + sd_1 ^2)/2), ## using pooled sd (for ATE)
 var_ratio = sd_1 ^2 / sd_0²,
 m_imbalance = if_else(norm_diff >= m_threshold, 'imbalance', paste0('balance: sme < ', m_threshold)),
 v_imbalance = if_else(var_ratio >= v_threshold, 'imbalance', paste0('balance: var ratio < ', v_threshold)),
 imbalance = if_else(m_imbalance == 'imbalance' | v_imbalance == 'imbalance', 'imbalance', 'balance')
 ) %>% 
 arrange(desc(norm_diff))
 
 ## plot the balance
 .plot <- .summary %>% 
 ggplot(aes(x = reorder(characteristic, norm_diff), y = norm_diff, color = imbalance)) + 
 geom_point() + 
 geom_hline(yintercept = m_threshold, linetype = 2) +
 coord_flip() + 
 scale_color_manual(values = c("imbalance" = "#F8766D", "balance" = "#00BFC4")) +
 labs(title = paste0("Balance between ", grp, " = 0 vs 1"))
 
 return(list(size=.size, summary=.summary, plot=.plot))
}
################################ View Package (Example available) ######################
#docstring::docstring(balance_eval)