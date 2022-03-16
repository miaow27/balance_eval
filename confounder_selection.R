library(glmnet) 
library(tidyverse)

## tbl is the dataframe contains exposure (A), confounders (W), and outcome (Y) 
## w_confounder is the list of confounders 

## Need to Normalize confounder prior LASSO


normalize <- function(x){
  if (min(x) == max(x)){return(as.numeric(length(x)))} 
  return((x - min(x))/(max(x) - min(x))) 
}

tbl_std <- tbl %>% 
  select(all_of(w_confounder)) %>% 
  mutate_all(normalize) %>% 
  as.matrix() ## Lasso only accept matrix form

# Run cross-validation & select lambda 

mod_cv <- cv.glmnet(
  tbl_std, tbl[, A] %>% as.matrix(), 
  family='binomial', 
  intercept = F, 
  alpha=1, ## lasso 
  nfolds = 30) 

plot(log(mod_cv$lambda), mod_cv$cvm) 
print(paste0("min lambda: ", mod_cv$lambda.min)) 
print(paste0("1se lambda: ", mod_cv$lambda.1se)) 

lambda = mod_cv$lambda.min

# Run optimal lasso with the lamda 
la.eq <- glmnet(tbl_std, 
  tbl[, A] %>% as.matrix(), 
  lambda=lambda, 
  family='binomial', 
  intercept = F, 
  alpha=1 ## lasso 
  ) 
  
  x <- la.eq$beta[,1]
  print(names(x[x > 0])) ## this will be the final selected confounders base on A ~ W models 
  
  ## you can also repeat the similar process for Y ~ A + W 
