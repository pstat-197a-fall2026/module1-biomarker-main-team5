# Question 3

#Experiment with the following modifications:

#repeat the analysis but carry out the entire selection procedure on a training partition -- 
#in other words, set aside some testing data at the very beginning and don't use it until you are 
#evaluating accuracy at the very end

#choose a larger number (more than ten) of top predictive proteins using each selection method

#use a fuzzy intersection instead of a hard intersection to combine the sets of top predictive proteins across selection methods

#How are results affected by each modification?

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('data/biomarker-clean.RData')


set.seed(123)

## Splitting Train/Testing Datasets

data_split <- initial_split(biomarker_clean, prop = 0.80, strata = group)
train_data <- training(data_split)
test_data  <- testing(data_split)

## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- train_data %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins 
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 25) %>% # selecting larger subset (25 proteins)
  pull(protein)



## RANDOM FOREST
##################

set.seed(123)

# store predictors and response separately for traning data
train_predictors <- train_data %>%
  select(-c(group, ados))

train_response <- train_data %>% pull(group) %>% factor()

rf_out <- randomForest(x = train_predictors, 
                       y = train_response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 25) %>% ## change to n=25
  pull(protein)
proteins_s2


## LOGISTIC REGRESSION
#######################
set.seed(123)

# Fuzzy intersection (union)
proteins_fuzzy <- union(proteins_s1, proteins_s2)

train_sstar <- train_data %>%
  select(group, any_of(proteins_fuzzy)) %>%
  mutate(class = factor(group == 'ASD', levels = c('FALSE', 'TRUE'))) %>%
  select(-group)

# Fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = train_sstar, 
           family = 'binomial')

test_sstar <- test_data %>%
  select(group, any_of(proteins_fuzzy)) %>%
  mutate(class = factor(group == 'ASD', levels = c('FALSE', 'TRUE'))) %>%
  select(-group)

# Evaluating Accuracy
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy, 
                            roc_auc)

test_sstar %>%
  add_predictions(fit, type = 'response') %>%
  mutate(
    pred_class = factor(pred > 0.5, levels = c('FALSE', 'TRUE')),
    truth_class = class
  ) %>%
  class_metrics(
    truth = truth_class,
    estimate = pred_class,
    pred,                    
    event_level = 'second'
  )

Intepretation: 
  
#The overall test set accuracy decreased from 83.9% to 75.0% and ROC AUC from 0.908 to 0.871. 
#The main reason performance dropped is that we fixed data leakage. By picking our proteins using 
#only the training data, the test set stayed completely unseen and provided a truly fair test. 
#Additionally, expanding the feature threshold to n = 25 and using a fuzzy union substantially enlarged 
#the candidate predictor pool, increasing model complexity and leading to slight overfitting in the logistic regression.
