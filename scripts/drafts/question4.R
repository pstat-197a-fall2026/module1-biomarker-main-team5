# Question 4

# Use any method to find:
# a simpler panel that achieves comparable classification accuracy

source('scripts/drafts/question3.R')

## Protein Selection

# The top 2 from multiple testing question-3
proteins_small <- proteins_s1[1:5]
proteins_small

## Logistic Regression
train_small <- train_data %>%
  select(group, any_of(proteins_small)) %>%
  mutate(class = factor(group == 'ASD', levels = c('FALSE','TRUE'))) %>%
  select(-group)

fit_small <- glm(class ~.,
                 data = train_small,
                 family = 'binomial')

test_small <- test_data %>%
  select(group, any_of(proteins_small)) %>%
  mutate(class = factor(group == 'ASD', levels = c('FALSE','TRUE'))) %>%
  select(-group)

# Evaluating accuracy for smaller panel

small_results <- test_small %>%
  add_predictions(fit_small, type = 'response') %>%
  mutate(
    pred_class = factor(pred > 0.5, levels = c('FALSE','TRUE')),
    truth_class = class
  ) %>%
  class_metrics(
    truth = truth_class,
    estimate = pred_class,
    pred,
    event_level = 'second'
  )

small_results