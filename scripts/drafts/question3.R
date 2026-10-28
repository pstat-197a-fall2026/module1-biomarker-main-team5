# Question 3

#Experiment with the following modifications:

#repeat the analysis but carry out the entire selection procedure on a training partition -- 
#in other words, set aside some testing data at the very beginning and don't use it until you are 
#evaluating accuracy at the very end

#choose a larger number (more than ten) of top predictive proteins using each selection method

#use a fuzzy intersection instead of a hard intersection to combine the sets of top predictive proteins across selection methods

library(tidyverse)
library(tidymodels)

set.seed(123)


# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

## read biomarker_raw to make split
biomarker_raw <- read_csv('data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group))

data_split <- initial_split(biomarker_raw, prop = 0.80, strata = group)
train_data <- training(data_split)
test_data  <- testing(data_split)

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

biomarker_train <- train_data %>% mutate(across(.cols = -c(group, ados), 
                                               ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  # reorder columns
  select(group, ados, everything())


train_nested <- biomarker_train %>%
  mutate(group = as.factor(group)) %>% # making sure the group is a factor
  select(-ados)%>% # removing ados
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  nest(data = c(level, group))


test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

tt_out <- train_nested %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value)

tt_out %>% head(5)


# multiple testing corrections
m <- nrow(tt_out)
hm <- log(m) + 1/(2*m) - digamma(1)

tt_corrected <- tt_out %>%
  select(data, protein, p_value) %>%
  mutate(rank = row_number()) %>%
  mutate(p_bh = p_value*m/rank,
         p_by = p_value*m*hm/rank,
         p_bonf = p_value*m)

tt_corrected %>% head(5)








# export as r binary
#save(list = 'biomarker_clean', 
#file = 'data/biomarker-clean.RData')
