#question 2

#Temporarily remove the outlier trimming from preprocessing and
#do some exploratory analysis of the outlying values.
#Are there specific subjects (not values) that seem to be outliers?
#If so, are outliers more frequent in one group or the other?
#(Hint: consider tabulating the number of outlying values per subject.)



#adjusted pre processing code
library(tidyverse)


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

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_untrimmed <- read_csv('data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~ scale(log10(.x))[, 1])) %>%
  # reorder columns
  select(group, ados, everything())



#Q2

outlier_counts <- biomarker_untrimmed %>%
  mutate(subject = row_number()) %>%
  rowwise() %>%
  mutate(n_outliers = sum(abs(c_across(-c(group, ados, subject))) > 3,
                          na.rm = TRUE)) %>%
  ungroup() %>%
  select(subject, group, n_outliers)


#top 10 outliers
outlier_counts %>%
  arrange(desc(n_outliers)) %>%
  head(10)


#more frequent in one group than another?
outlier_summary <- outlier_counts %>%
  group_by(group) %>%
  summarise(
    n_subjects = n(),
    mean_outliers = mean(n_outliers),
    median_outliers = median(n_outliers),
    max_outliers = max(n_outliers)
  )

outlier_summary


#boxplot for visual analysis
ggplot(outlier_counts, aes(x = group, y = n_outliers)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5) +
  labs(
    title = "Number of Outlying Protein Measurements by Group",
    x = "Group",
    y = "Number of Outliers per Subject"
  ) +
  theme_minimal()


#summary
#Analyzed 154 subjects: 76 ASD and 78 TD
#TD subjects had a higher mean number of outliers (17.6) than ASD subjects (13.2)
#Both groups had the same median of 8.5 outliers
#Boxplot showed that most subjects had relatively few outliers, with a small number of extreme cases
#Outliers appear concentrated in certain subjects rather than evenly distributed
#Although TD had more outliers on average, the identical medians suggest the difference may be driven by a few extreme subjects
#Trimming may reduce the influence of extreme measurements but could also remove meaningful biological variation
#Overall, there is no clear difference in typical outlier frequency between ASD and TD groups
