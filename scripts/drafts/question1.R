#What do you imagine is the reason for log-transforming the protein levels in biomarker-raw.csv? (Hint: look at the distribution of the raw values for a sample of proteins.)

library(ggplot2)
library(tidyverse)
raw_data <- read.csv("data/biomarker-raw.csv")

numeric_data <- raw_data[sapply(raw_data, is.numeric)]

biggest_deviation <- function(x) {
  return(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

highest_deviation_protein <- sapply(numeric_data, biggest_deviation) %>% which.max() %>% names()

# plotting dist
ggplot(raw_data, aes(x = .data[[highest_deviation_protein]])) +
  geom_histogram(bins = 30) + 
  labs(title = paste("Dist of", highest_deviation_protein),
       x = "Protein Level",
       y = "Frequency") +
  theme_minimal()

# plotting same dist but log-transformed
ggplot(raw_data, aes(x = log(.data[[highest_deviation_protein]]))) +
  geom_histogram(bins = 30) + 
  labs(title = paste("Dist of", highest_deviation_protein), 
       x = "Protein Level", 
       y = "Frequency") +
  theme_minimal()


# The thing about the proteins is that they can have a super high range of values, with basically no limit to how high they can go. This makes the dist of the raw values right skewed, with a super long tail, when really the majority of the values are concentrated in a much smaller range. This is super clearly shown in our normal vs log-transformed plot of Protein 4.1. (Although that protein was selected for having a high deviation). 

#Log-transforming the values helps to compress the scale, making the distribution more symmetric and easier to analyze. It also helps to reduce the impact of extreme outliers. With big scales like this, the log transformation can make it easier to visualize and interpret the data, but also make models fit better should they assume normality.