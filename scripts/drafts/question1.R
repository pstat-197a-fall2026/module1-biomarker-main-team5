#What do you imagine is the reason for log-transforming the protein levels in biomarker-raw.csv? (Hint: look at the distribution of the raw values for a sample of proteins.)

library(ggplot2)
raw_data <- read.csv("data/biomarker-raw.csv")

#Something about how the scales can be super massive or skewed. Also how there might be a big gap in the distribution.