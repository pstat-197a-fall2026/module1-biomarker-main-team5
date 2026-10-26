library(ggplot2)
library(tidyverse)
raw_data <- read.csv("data/biomarker-raw.csv")
numeric_data <- raw_data[sapply(raw_data, is.numeric)]


biggest_deviation <- function(x) {
  return(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

highest_deviation_protein <- sapply(numeric_data, biggest_deviation) %>% which.max() %>% names()