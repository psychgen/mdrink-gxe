#02_example_analysis.R

# This script is an example analytic script to illustrate R Markdown reporting from external scripts

require(tidyverse)
require(knitr)

# Read in processed data from 01_data_prep.R

load("./data/processed_data.RData")

# Run an analysis

a <- lm(formula = "mpg ~ disp", data=mtcars)

a$coefficients
a$effects
b <-summary(a)

