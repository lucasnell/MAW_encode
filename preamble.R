# Preamble for all ENCODE analyses


library('dplyr')
library('ggplot2')
library('stringr')
library('readr')
library('tidyr')
library('parallel')


options(stringsAsFactors = FALSE)

# Number of cores to use for parallel processes
numCores <- 3


