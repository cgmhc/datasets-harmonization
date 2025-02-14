setwd("/Users/henrylian/Desktop")

#load libraries
library(dplyr)
library(tidyverse)
library(haven)
library(ggplot2)
library(readxl)
library(writexl)
library(lavaan)
library(httpgd)

#sample data
path1 <- read_sav("/Users/henrylian/Documents/UNSW PhD/PATH Data/PATH Data/PATH_Wave1_Data/PATH_wave 1_220910.sav")
view(path1)

##data preparation
#columns should be ID and every symptom. rows should be each participant
#all outcomes should be coded so that higher values means worse symptomology 
#transform variables to start at 0 (optional)

#make all outcomes ordinal factors
ordered_vars <- c("mood", "anhe", "guil", "deat", "hope", "fati", 
                  "inso", "weig", "apat", "strs", "agit", "cogi",
                  "worr", "ganx", "irri", "panc", "impr", "hlth", 
                  "rest", "phob", "soma", "lone")

data[ordered_vars] <- lapply(data[ordered_vars], function(x) {
  if (!is.factor(x)) {  # Convert only if not already a factor
    return(ordered(x))  # Convert to ordered factor
  } else {
    return(x)  # Keep as is if already a factor
  }
})

str(data[ordered_vars])

# Read in data with missing values coded as 999
data <- read.csv("/Users/HenryLian/Desktop/harmonisation_example.csv")

# Define the list of categorical (ordered) variables
ordered_vars <- c("mood", "anhe", "guil", "deat", "hope", "fati", 
                  "inso", "weig", "apat", "strs", "agit", "cogi",
                  "worr", "ganx", "irri", "panc", "impr", "hlth", 
                  "rest", "phob", "soma", "lone")

# Specify the measurement model
model <- '
  # Depression latent variable
  dep =~ mood + anhe + guil + deat + hope + fati + inso + weig + apat + stre + agit + cogn + inde
  
  # Anxiety latent variable
  anx =~ worr + ganx + irri + pani + func + hlth + rest + phob + soma + situ + over + brkd + lone
'

### ------------------------------
### 1. Configural Invariance Model
### ------------------------------
# No equality constraints: the factor structure is the same across groups,
# but loadings and thresholds are freely estimated.
fit_configural <- cfa(model, 
                      data = data, 
                      group = "cohort", #e.g. 1 is elsa-brasil and 2 is elsa-uk
                      ordered = ordered_vars, 
                      estimator = "WLSMV")

# Display summary for the configural model
summary(fit_configural, standardized = TRUE, modindices = TRUE)


### ------------------------------
### 2. Metric Invariance Model
### ------------------------------
# Constrain factor loadings to be equal across groups.
fit_metric <- cfa(model, 
                  data = data, 
                  group = "cohort", 
                  ordered = ordered_vars, 
                  estimator = "WLSMV",
                  group.equal = "loadings")

# Display summary for the metric invariance model
summary(fit_metric, standardized = TRUE, modindices = TRUE)


### ------------------------------
### 3. Scalar Invariance Model
### ------------------------------
# Constrain both factor loadings and thresholds to be equal across groups.
fit_scalar <- cfa(model, 
                  data = data, 
                  group = "cohort", 
                  ordered = ordered_vars, 
                  estimator = "WLSMV",
                  group.equal = c("loadings", "thresholds"))

# Display summary for the scalar invariance model
summary(fit_scalar, standardized = TRUE, modindices = TRUE)



########################################
#fit indices evaluation
#CFI/TLI >.90 (ideally >.95)
#RMSEA <.08 (ideally .06)
#SRMR <.08
#Chi-square p>.05 (doesn't matter if large sample size)