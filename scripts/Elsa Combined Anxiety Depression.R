setwd("S:/0859 - International evaluation of modifiable social determinants of health/datasets_harmonization")

library(tidyverse)
library(lavaan)
library(haven)
library(foreign)
library(mice)
library(readr)
library(lavaan)
library(semTools)
library(semPlot)
library(mice)

#read datasets
elsauk_dep <- readRDS("elsauk_recoded_dep.rds")
elsauk_anx <- readRDS("elsauk_recoded_anx.rds")
elsabr_dep <- readRDS("elsabr_recoded_dep.rds")
elsabr_anx <- readRDS("elsabr_recoded_anx.rds")

#change colnames to "id"
colnames(elsauk_dep)[1] <- "id"
colnames(elsauk_anx)[1] <- "id"
colnames(elsabr_anx)[1] <- "id"

#change id for ulsauk to chr instead of dbl
elsauk_dep$id <- as.character(elsauk_dep$id)
elsauk_anx$id <- as.character(elsauk_anx$id)

#check str
str(elsauk_dep)

#add cohort column (1=elsauk; 2=elsabr)
elsauk_dep$cohort <- factor(1)
elsauk_anx$cohort <- factor(1)
elsabr_anx$cohort <- factor(2)

#remove W1 of elsabr_anx to match elsauk_dep
elsabr_anx_3w <- elsabr_anx %>% select(-contains("1"))

#remove sex for elsabr
elsabr_anx_3w <- subset(elsabr_anx_3w, select = -SEXO)

#merge
elsa_merge_dep <- bind_rows(elsauk_dep, elsabr_dep)
elsa_merge_anx <- bind_rows(elsauk_anx, elsabr_anx_3w)

table(elsauk_anx$worr3)

#MGCFA invariance test

#configural model for anx
anx_configural_model <- ' 
  # Configural model: No constraints on loadings
  anx2 =~ worr2 + ganx2 + irri2 + rest2
  anx3 =~ worr3 + ganx3 + irri3 + rest3
  anx4 =~ worr4 + ganx4 + irri4 + rest4


  # Allow covariance between latent factors throughout time 
  anx2 ~~ anx3 + anx4  
  anx3 ~~ anx4 

  # Allow covariance between errors of the same variables throughout time
  worr2 ~~ worr3 + worr4
  ganx2 ~~ ganx3 + ganx4
  irri2 ~~ irri3 + irri4
  rest2 ~~ rest3 + rest4
'

anx_fit_configural <- cfa(anx_configural_model, data = elsa_merge_anx,
                          group = "cohort",
                          ordered = c("worr2", "ganx2", "irri2", "rest2",
                                      "worr3", "ganx3", "irri3", "rest3",
                                      "worr4", "ganx4", "irri4", "rest4"), 
                          parameterization = "theta", estimator = "wlsmv",
                          missing = "pairwise")

summary(anx_fit_configural, fit.measures = TRUE, standardized = TRUE)

fitMeasures(anx_fit_configural, c("cfi", "tli", "rmsea", "srmr"))


#metric model for anx
anx_metric_model <- '
  anx2 =~ W1Loading*worr2 + W2Loading*ganx2 + W3Loading*irri2 + W4Loading*rest2
  anx3 =~ W1Loading*worr3 + W2Loading*ganx3 + W3Loading*irri3 + W4Loading*rest3
  anx4 =~ W1Loading*worr4 + W2Loading*ganx4 + W3Loading*irri4 + W4Loading*rest4

  # Allow covariance between latent factors throughout time 
  anx2 ~~ anx3 + anx4  
  anx3 ~~ anx4

  # Allow covariance between errors of the same variables throughout time
  worr2 ~~ worr3 + worr4
  ganx2 ~~ ganx3 + ganx4
  irri2 ~~ irri3 + irri4
  rest2 ~~ rest3 + rest4
'

anx_fit_metric <- cfa(anx_metric_model, data = elsa_merge_anx,
                      group = "cohort",
                      group.equal = "loadings",
                      ordered = c("worr2", "ganx2", "irri2", "rest2",
                                  "worr3", "ganx3", "irri3", "rest3",
                                  "worr4", "ganx4", "irri4", "rest4"), 
                      parameterization = "theta", estimator = "wlsmv",
                      missing = "pairwise",
                      std.lv = TRUE)

summary(anx_fit_metric, fit.measures = TRUE, standardized = TRUE)

fitMeasures(anx_fit_metric, c("cfi", "tli", "rmsea", "srmr"))

#scalar model for anx
anx_scalar_model <- '
  anx2 =~ W1Loading*worr2 + W2Loading*ganx2 + W3Loading*irri2 + W4Loading*rest2
  anx3 =~ W1Loading*worr3 + W2Loading*ganx3 + W3Loading*irri3 + W4Loading*rest3
  anx4 =~ W1Loading*worr4 + W2Loading*ganx4 + W3Loading*irri4 + W4Loading*rest4

  # Allow covariance between latent factors throughout time 
  anx2 ~~ anx3 + anx4  
  anx3 ~~ anx4

  # Allow covariance between errors of the same variables throughout time
  worr2 ~~ worr3 + worr4
  ganx2 ~~ ganx3 + ganx4
  irri2 ~~ irri3 + irri4
  rest2 ~~ rest3 + rest4
  
  # Constraints intercepts to be equal across time points
  worr2 | t_worr*t1
  worr3 | t_worr*t1
  worr4 | t_worr*t1

  ganx2 | t_ganx*t1
  ganx3 | t_ganx*t1
  ganx4 | t_ganx*t1

  irri2 | t_irri*t1
  irri3 | t_irri*t1
  irri4 | t_irri*t1

  rest2 | t_rest*t1
  rest3 | t_rest*t1
  rest4 | t_rest*t1
'

anx_fit_scalar <- cfa(anx_scalar_model, data = elsa_merge_anx,
                      group = "cohort",
                      group.equal = c("loadings", "thresholds"),
                      ordered = c("worr2", "ganx2", "irri2", "rest2",
                                  "worr3", "ganx3", "irri3", "rest3",
                                  "worr4", "ganx4", "irri4", "rest4"), 
                      parameterization = "theta", estimator = "wlsmv",
                      missing = "pairwise",
                      check.gradient = FALSE)#,
                      #std.lv = TRUE)

summary(anx_fit_scalar, fit.measures = TRUE, standardized = TRUE)

fitMeasures(anx_fit_scalar, c("cfi", "tli", "rmsea", "srmr"))

#extract factor scores (, transform = TRUE)
factor_scores_anx <- lavPredict(anx_fit_scalar)
#factor_scores <- lavPredict(fit_scalar)
head(factor_scores_anx)

#try multiple imputation [it didn't help]
#imputed_anx <- mice(elsa_merge_anx, method = "cart", m=5, seed=123)
#imputed_anx_complete <- complete(imputed_anx, action=1)

###measurement invariance for depression###
dep_config_model <- '
dep1 =~ mood1 + anhe1 + guil1 + hope1 + apat1
dep2 =~ mood2 + anhe2 + guil2 + hope2 + apat2 
dep3 =~ mood3 + anhe3 + guil3 + hope3 + apat3
dep4 =~ mood4 + anhe4 + guil4 + hope4 + apat4

dep1 ~~ dep1 + dep2 + dep3 + dep4
dep2 ~~ dep2 + dep3 + dep4
dep3 ~~ dep3 + dep4
dep4 ~~ dep4

mood1 ~~ mood2 + mood3 + mood4 
mood2 ~~ mood3 + mood4
mood3 ~~ mood4

anhe1 ~~ anhe2 + anhe3 + anhe4
anhe2 ~~ anhe3 + anhe4
anhe3 ~~ anhe4

guil1 ~~ guil2 + guil3 + guil4
guil2 ~~ guil3 + guil4
guil3 ~~ guil4

hope1 ~~ hope2 + hope3 + hope4
hope2 ~~ hope3 + hope4
hope3 ~~ hope4

apat1 ~~ apat2 + apat3 + apat4
apat2 ~~ apat3 + apat4
apat3 ~~ apat4
'

# Fit
fit_configural_dep <- cfa(dep_config_model, data = elsa_merge_dep,
                      group = "cohort", #group by cohort (1=elsauk; 2=elsabr)
                      ordered = c("mood1", "anhe1", "guil1", "hope1", "apat1",
                                  "mood2", "anhe2", "guil2", "hope2", "apat2",
                                  "mood3", "anhe3", "guil3", "hope3", "apat3",
                                  "mood4", "anhe4", "guil4", "hope4", "apat4"),
                      parameterization = "theta", estimator = "WLSMV",
                      missing = "pairwise")

summary(fit_configural_dep, fit.measures = TRUE, standardized = TRUE)


fitMeasures(fit_configural_dep, c("cfi", "tli", "rmsea", "srmr"))


#metric invariance
dep_metric_model <- '
  # Constraints factorial loadings to be equal across time points
dep1 =~ V1Loading*mood1 + V2Loading*anhe1 + V3Loading*guil1 + V4Loading*hope1 + V5Loading*apat1

dep2 =~ V1Loading*mood2 + V2Loading*anhe2 + V3Loading*guil2 + V4Loading*hope2 + V5Loading*apat2

dep3 =~ V1Loading*mood3 + V2Loading*anhe3 + V3Loading*guil3 + V4Loading*hope3 + V5Loading*apat3

dep4 =~ V1Loading*mood4 + V2Loading*anhe4 + V3Loading*guil4 + V4Loading*hope4 + V5Loading*apat4


  # Allow covariance between latent factors throughtout time 
  dep1 ~~ dep2 + dep3 + dep4
  dep2 ~~ dep3 + dep4
  dep3 ~~ dep4
 

  # Allow covariance between errors of the same variables throughout time
  mood1 ~~ mood2 + mood3 + mood4
  anhe1 ~~ anhe2 + anhe3 + anhe4
  guil1 ~~ guil2 + guil3 + guil4
  hope1 ~~ hope2 + hope3 + hope4
  apat1 ~~ apat2 + apat3 + apat4
'

fit_metric_dep <- cfa(dep_metric_model, data = elsa_merge_dep,
                  group = "cohort",
                  group.equal = "loadings",
                  ordered = c("mood1", "anhe1", "guil1", "hope1", "apat1",
                              "mood2", "anhe2", "guil2", "hope2", "apat2",
                              "mood3", "anhe3", "guil3", "hope3", "apat3",
                              "mood4", "anhe4", "guil4", "hope4", "apat4"),
                  parameterization = "theta", estimator = "WLSMV",
                  missing = "pairwise")

summary(fit_metric_dep, fit.measures = TRUE, standardized = TRUE)

fitMeasures(fit_metric_dep, c("cfi", "tli", "rmsea", "srmr"))

#scalar invariance
dep_scalar_model <- '
# Constraints factorial loadings to be equal across time points
dep1 =~ V1Loading*mood1 + V2Loading*anhe1 + V3Loading*guil1 + V4Loading*hope1 + V5Loading*apat1

dep2 =~ V1Loading*mood2 + V2Loading*anhe2 + V3Loading*guil2 + V4Loading*hope2 + V5Loading*apat2

dep3 =~ V1Loading*mood3 + V2Loading*anhe3 + V3Loading*guil3 + V4Loading*hope3 + V5Loading*apat3

dep4 =~ V1Loading*mood4 + V2Loading*anhe4 + V3Loading*guil4 + V4Loading*hope4 + V5Loading*apat4


# Allow covariance between latent factors throughtout time 
dep1 ~~ dep2 + dep3 
dep2 ~~ dep3 
dep3 ~~ dep4

# Allow covariance between errors of the same variables throughout time
mood1 ~~ mood2 + mood3 + mood4
anhe1 ~~ anhe2 + anhe3 + anhe4
guil1 ~~ guil2 + guil3 + guil4
hope1 ~~ hope2 + hope3 + hope4
apat1 ~~ apat2 + apat3 + apat4


# Constrain thresholds to be equal across time points
mood1 | t_mood*t1
mood2 | t_mood*t1
mood3 | t_mood*t1
mood4 | t_mood*t1

anhe1 | t_anhe*t1
anhe2 | t_anhe*t1
anhe3 | t_anhe*t1
anhe4 | t_anhe*t1

guil1 | t_guil*t1
guil2 | t_guil*t1
guil3 | t_guil*t1
guil4 | t_guil*t1

hope1 | t_hope*t1
hope2 | t_hope*t1
hope3 | t_hope*t1
hope4 | t_hope*t1

apat1 | t_apat*t1
apat2 | t_apat*t1
apat3 | t_apat*t1
apat4 | t_apat*t1
'

fit_scalar_dep <- cfa(dep_scalar_model, data = elsa_merge_dep,
                  group = "cohort",
                  group.equal = c("loadings", "thresholds"),
                  ordered = c("mood1", "anhe1", "guil1", "hope1", "apat1",
                              "mood2", "anhe2", "guil2", "hope2", "apat2",
                              "mood3", "anhe3", "guil3", "hope3", "apat3",
                              "mood4", "anhe4", "guil4", "hope4", "apat4"),
                  parameterization = "theta", estimator = "WLSMV",
                  missing = "pairwise",
                  std.lv = TRUE)

summary(fit_scalar_dep, fit.measures = TRUE, standardized = TRUE)

fitMeasures(fit_scalar_dep, c("cfi", "tli", "rmsea", "srmr"))


#extract factor scores
#factor_scores <- lavPredict(fit_scalar, transform = TRUE)
factor_scores_dep <- lavPredict(fit_scalar_dep)
head(factor_scores_dep)

###visualisations for all models###
semPaths(fit_configural_dep, 
         what = "std", 
         whatLabels = "est", 
         layout = "tree2", 
         rotation = 1, 
         edge.label.cex = 0.45, 
         nCharNodes = 10, 
         sizeMan = 3.2,         # size of manifest variables
         label.cex = 1.2,     # font size in manifest variables
         edge.color = "black" # lines color
)

semPaths(fit_metric_dep, 
         what = "std", 
         whatLabels = "est", 
         layout = "tree2", 
         rotation = 1, 
         edge.label.cex = 0.45, 
         nCharNodes = 10, 
         sizeMan = 3.2,         # size of manifest variables
         label.cex = 1.2,     # font size in manifest variables
         edge.color = "black" # lines color
)

semPaths(fit_scalar_dep, 
         what = "std", 
         whatLabels = "est", 
         layout = "tree2", 
         rotation = 1, 
         edge.label.cex = 0.45, 
         nCharNodes = 10, 
         sizeMan = 3.2,         # size of manifest variables
         label.cex = 1.2,     # font size in manifest variables
         edge.color = "black" # lines color
)