#MG-CFA 

library(lavaan)
library(foreign)
library(tidyverse)
library(lavaanPlot)
library(effectsize)
library(semTools)

#set wd
setwd("/Users/pedrozuccolo/Desktop/harmonization/data")

#data
#dat <- read.spss("http://stats.idre.ucla.edu/wp-content/uploads/2018/05/SAQ.sav", to.data.frame = TRUE, use.value.labels = TRUE)

#check data
#str(dat)

#change to numeric
#dat2 <- dat %>%
#  mutate(across(where(is.factor), ~ as.numeric(.))) 

#corr structure
#round(cor(dat2[, 1:8]), 2)


#with ELSA Brasil
elsabr <- read.csv("elsabr_covid.csv")


View(elsabr)

dass_items <- elsabr %>%
  select(IDELSA, SEXO, matches("dass.*rec_c1")) %>%
  mutate(across(where(is.integer), ~ as.numeric(.))) %>%
  #keep only depression items
  select(-dass1_rec_c1, 
         -dass2_rec_c1,
         -dass4_rec_c1,
         -dass6_rec_c1,
         -dass7_rec_c1,
         -dass8_rec_c1,
         -dass9_rec_c1,
         -dass11_rec_c1,
         -dass12_rec_c1,
         -dass14_rec_c1,
         -dass15_rec_c1,
         -dass18_rec_c1,
         -dass19_rec_c1,
         -dass20_rec_c1
         )

str(dass_items)

#corr matrix
round(cor(dass_items[, 2:ncol(dass_items)], use = "pairwise.complete.obs"), 2)


#cov matrix
#round(cov(dat2[, 3:5]), 2)

round(cov(dass_items[, 2:ncol(dass_items)], use = "pairwise.complete.obs"), 2)



#running measurement model in lavaan

#marker method - setting the first loading to 1
m1a <- 'f =~ q03 + q04 +q05'
onefac3items_a <- cfa(m1a, data=dat2)
summary(onefac3items_a)


#variance standardization method
m1b <- 'f =~ NA*q03 + q04 +q05
       f~~1*f'
onefac3items_b <- cfa(m1b, data=dat2)
summary(onefac3items_b)


#or to use standardization of both predictor and outcome
m1a <- 'f =~ q03 + q04 +q05'
onefac3items_a <- cfa(m1a, data=dat2)
summary(onefac3items_a, standardized= TRUE)


colnames(dass_items)

#In Elsa
dep_model <- 'f =~ dass3_rec_c1 + dass5_rec_c1 + dass10_rec_c1 + dass13_rec_c1 + dass16_rec_c1 + dass17_rec_c1 + dass21_rec_c1'
dep_cfa <- cfa(dep_model, data=dass_items)
summary(dep_cfa, standardized= TRUE, fit.measures = TRUE)

# see indices that could be improved
modindices(dep_cfa, sort = TRUE)


dep_model_mod <- 'f =~ dass3_rec_c1 + dass5_rec_c1 + dass10_rec_c1 + dass13_rec_c1 + dass16_rec_c1 + dass17_rec_c1 + dass21_rec_c1
                  dass17_rec_c1 ~~ dass21_rec_c1'
dep_cfa_mod <- cfa(dep_model_mod, data=dass_items)
summary(dep_cfa_mod, standardized= TRUE, fit.measures = TRUE)


#CFI (confirmatory factor index) and TFI (Tucker Lewis Index) *need to approach 1
#RMSEA (root mean square error of approximation) < 0.05 close fit
interpret(dep_cfa_mod)


# Criar o gráfico CFA
lavaanPlot(model = dep_cfa_mod,
coefs = TRUE)


#configural invariance
#cut missing data from dass variables
dass_items <- dass_items %>%
  filter(complete.cases(.))

dep_model <- 'f =~ dass3_rec_c1 + dass5_rec_c1 + dass10_rec_c1 + dass13_rec_c1 + dass16_rec_c1 + dass17_rec_c1 + dass21_rec_c1'
dep_cfa_configural <- cfa(dep_model, data = dass_items, group="SEXO")
summary(dep_cfa_configural, fit.measures=T, standardized = T)

#metric invariance
dep_cfa_metric <- cfa(dep_model, data = dass_items, group="SEXO", group.equal = c("loadings"))
summary(dep_cfa_metric, fit.measures=T, standardized = T)

#scalar invariance
dep_cfa_scalar <- cfa(dep_model, data = dass_items, group="SEXO", group.equal = c("loadings", "intercepts"))
summary(dep_cfa_scalar, fit.measures=T, standardized = T)

#compare models
results <- compareFit(dep_cfa_configural, dep_cfa_metric, dep_cfa_scalar)
summary(results)

#alternative: compare all model at once
semTools::measurementInvariance(model = dep_model, data = dass_items, group = "SEXO")


### tutorial q usei: https://bookdown.org/content/5737/lavaan.html