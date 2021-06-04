# Analysis of synthetic data
# For Psych Methods paper, review round 2


# Tianchen Qian
# 2021.05.30


# data loading and preprocessing is copied from Tianchen's code for LMM paper:
# linear_mixed_model-continuous_outcome.R

rm(list = ls())

library(tidyverse)
library(xtable)
library(geepack)
source("xgeepack.R")
source("estimate.R")

synthetic_data <- read.csv("synthetic_data_37subject_210time.csv")

# create two variables to be used in the WCLS regression below
synthetic_data$"(Intercept)" <- 1
synthetic_data$"I(send - 0.6)" <- synthetic_data$send - 0.6



# Using WCLS to analyze the data ------------------------------------------


##### Model 1. marginal effect #####

# Question 1: What is the effect of delivering activity suggestions on individuals’ subsequent 30-minute step counts?

xmat <- synthetic_data %>%
    transmute("(Intercept)" = .$"(Intercept)",
              "jbsteps30pre.log" = .$"jbsteps30pre.log",
              "I(send - 0.6)" = .$"I(send - 0.6)")

fit_model1 <- geese.glm(x = as.matrix(xmat),
                        y = synthetic_data$jbsteps30.log, w = synthetic_data$avail, id = as.factor(synthetic_data$user),
                        family = gaussian(), corstr = "independence")
estimate(fit_model1)


#                  Estimate  95% LCL  95% UCL       SE Hotelling      df1 df2 p-value    
# (Intercept)      2.01e+00 1.92e+00 2.10e+00 4.57e-02  1.94e+03 1.00e+00  34  <1e-04 ***
# jbsteps30pre.log 3.40e-01 3.00e-01 3.80e-01 1.97e-02  2.98e+02 1.00e+00  34  <1e-04 ***
# I(send - 0.6)    1.57e-01 3.10e-02 2.84e-01 6.22e-02  6.40e+00 1.00e+00  34  0.0162 *  


xmat <- synthetic_data %>%
    transmute("(Intercept)" = .$"(Intercept)",
              "jbsteps30pre.log" = .$"jbsteps30pre.log",
              "jbsteps30.log.lag1" = .$"jbsteps30.log.lag1",
              "I(send - 0.6)" = .$"I(send - 0.6)")

fit_model1.1 <- geese.glm(x = as.matrix(xmat),
                          y = synthetic_data$jbsteps30.log,
                          w = synthetic_data$avail,
                          id = as.factor(synthetic_data$user),
                          family = gaussian(), corstr = "independence")
estimate(fit_model1.1)

#                    Estimate  95% LCL  95% UCL       SE Hotelling      df1 df2 p-value    
# (Intercept)        1.90e+00 1.80e+00 2.00e+00 4.91e-02  1.50e+03 1.00e+00  33 < 1e-04 ***
# jbsteps30pre.log   3.41e-01 3.01e-01 3.81e-01 1.98e-02  2.98e+02 1.00e+00  33 < 1e-04 ***
# jbsteps30.log.lag1 3.97e-02 1.69e-02 6.25e-02 1.12e-02  1.25e+01 1.00e+00  33 0.00121 ** 
# I(send - 0.6)      1.61e-01 3.80e-02 2.85e-01 6.07e-02  7.08e+00 1.00e+00  33 0.01197 *  


##### Model 2. effect change over time #####

# Question 2: How does the effect of activity suggestions change with each additional day in the study?

xmat <- synthetic_data %>%
    transmute("(Intercept)" = .$"(Intercept)",
              "jbsteps30pre.log" = .$"jbsteps30pre.log",
              "study.day.nogap" = .$"study.day.nogap",
              "I(send - 0.6)" = .$"I(send - 0.6)",
              "I(send - 0.6):study.day.nogap" = .$"I(send - 0.6)" * .$"study.day.nogap")

fit_model2 <- geese.glm(x = as.matrix(xmat),
                        y = synthetic_data$jbsteps30.log, w = synthetic_data$avail, id = as.factor(synthetic_data$user),
                        family = gaussian(), corstr = "independence")
estimate(fit_model2)

#                                Estimate   95% LCL   95% UCL        SE Hotelling       df1 df2 p-value    
# (Intercept)                     2.18752   2.04533   2.32971   0.06981 982.04017   1.00000  32  <1e-04 ***
# jbsteps30pre.log                0.33967   0.30003   0.37930   0.01946 304.75364   1.00000  32  <1e-04 ***
# study.day.nogap                -0.00852  -0.01398  -0.00305   0.00268  10.08379   1.00000  32  0.0033 ** 
# I(send - 0.6)                   0.64860   0.43050   0.86670   0.10707  36.69331   1.00000  32  <1e-04 ***
# I(send - 0.6):study.day.nogap  -0.02374  -0.03279  -0.01469   0.00444  28.55599   1.00000  32  <1e-04 ***


beta_index <- 4:5

beta_hat <- coef(fit_model2)[beta_index]
vcov <- fit_model2$geese$vbeta[beta_index, beta_index]

newdta <- df_tx <- data.frame(Intercept = 1, study.day.nogap = 0:41)
df_tx$treatment_effect <- as.matrix(newdta) %*% beta_hat
df_tx$tx_se <- NA
for (i in 1:nrow(df_tx)) {
    f_t <- as.numeric(newdta[i, ]) # feature
    df_tx$tx_se[i] <- sqrt(t(f_t) %*% vcov %*% f_t)
}
df_tx$left_ci <- df_tx$treatment_effect - 1.96 * df_tx$tx_se
df_tx$right_ci <- df_tx$treatment_effect + 1.96 * df_tx$tx_se

df_tx_linear <- df_tx


ggplot(df_tx) + 
    geom_line(aes(x = study.day.nogap, y = treatment_effect), color = "blue") +
    geom_ribbon(aes(ymin = left_ci, ymax = right_ci, x = study.day.nogap), alpha = 0.3, fill = "blue") +
    geom_hline(yintercept = 0, color = "black") +
    xlab(label = "day in the study") +
    ylab(label = "treatment effect\n(on log step count)") +
    ggtitle(paste0("Effect of activity suggestion over time")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))



##### Model 3. effect moderated by outcome at previous time point #####

# Question 3: How does the effect of activity suggestions depend on the step count at previous decision point?
# Note: this analysis is not included in the real data analysis in the original Psych Methods paper.
# I included it here to illustrate the analysis of effect moderation by a time-varying covariate.

xmat <- synthetic_data %>%
    transmute("(Intercept)" = .$"(Intercept)",
              "jbsteps30pre.log" = .$"jbsteps30pre.log",
              "jbsteps30.log.lag1" = .$"jbsteps30.log.lag1",
              "location.homework" = .$"location.homework",
              "I(send - 0.6)" = .$"I(send.active - 0.6)",
              "I(send - 0.6):jbsteps30.log.lag1" = .$"I(send - 0.6)" * .$"jbsteps30.log.lag1")

fit_model3 <- geese.glm(x = as.matrix(xmat),
                        y = synthetic_data$jbsteps30.log, w = synthetic_data$avail, id = as.factor(synthetic_data$user),
                        family = gaussian(), corstr = "independence")
estimate(fit_model3)

#                                   Estimate   95% LCL   95% UCL        SE Hotelling       df1 df2 p-value    
# (Intercept)                       1.85e+00  1.74e+00  1.96e+00  5.34e-02  1.20e+03  1.00e+00  32 < 1e-04 ***
# jbsteps30pre.log                  3.41e-01  3.01e-01  3.81e-01  1.97e-02  3.01e+02  1.00e+00  32 < 1e-04 ***
# jbsteps30.log.lag1                3.87e-02  1.55e-02  6.19e-02  1.14e-02  1.16e+01  1.00e+00  32 0.00183 ** 
# location.homework                 1.51e-01  3.93e-02  2.63e-01  5.48e-02  7.58e+00  1.00e+00  32 0.00964 ** 
# I(send - 0.6):jbsteps30.log.lag1  2.83e-02 -1.66e-03  5.82e-02  1.47e-02  3.70e+00  1.00e+00  32 0.06326 .  

