# Heartsteps data set
# WCLS analysis using xgeepack.R
# (require to source("xgeepack.R"))

# For the MRT paper submitted to Psychology Methods

# Tianchen Qian
# 2020.11.25

rm(list = ls())

library(tidyverse)
library(xtable)

DIAGNOSTIC_CODE <- FALSE # run diagnostic code (plots/view data frames) throughout the code?


# Preprocessing -----------------------------------------------------------

##### load data

# loading data from mounted mBox "HeartSteps" folder
sys.var <- switch(Sys.info()["sysname"],
                  "Windows" = list(locale = "English",
                                   mbox = "Z:/HeartSteps/"),
                  "Darwin" = list(locale = "en_US",
                                  mbox = "/Volumes/dav/HeartSteps/"),
                  "Linux" = list(locale = "en_US.UTF-8",
                                 mbox = "~/mbox/HeartSteps/"))
sys.var$mbox.data <- paste0(sys.var$mbox, "Data/")
Sys.setlocale("LC_ALL", sys.var$locale)
Sys.setenv(TZ = "GMT")
load(paste0(sys.var$mbox.data, "analysis.RData"))


# Ordering the list in terms of user and then decision index
suggest <- suggest[order(suggest$user,suggest$decision.index),]

suggest$at.home.or.work <- (suggest$location.category %in% c("home", "work"))


##### Exclude 734 decision points in suggest data.frame

# This code chunk is copied from Brook's "main effects replication.R":
# git_HeartstepsV1Code/brook/maineff_replaction.Rmd

# This code by Brook code removed the following decision points from the analysis:
# (a) travel days (ntravel = 390 decision points)
# (b) study.day.nogap past 41 (npost41 = 340 decision points)
# (c) unavail_sent_slots (user deemed unavailable but a message was sent: 3 decision points)
# (d) no_message_tag (1 decision point)

unavail_sent_slots <-
    filter(suggest, !avail & !is.na(notification.message)) %>%
    dplyr::select(user, study.day.nogap, slot, avail, notification.message) %>%
    mutate('Message sent' = !is.na(notification.message)) %>% dplyr::select(-notification.message) %>%
    rename('User' = user, 'Non-travel study day'= study.day.nogap,
           'Decision slot' = slot, 'Available' = avail)
no_message_tag <- 
    filter(suggest, send & !travel & study.day.nogap <= 42 & is.na(send.sedentary)) %>%
    dplyr::select(user, study.day.nogap, slot, 
           notification.message) %>%
    rename('User' = user,
           'Non-travel study day' = study.day.nogap,
           'Decision slot' = slot,
           'Sent message' = notification.message)

ntravel <- sum(suggest$travel)
npost41 <- with(suggest, sum(study.day.nogap > 41, na.rm=T))
nexclude <- 
    ntravel + npost41 + nrow(unavail_sent_slots) + nrow(no_message_tag)

suggest.included <- 
    suggest %>%
    filter(!travel & study.day.nogap <= 41) %>%
    anti_join(mutate(no_message_tag, no_message_tag = T),
              by=c('user'='User','study.day.nogap'='Non-travel study day',
                   'slot'='Decision slot')) %>%
    anti_join(mutate(unavail_sent_slots, unavail_sent_slots = T),
              by=c('user'='User','study.day.nogap'='Non-travel study day',
                   'slot'='Decision slot'))
suggest.analysis <-
    suggest.included %>%
    arrange(user, study.day.nogap, decision.index.nogap)

navail <- sum(suggest.included$avail)


if (DIAGNOSTIC_CODE) {
    plot(suggest.included$study.day, cex = 0.1)
    plot(suggest.included$study.day.nogap, cex = 0.1)
}


##### Add indicator of planning on the previous day to the suggest data set

# summary(suggest.included$study.date)
# summary(daily$study.date)

for (i in 1:nrow(suggest.included)) {
    if (i %% 100 == 0) {
        cat(i, "")
    }
    row_index <- which((daily$study.date == suggest.included$study.date[i] - 1) & (daily$user == suggest.included$user[i]))
    if (length(row_index) == 1) {
        if (daily$planning[row_index] %in% c("structured", "unstructured")) {
            suggest.included$planning_previousday[i] <- 1
        } else {
            suggest.included$planning_previousday[i] <- 0
        }
    } else if (length(row_index) == 0) { # planning indicator of previous day not found, likely because it is the first day of a user
        suggest.included$planning_previousday[i] <- 0
    } else {
        stop("Multiple instances found in daily.")
    }
}



##### create variables for suggest.included

suggest.included <- as_tibble(suggest.included)

if (DIAGNOSTIC_CODE) {
    plot(suggest.included$decision.index, cex = 0.05)
    plot(suggest.included$decision.index.nogap, cex = 0.05)
}

# create lagged 30min step
suggest.included$jbsteps30.log.lag1 <- c(0, suggest.included$jbsteps30.log[1 : (nrow(suggest.included) - 1)])
suggest.included$jbsteps30.log.lag1[suggest.included$decision.index.nogap == 0] <- 0

# create variables to pass into regression
suggest.included$"(Intercept)" <- 1
suggest.included$"I(send.active - 0.3)" <- suggest.included$send.active - 0.3
suggest.included$"I(send.sedentary - 0.3)" <- suggest.included$send.sedentary - 0.3
suggest.included$"I(send - 0.6)" <- suggest.included$send - 0.6

##### check number of 0's in 30-min step count

mean(suggest.included$jbsteps30 == 0, na.rm = TRUE) # 30% are 0's



# Fit WCLS ----------------------------------------------------------------


library(geepack)
source("xgeepack.R")

# The following "estimate" function is copied from xgeepack.R,
# and I added the degrees of freedom in the output table.
estimate <- function(x, combos = NULL, omnibus = FALSE, null = 0,
                     small = TRUE, conf.int = 0.95, normal = FALSE, ...) {
    if (is.null(combos)) {
        combos <- diag(length(coef(x)))
        rownames(combos) <- names(coef(x))
        omnibus <- FALSE
    }
    est <- combos %*% coef(x)
    if (nrow(est) != length(null)) null <- rep(null[1], nrow(est))
    ## apply Mancl and DeRouen's (2001) small sample correction
    if (is.logical(small)) small <- small * 50
    n <- cluster.number(x, overall = FALSE)
    d1 <- if (omnibus) nrow(combos)
    else apply(combos != 0, 1, sum)
    d2 <- n - length(coef(x))
    ## apply Hotelling's T-squared test, following Liao et al. (2016)
    if (n <= small & !normal) {
        type <- "Hotelling"
        adj <- d1 * (d1 + d2 - 1) / d2
        qfun <- function(p) mapply(qf, p = p, df1 = d1, df2 = d2) / adj
        pfun <- function(q) 1 - mapply(pf, q = q * adj, df1 = d1, df2 = d2)
    }
    else {
        type <- "Wald"
        qfun <- if (normal) function(p) qnorm((1 + p) / 2)
        else function(p) mapply(qf, p = p, df1 = d1, df2 = d2)
        pfun <- if (normal) function(q) 1 - mapply(pchisq, q = q, df = d1)
        else function(q) 1 - mapply(pf, q = q, df1 = d1, df2 = d2)
    }
    var.est <- combos %*% vcov(x, small = small, ...) %*% t(combos)
    se.est <- sqrt(diag(var.est))
    crit <- sqrt(qfun(conf.int))
    lcl <- est - se.est * crit
    ucl <- est + se.est * crit
    stat <- if (omnibus) rep(t(est - null) %*% solve(var.est) %*% (est - null), d1)
    else (est - null)^2 / diag(var.est)
    pvalue <- pfun(stat)
    out <- cbind(est, lcl, ucl, se.est, stat, d1, d2, pvalue)
    rownames(out) <- rownames(combos)
    colnames(out) <- c("Estimate",
                       paste0(round(conf.int * 100), "% ", c("LCL", "UCL")),
                       "SE", type, "df1", "df2", "p-value")
    class(out) <- c("estimate", "matrix")
    out
}


##### Model 1. marginal effect #####

# Question 1: What is the effect of delivering activity suggestions on individuals’ subsequent 30-minute step counts?

fit_model1 <- geese.glm(x = as.matrix(suggest.included[, c("(Intercept)", "jbsteps30pre.log", "I(send - 0.6)")]),
                 y = suggest.included$jbsteps30.log, w = suggest.included$avail, id = as.factor(suggest.included$user),
                 family = gaussian(), corstr = "independence")
estimate(fit_model1)

output <- estimate(fit_model1)
write.csv(round(output, 3), file = "output.csv")

# > output
#                   Estimate   95% LCL   95% UCL        SE Hotelling       df1 df2 p-value    
# (Intercept)        1.78314   1.53732   2.02896   0.12096 217.31012   1.00000  34  <1e-04 ***
# jbsteps30pre.log   0.41360   0.35116   0.47604   0.03072 181.21118   1.00000  34  <1e-04 ***
# I(send - 0.6)      0.13120  -0.00576   0.26815   0.06739   3.78979   1.00000  34  0.0599 .  




##### Model 2. effect change over time #####

# Question 2: How does the effect of activity activity suggestions change with each additional day in the study?

xmat <- suggest.included %>%
    transmute("(Intercept)" = .$"(Intercept)",
              "jbsteps30pre.log" = .$"jbsteps30pre.log",
              "study.day.nogap" = .$"study.day.nogap",
              "I(send - 0.6)" = .$"I(send - 0.6)",
              "I(send - 0.6):study.day.nogap" = .$"I(send - 0.6)" * .$"study.day.nogap")

fit_model2 <- geese.glm(x = as.matrix(xmat),
                        y = suggest.included$jbsteps30.log, w = suggest.included$avail, id = as.factor(suggest.included$user),
                        family = gaussian(), corstr = "independence")
estimate(fit_model2)

output <- estimate(fit_model2)
write.csv(round(output, 3), file = "output.csv")


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


p <- ggplot(df_tx) + 
    geom_line(aes(x = study.day.nogap, y = treatment_effect), color = "blue") +
    geom_ribbon(aes(ymin = left_ci, ymax = right_ci, x = study.day.nogap), alpha = 0.3, fill = "blue") +
    geom_hline(yintercept = 0, color = "black") +
    xlab(label = "day in the study") +
    ylab(label = "treatment effect\n(on log step count)") +
    ggtitle(paste0("Effect of activity suggestion over time")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

dir.create("plot", showWarnings = FALSE, recursive = "TRUE")
ggsave(p, filename = paste0("plot/model-2.png"), width = 5, height = 3)

# > output
#                                Estimate   95% LCL   95% UCL        SE Hotelling       df1 df2 p-value    
# (Intercept)                     2.00283   1.76518   2.24048   0.11667 294.69069   1.00000  32 < 1e-04 ***
# jbsteps30pre.log                0.41200   0.35104   0.47295   0.02992 189.56395   1.00000  32 < 1e-04 ***
# study.day.nogap                -0.01058  -0.02013  -0.00103   0.00469   5.09091   1.00000  32 0.03102 *  
# I(send - 0.6)                   0.50744   0.20086   0.81402   0.15051  11.36666   1.00000  32 0.00197 ** 
# I(send - 0.6):study.day.nogap  -0.01848  -0.03090  -0.00607   0.00610   9.19229   1.00000  32 0.00479 ** 




##### Model 3. effect moderated by location (home/work vs. other) #####

# Question 3: How does the effect of each type of activity suggestions depend on the location of the individual?

suggest.included$location.other <- !(suggest.included$location.category %in% c("home", "work"))

suggest.included$location.homework <- (suggest.included$location.category %in% c("home", "work"))

xmat <- suggest.included %>%
    transmute("(Intercept)" = .$"(Intercept)",
              "jbsteps30pre.log" = .$"jbsteps30pre.log",
              "location.homework" = .$"location.homework",
              "I(send.active - 0.3)" = .$"I(send.active - 0.3)",
              "I(send.active - 0.3):location.homework" = .$"I(send.active - 0.3)" * .$"location.homework",
              "I(send.sedentary - 0.3)" = .$"I(send.sedentary - 0.3)",
              "I(send.sedentary - 0.3):location.homework" = .$"I(send.sedentary - 0.3)" * .$"location.homework")

fit_model3 <- geese.glm(x = as.matrix(xmat),
                        y = suggest.included$jbsteps30.log, w = suggest.included$avail, id = as.factor(suggest.included$user),
                        family = gaussian(), corstr = "independence")
estimate(fit_model3)

#                                            Estimate   95% LCL   95% UCL        SE Hotelling       df1 df2 p-value    
# (Intercept)                                1.71e+00  1.46e+00  1.97e+00  1.24e-01  1.91e+02  1.00e+00  30  <1e-04 ***
# jbsteps30pre.log                           4.14e-01  3.51e-01  4.77e-01  3.07e-02  1.82e+02  1.00e+00  30  <1e-04 ***
# location.homework                          1.43e-01 -8.26e-02  3.68e-01  1.10e-01  1.67e+00  1.00e+00  30  0.2055    
# I(send.active - 0.3)                       5.02e-02 -1.67e-01  2.67e-01  1.06e-01  2.23e-01  1.00e+00  30  0.6398    
# I(send.active - 0.3):location.homework     3.77e-01  3.99e-04  7.53e-01  1.84e-01  4.18e+00  1.00e+00  30  0.0498 *  
# I(send.sedentary - 0.3)                    9.22e-02 -1.66e-01  3.51e-01  1.27e-01  5.30e-01  1.00e+00  30  0.4721    
# I(send.sedentary - 0.3):location.homework -1.42e-01 -5.40e-01  2.56e-01  1.95e-01  5.31e-01  1.00e+00  30  0.4720   



##### Model 4. effect moderated by planning on previous day #####

# Question 4: How does the effect of the activity suggestion depend on whether the individual received a planning prompt on the previous day?

xmat <- suggest.included %>%
    transmute("(Intercept)" = .$"(Intercept)",
              "jbsteps30pre.log" = .$"jbsteps30pre.log",
              "planning_previousday" = .$"planning_previousday",
              "I(send - 0.6)" = .$"I(send - 0.6)",
              "I(send - 0.6):planning_previousday" = .$"I(send - 0.6)" * .$"planning_previousday")

fit_model4 <- geese.glm(x = as.matrix(xmat),
                        y = suggest.included$jbsteps30.log, w = suggest.included$avail, id = as.factor(suggest.included$user),
                        family = gaussian(), corstr = "independence")
estimate(fit_model4)

output <- estimate(fit_model4)
write.csv(round(output, 3), file = "output.csv")

# > output
#
#                                    Estimate  95% LCL  95% UCL       SE Hotelling      df1 df2 p-value    
# (Intercept)                          1.7642   1.5109   2.0174   0.1243  201.3204   1.0000  32  <1e-04 ***
# jbsteps30pre.log                     0.4137   0.3509   0.4764   0.0308  180.5027   1.0000  32  <1e-04 ***
# planning_previousday                 0.0499  -0.1055   0.2053   0.0763    0.4273   1.0000  32   0.518    
# I(send - 0.6)                        0.1133  -0.0348   0.2614   0.0727    2.4290   1.0000  32   0.129    
# I(send - 0.6):planning_previousday   0.0460  -0.2275   0.3196   0.1343    0.1175   1.0000  32   0.734    
