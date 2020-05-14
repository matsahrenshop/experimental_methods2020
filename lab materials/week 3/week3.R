############# Week 3 ##############
## Programmed by: Mats Ahrenshop ##
####### Date: 12 May 2020 #########

# COMPLIANCE #

library(sandwich)
library(lmtest)
library(kableExtra)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(plyr)
library(ebal)

data <- read.csv("kalla20.csv", stringsAsFactors = FALSE)

#---- INFORMATION ON DATA AND DESIGN ----

## Different stages of reaching subjects:
# Baseline survey
# Contact = coming to the door
# Giving initial rating conditional on coming to the door
# Canvass completion rate 

## Corresponding variables in the data set
# responded == 1 -- voters responded to baseline survey
# treat == {0, 1, 2} -- assignment to placebo, abbreviated, full intervention
# treat_full == 1 -- dummy for [treat == 2]
# canvassed == 1 -- voters came to the door
# has_first_rating == 1 -- voter responded to initial question of intervention
# tALL_factor_ -- is the pooled outcome index


#---- HOUSEKEEPING AND FUNCTIONS ----

t0.covariate.names <- c('t0_imm_better_worse', 't0_imm_police',
                        't0_imm_driverslicense', 't0_imm_daca',
                        't0_imm_citizenship', 't0_imm_deportall',
                        't0_imm_attorney', 't0_imm_prej_living',
                        't0_imm_prej_neighbor', 't0_imm_prej_speaking',
                        't0_imm_prej_workethic', 't0_imm_prej_fit',
                        't0_imm_know', 't0_social_distance_immigrant',
                        't0_therm_illegal_immigrant', 't0_therm_legal_immigrant',
                        't0_college_educ', 't0_asian', 't0_latino',
                        't0_black', 't0_white', 't0_born_in_us',
                        't0_factor_undoc_immigrant', 't0_factor_lgbt',
                        't0_factor_trump', 'vf_age',
                        'vf_voted08', 'vf_voted10', 'vf_voted12',
                        'vf_voted14', 'vf_voted16', 'vf_female',
                        'site_oc', 'site_tn')

x <- data[,c(t0.covariate.names)]
x <- as.matrix(x, dimnames = list(NULL, names(x)))

# Function to compute clustered standard errors, from Mahmood Arai.
cl <- function(fm, cluster){
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm), 2, function(x) tapply(x, cluster, sum))
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  coeftest(fm, vcovCL)
}

# Function to extract the ATE from OLS with clustered SEs.
est.ate <- function(dv, include.obs = NULL, include.covariates = TRUE, include.placebo = TRUE){
  if(is.null(include.obs)){
    include.obs <- !is.na(dv) 
  }
  include.obs <- which(include.obs & !is.na(dv))
  
  if(include.covariates & include.placebo){
    lm.obj <- lm(dv[include.obs] ~ data$treat_full[include.obs] +
                   data$treat_mod[include.obs] +
                   x[include.obs,])
    # Calculate cluster-robust standard errors.
    result <- cl(lm.obj, data$hh_id[include.obs])[2:3,]
    result <- data.frame(result)
    names(result) <- c("Effect", "SE", "t-stat", "p")
    rownames(result) <- c("Full", "Abbrev")
  }
  if(!include.covariates & include.placebo){
    lm.obj <- lm(dv[include.obs] ~ data$treat_full[include.obs] +
                   data$treat_mod[include.obs])
    # Calculate cluster-robust standard errors.
    result <- cl(lm.obj, data$hh_id[include.obs])[2:3,]
    result <- data.frame(result)
    names(result) <- c("Effect", "SE", "t-stat", "p")
    rownames(result) <- c("Full", "Abbrev")
  }
  if(include.covariates & !include.placebo){
    lm.obj <- lm(dv[include.obs] ~ data$treat_full[include.obs] + 
                   x[include.obs,])
    # Calculate cluster-robust standard errors.
    result <- cl(lm.obj, data$hh_id[include.obs])[2,]
    result <- data.frame(t(result))
    names(result) <- c("Effect", "SE", "t-stat", "p")
    rownames(result) <- c("Full vs Abbrev")
  }
  if(!include.covariates & !include.placebo){
    lm.obj <- lm(dv[include.obs] ~ data$treat_full[include.obs])
    # Calculate cluster-robust standard errors.
    result <- cl(lm.obj, data$hh_id[include.obs])[2,]
    result <- data.frame(t(result))
    names(result) <- c("Effect", "SE", "t-stat", "p")
    rownames(result) <- c("Full vs Abbrev")
  }
  
  return(result)
}

#---- ITT ----

make.results.table.ate <- function(dv, caption) {
  t1.ate.nocovars <- est.ate(data[,paste0("t1_",dv)], data$t1_respondent==1, include.covariates = FALSE)
  t1.ate.covars <- est.ate(data[,paste0("t1_",dv)], data$t1_respondent==1, include.covariates = TRUE)
  t1.ate.treatVS.nocovars <- est.ate(data[,paste0("t1_",dv)], 
                                     data$t1_respondent==1 & data$treat_string != "Placebo",
                                     include.covariates = FALSE,
                                     include.placebo = FALSE)
  t1.ate.treatVS.covars <- est.ate(data[,paste0("t1_",dv)], 
                                   data$t1_respondent==1 & data$treat_string != "Placebo",
                                   include.covariates = TRUE,
                                   include.placebo = FALSE)
  t1.covars <- rbind(t1.ate.covars, t1.ate.treatVS.covars)
  t1.nocovars <- rbind(t1.ate.nocovars, t1.ate.treatVS.nocovars)
  t1 <- cbind(t1.covars, t1.nocovars)
  
  t2.ate.nocovars <- est.ate(data[,paste0("t2_",dv)], data$t2_respondent==1, include.covariates = FALSE)
  t2.ate.covars <- est.ate(data[,paste0("t2_",dv)], data$t2_respondent==1, include.covariates = TRUE)
  t2.ate.treatVS.nocovars <- est.ate(data[,paste0("t2_",dv)], 
                                     data$t2_respondent==1 & data$treat_string != "Placebo",
                                     include.covariates = FALSE,
                                     include.placebo = FALSE)
  t2.ate.treatVS.covars <- est.ate(data[,paste0("t2_",dv)], 
                                   data$t2_respondent==1 & data$treat_string != "Placebo",
                                   include.covariates = TRUE,
                                   include.placebo = FALSE)
  t2.covars <- rbind(t2.ate.covars, t2.ate.treatVS.covars)
  t2.nocovars <- rbind(t2.ate.nocovars, t2.ate.treatVS.nocovars)
  t2 <- cbind(t2.covars, t2.nocovars)
  
  t3.ate.nocovars <- est.ate(data[,paste0("t3_",dv)], data$t3_respondent==1, include.covariates = FALSE)
  t3.ate.covars <- est.ate(data[,paste0("t3_",dv)], data$t3_respondent==1, include.covariates = TRUE)
  t3.ate.treatVS.nocovars <- est.ate(data[,paste0("t3_",dv)], 
                                     data$t3_respondent==1 & data$treat_string != "Placebo",
                                     include.covariates = FALSE,
                                     include.placebo = FALSE)
  t3.ate.treatVS.covars <- est.ate(data[,paste0("t3_",dv)], 
                                   data$t3_respondent==1 & data$treat_string != "Placebo",
                                   include.covariates = TRUE,
                                   include.placebo = FALSE)
  t3.covars <- rbind(t3.ate.covars, t3.ate.treatVS.covars)
  t3.nocovars <- rbind(t3.ate.nocovars, t3.ate.treatVS.nocovars)
  t3 <- cbind(t3.covars, t3.nocovars)
  
  tALL.ate.nocovars <- est.ate(data[,paste0("tALL_",dv)], include.covariates = FALSE)
  tALL.ate.covars <- est.ate(data[,paste0("tALL_",dv)], include.covariates = TRUE)
  tALL.ate.treatVS.nocovars <- est.ate(data[,paste0("tALL_",dv)],
                                       data$treat_string != "Placebo",
                                       include.covariates = FALSE,
                                       include.placebo = FALSE)
  tALL.ate.treatVS.covars <- est.ate(data[,paste0("tALL_",dv)],
                                     data$treat_string != "Placebo",
                                     include.covariates = TRUE,
                                     include.placebo = FALSE)
  tALL.covars <- rbind(tALL.ate.covars, tALL.ate.treatVS.covars)
  tALL.nocovars <- rbind(tALL.ate.nocovars, tALL.ate.treatVS.nocovars)
  tALL <- cbind(tALL.covars, tALL.nocovars)
  
  overall <- round(rbind(t1, t2, t3, tALL), 4)
  
  names(overall) <- rep(c("Effect", "SE", "t.stat", "p"), 2)
  overall <- as.matrix(overall)
  rownames(overall) <- rep(c("Full vs. Placebo", "Abbrev. vs. Placebo", "Full vs. Abbrev."), 4)
  return(kable(overall, digits=3, caption = caption) %>%
           add_header_above(c(" " = 1, "With Covariates" = 4, "Without Covariates" = 4)) %>%
           kableExtra::group_rows("1 Week", 1, 3) %>%
           kableExtra::group_rows("1 Month", 4, 6) %>%
           kableExtra::group_rows("3-6 Months", 7, 9) %>%
           kableExtra::group_rows("Pooled", 10, 12) %>%
           kable_styling(latex_options = c("striped", "HOLD_position")))
}

make.results.table.ate("factor_overall", "ATE effects on overall index")


#---- CACE ----

balance.vars <- c('vf_age', 'vf_female', 
                  'vf_latino', 't0_therm_legal_immigrant',
                  't0_therm_illegal_immigrant', 't0_factor_undoc_immigrant')
balance.vars.names <- c("Age", "Female", "Latino", 
                        "Legal Immigrant Feeling Thermometer t0", 
                        "Illegal Immigrant Feeling Thermometer t0", 
                        "Baseline Factor of Support")

make.balance.table <- function(subset, varlist, names, caption){
  file <- subset(data, subset)
  balance <- matrix(ncol=3, nrow=length(varlist)+1)
  for(i in 1:length(varlist)){
    balance[i,1] <- mean(file[file$treat == 0,balance.vars[i]])
    balance[i,2] <- mean(file[file$treat == 1,balance.vars[i]])
    balance[i,3] <- mean(file[file$treat == 2,balance.vars[i]])
  }
  balance[length(varlist)+1,1] <- nrow(file[file$treat==0,])
  balance[length(varlist)+1,2] <- nrow(file[file$treat==1,])
  balance[length(varlist)+1,3] <- nrow(file[file$treat==2,])
  balance <- data.frame(round(balance, digits=2))
  rownames(balance) <- c(names, "N")
  anova.test.vector <- matrix(ncol=1, nrow=nrow(balance))
  for(i in 1:length(varlist)){
    anova.test.vector[i,1] <- round(summary(aov(file[,varlist[i]] ~ 
                                                  as.factor(treat), data = file))[[1]][["Pr(>F)"]][[1]], 2)
  }
  anova.test.vector[nrow(balance),1] <- "-"
  balance <- cbind(balance,anova.test.vector)
  colnames(balance) <- c("Placebo", "Abbrev Intervention", "Full Intervention", "p-value")
  return(kable(balance, caption = caption)  %>%
           kable_styling(latex_options = c("striped", "HOLD_position")))
} 

make.balance.table(data$canvassed == 1, balance.vars, balance.vars.names,
                   "Covariate Balance among Compliers.")

## contact rates
  
# Placebo:
round(mean(subset(data, data$treat == 0)$canvassed, na.rm = TRUE), 2)

# Abbreviated Intervention:
round(mean(subset(data, data$treat == 1)$canvassed, na.rm = TRUE), 2)

# Full Intervention:
round(mean(subset(data, data$treat == 2)$canvassed, na.rm = TRUE), 2)


## compute share of compliers manually 

# Abbreviated Intervention:
round(mean(subset(data, data$treat == 1 & data$canvassed ==1)$has_first_rating, na.rm = TRUE), 2)

# Full Intervention:
round(mean(subset(data, data$treat == 2 & data$canvassed ==1)$has_first_rating, na.rm = TRUE), 2)

## compute share of compliers in regression framework

data$responded.to.any <- data$t1_respondent | data$t2_respondent | data$t3_respondent
data$responded.to.any[is.na(data$responded.to.any)] <- FALSE

itt_d_full <- lm(has_first_rating ~ treat_full + treat_mod,
                 data[data$responded.to.any,])$coefficients[2]

itt_d_mod <- lm(has_first_rating ~ treat_full + treat_mod,
                data[data$responded.to.any,])$coefficients[3]

## CACE
round(est.ate(data$tALL_factor_overall, include.covariates = TRUE)[1,1] / itt_d_full, 3)
round(est.ate(data$tALL_factor_overall, include.covariates = TRUE)[2,1] / itt_d_mod, 3)
