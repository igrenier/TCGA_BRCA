# Title: AZ Interview Question 1
# Author: @igrenier
# Date: 04/28/2023

# Description: Modeling the 5-year survival rate in the database


# libraries
library(data.table)
library(ggplot2)
library(survival)

# Questions:
# 1. am i using the right age variable
# 2. Disease Stage: Remove the missing values, have a description of the
# subgroups and T,N and M rating and how we are going back a level.
# 3. In the Cox model, how to handle the groups with no events?
# 4. Check the assumptions of the regressions models
# 5. Plot the results and make some conclusion about the impact of age
#   for example, someone who is diagnosed at 50 is X more likely to survive
#   than someone diagnosed at 60 for example.

# Read the data
mydata <-  readxl::read_xlsx("C:/Work/CodeIG/Year_2_archive/clinical.xlsx")
mydata <- as.data.table(mydata)

# Reduce to needed variables for project
mydata_red <- mydata[, .(case_id, age_at_index, ethnicity, race, gender, ajcc_pathologic_stage,
                         days_to_birth, days_to_death, age_at_diagnosis,
                         days_to_last_follow_up)]

# remove duplicate patient information due to treatment info
mydata_red <- unique(mydata_red)

# did the patient live for 5 years, yes or no:
mydata_red[days_to_death == "'--", days_to_death := NA]
mydata_red[days_to_birth == "'--", days_to_birth := NA]
mydata_red[age_at_diagnosis == "'--", age_at_diagnosis := NA]
mydata_red[days_to_last_follow_up == "'--", days_to_last_follow_up := NA]
mydata_red[, days_to_death := as.numeric(days_to_death)]
mydata_red[, days_to_birth := as.numeric(days_to_birth)]
mydata_red[, age_at_diagnosis := as.numeric(age_at_diagnosis)]
mydata_red[, days_to_last_follow_up := as.numeric(days_to_last_follow_up)]
mydata_red[, survival_5 := "alive5"]
mydata_red[days_to_death <= 365 * 5, survival_5 := "dead"]
mydata_red[days_to_last_follow_up <= 365 * 5, survival_5 := "alive_under5"]
mydata_red[, survival_5 := factor(survival_5,
                                  levels = c("dead", "alive5", "alive_under5"))]

# issues with gender: unbalanced data, very few male cases for the disease
# I researched it and found that on average 1 in 100 breast cancer case is
# male therefore the 12 to 1000 ratio holds up.

# EDA: Some visual stuff
ggplot(mydata_red[survival_5 %in% c("dead", "alive5")],
       aes(x = survival_5, group = survival_5, y = age_at_index)) +
  geom_boxplot()

ggplot(mydata_red[survival_5 %in% c("dead", "alive5")],
       aes(x = survival_5, group = survival_5, y = ajcc_pathologic_stage)) +
  geom_boxplot()

mydata_red[survival_5 %in% c("dead", "alive5"), .N, by = .(gender, survival_5)]
mydata_red[survival_5 %in% c("dead", "alive5"), .N, keyby = .(race, survival_5)]

# EDA: Logistic regression, we only keep the 0 and 1:
mydata_log <- mydata_red[survival_5 %in% c(0, 1)]

model_0 <- glm(survival_5 ~ age_at_index + gender + race + ajcc_pathologic_stage,
               data = mydata_log)
summary(model_0)

model_1 <- glm(survival_5 ~ age_at_index + race + ajcc_pathologic_stage,
               data = mydata_log)
summary(model_1)

model_2 <- glm(survival_5 ~ age_at_index + ajcc_pathologic_stage,
               data = mydata_log)
summary(model_2)

# interaction not significant (delete)
model_3 <- glm(survival_5 ~ age_at_index * ajcc_pathologic_stage,
               data = mydata_log)
summary(model_3)

# we choose model 2 from AIC, BIC and from backward selection as well
AIC(model_0, model_1, model_2, model_3)
BIC(model_0, model_1, model_2, model_3)

# Cox model
mydata_cox <- copy(mydata_red)
mydata_cox[, survival_time := days_to_death]
mydata_cox[, status := 2]
mydata_cox[is.na(survival_time), status := 1]
mydata_cox[is.na(survival_time), survival_time := days_to_last_follow_up]

# clean up?
mydata_cox <- mydata_cox[survival_time > 0]

model_cox_0 <- coxph(Surv(time = survival_time, event = status) ~ age_at_index, data = mydata_cox)

summary(model_cox_0)

model_cox_1 <-
  coxph(Surv(time = survival_time, event = status) ~
          ajcc_pathologic_stage, data = mydata_cox)

summary(model_cox_1)

# this leads to issues because two stages don't have any "events" or death.
# One way to remedy is to group the stages into two groups Stage 1 and 2,
# Stage 3+.

mydata_cox[, stage := "Stage I"]
mydata_cox[ajcc_pathologic_stage %in% c("Stage II",
                                        "Stage IIA",
                                        "Stage IIB",
                                        "Stage IIC"),
           stage := "Stage II"]
mydata_cox[ajcc_pathologic_stage %in% c("Stage III",
                                        "Stage IIIA",
                                        "Stage IIIB",
                                        "Stage IIIC"),
           stage := "Stage III"]
mydata_cox[ajcc_pathologic_stage %in% c("Stage IV"),
            stage := "Stage IV"]
mydata_cox[ajcc_pathologic_stage %in% c("Stage X"),
           stage := "Stage X"]

model_cox_2 <-
  coxph(Surv(time = survival_time, event = status) ~
          age_at_index + stage, data = mydata_cox)

summary(model_cox_2)

# Could be done in the future
# 1. look at the correlation between therapy type and survival time
# 2. look into the Bayesian version of the model to allow for
# uncertainty quantification of the posterior parameters
