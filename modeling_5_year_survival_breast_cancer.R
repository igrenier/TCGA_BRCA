# Title: Understanding the factors associated with patient survival
# Author: @igrenier
# Date: 04/28/2023

# Description: Modeling the 5-year survival rate and understand the association
# with age, gender, ethnicity, disease stage, year of diagnosis.

# libraries
library(data.table)
library(ggplot2)
library(survival)
library(broom)
library(ciTools)

# Questions:
# 4. Check the assumptions of the regressions models
# 5. Plot the results and make some conclusion about the impact of age
#   for example, someone who is diagnosed at 50 is X more likely to survive
#   than someone diagnosed at 60 for example.

# Read the data ----
clinical_data <- fread("data/clinical.tsv", sep ="\t", header = TRUE)
clinical_data <- as.data.table(clinical_data)

## reduce to needed variables for project
clinical_data <- clinical_data[, .(case_id, age_at_index, ethnicity, race,
                                   gender, ajcc_pathologic_stage,
                                   year_of_diagnosis, days_to_birth,
                                   days_to_death, age_at_diagnosis,
                                   days_to_last_follow_up,
                                   vital_status)]

## remove duplicate patient information due to treatment info
clinical_data <- unique(clinical_data)

## change the missing values to NAs and change variables from string to numeric
clinical_data[days_to_death == "'--", days_to_death := NA]
clinical_data[days_to_birth == "'--", days_to_birth := NA]
clinical_data[year_of_diagnosis == "'--", year_of_diagnosis := NA]
clinical_data[age_at_diagnosis == "'--", age_at_diagnosis := NA]
clinical_data[days_to_last_follow_up == "'--", days_to_last_follow_up := NA]
clinical_data[, days_to_death := as.numeric(days_to_death)]
clinical_data[, days_to_birth := as.numeric(days_to_birth)]
clinical_data[, year_of_diagnosis := as.numeric(year_of_diagnosis)]
clinical_data[, age_at_diagnosis := as.numeric(age_at_diagnosis)]
clinical_data[, days_to_last_follow_up := as.numeric(days_to_last_follow_up)]

# EDA and data clean-up ----
## Age at diagnosis ----
sum(is.na(clinical_data$age_at_index))
range(clinical_data$age_at_index, na.rm = TRUE)

ggplot(clinical_data,
       aes(x = vital_status, group = vital_status, y = age_at_index)) +
  geom_boxplot() +
  labs(y = "Age at Diagnosis", x = "Status") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/boxplot_age_vs_status.png",
       units = "in", height = 5, width = 6)

## Gender ----
table(clinical_data$gender)
ggplot(clinical_data,
       aes(x = gender)) +
  geom_bar(color = "cornflowerblue", fill = "cornflowerblue") +
  labs(y = "Counts", x = "Gender") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/barplot_gender.png",
       units = "in", height = 5, width = 6)
# Issues: unbalanced data, very few male cases for the disease
# I researched it and found that on average 1 in 100 breast cancer case is
# male therefore the 12 to 1085 ratio holds up.
# Decision: remove the male cases
clinical_data <- clinical_data[gender != "male"]

## Ethnicity and race ----
table(clinical_data$ethnicity)
table(clinical_data$race)
# Issues: american indian and alaska native contain 1 case only. Asian has few
# cases with only one event making it an unbalanced category.
# Decision: merge it with not reported and relabel as "other"
clinical_data[, race_orig := race]
clinical_data[race == "american indian or alaska native", race := "other"]
clinical_data[race == "not reported", race := "other"]
clinical_data[race == "asian", race := "other"]

ggplot(clinical_data,
       aes(x = ethnicity)) +
  geom_bar(color = "cornflowerblue", fill = "cornflowerblue") +
  labs(y = "Counts", x = "Ethnicity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/barplot_ethnicity.png",
       units = "in", height = 5, width = 6)

clinical_data[race_orig == "american indian or alaska native",
              race_orig := "american indian \n alaska native"]
clinical_data[race_orig == "black or african american",
              race_orig := "black \n african american"]
ggplot(clinical_data,
       aes(x = race_orig)) +
  geom_bar(color = "cornflowerblue", fill = "cornflowerblue") +
  labs(y = "Counts", x = "Race") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/barplot_race.png",
       units = "in", height = 5, width = 6)

## Year of Diagnosis ----
sum(is.na(clinical_data$year_of_diagnosis))
# Issues: 2 cases with NA year of diagnosis
# Decision: remove them
clinical_data <- clinical_data[!is.na(year_of_diagnosis)]
range(clinical_data$year_of_diagnosis)

ggplot(clinical_data,
       aes(x = vital_status, group = vital_status, y = year_of_diagnosis)) +
  geom_boxplot() +
  labs(y = "Year of Diagnosis", x = "Status") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/boxplot_year_diagnosis_vs_status.png",
       units = "in", height = 5, width = 6)

## Disease Stage ----
table(clinical_data$ajcc_pathologic_stage)
# Issues: missing values and too many stages with few cases.
# Decision: remove the missing values and stage X. Group the data into
# the level before the subdivision into A, B and C.
clinical_data <- clinical_data[clinical_data$ajcc_pathologic_stage != "'--"]
clinical_data <- clinical_data[clinical_data$ajcc_pathologic_stage != "Stage X"]
clinical_data[, stage := "Stage I"]
clinical_data[ajcc_pathologic_stage %in% c("Stage II",
                                           "Stage IIA",
                                           "Stage IIB",
                                           "Stage IIC"),
              stage := "Stage II"]
clinical_data[ajcc_pathologic_stage %in% c("Stage III",
                                           "Stage IIIA",
                                           "Stage IIIB",
                                           "Stage IIIC"),
              stage := "Stage III"]
clinical_data[ajcc_pathologic_stage %in% c("Stage IV"),
              stage := "Stage IV"]

ggplot(clinical_data,
       aes(x = vital_status, group = vital_status, y = stage)) +
  geom_boxplot() +
  labs(y = "Disease Stage", x = "Status") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/boxplot_stage_vs_status.png",
       units = "in", height = 5, width = 6)

# Logistic regression -----
## add a column with 5-year survival logical
clinical_data[, survival_5 := 1]
clinical_data[days_to_death <= 365 * 5, survival_5 := 0]
clinical_data[days_to_last_follow_up <= 365 * 5, survival_5 := 2]
clinical_data_log <- clinical_data[survival_5 %in% c(0, 1)]

## Last check that we have events in all categorical variables
clinical_data_log[, .N, by = .(gender, survival_5)]
clinical_data_log[, .N, keyby = .(race, survival_5)]
clinical_data_log[, .N, keyby = .(stage, survival_5)]

# Standardized continuous variables
clinical_data_log[, age_at_index := scale(age_at_index)]
clinical_data_log[, year_of_diagnosis := scale(year_of_diagnosis)]

## Backward feature selection ----
### Model 1: Full model (age, race, stage and year of diagnosis) ----
model_1 <- glm(survival_5 ~ age_at_index + race + stage + year_of_diagnosis,
               data = clinical_data_log, family = binomial)
summary(model_1)

### Model 2: Remove race since it is not significant ----
model_2 <- glm(survival_5 ~ age_at_index + stage + year_of_diagnosis,
               data = clinical_data_log, family = binomial)
summary(model_2)

### Model 3: Remove year of diagnosis since it is not significant ----
model_3 <- glm(survival_5 ~ age_at_index + stage,
               data = clinical_data_log, family = binomial)
summary(model_3)

## Model Selection ----
# Model 3, the simplest model comes out as the best model
AIC(model_1, model_2, model_3)
BIC(model_1, model_2, model_3)

## Assumptions checks and Results ----
clinical_data_log[, probabilities := predict(model_3, type = "response")]
clinical_data_log[, logit_predict := qlogis(probabilities)]

# results
fp <- clinical_data_log[survival_5 == 0 & probabilities > 0.5, .N] # false positive
tn <- clinical_data_log[survival_5 == 0 & probabilities <= 0.5, .N] # true negative
tp <- clinical_data_log[survival_5 == 1 & probabilities > 0.5, .N] # true positive
fn <- clinical_data_log[survival_5 == 1 & probabilities <= 0.5, .N] # false negative

precision <- tp / (tp + fp)
recall <- tp / (tp + fn)

ggplot(clinical_data_log, aes(x = 1:300,
                              y = probabilities,
                              color = factor(survival_5))) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/prediction_results.png",
       units = "in", height = 5, width = 6)

# linearity of continuous covariates
ggplot(clinical_data_log, aes(x = logit_predict, y = age_at_index)) +
  geom_point(size = 0.5, alpha = 0.5, ) +
  geom_smooth(method = "loess", color = "cornflowerblue") +
  labs(y = "Standardixed Age at Diagnosis") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/linearity_assumption_age.png",
       units = "in", height = 5, width = 6)

# influential points
model_data <- as.data.table(augment(model_3))
model_data[, index := 1:.N]
ggplot(model_data, aes(index, .std.resid)) +
  geom_point(aes(color = factor(survival_5)), alpha = .5) +
  theme_bw() +
  labs(y = "Standardized Residuals") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/influential_points_std_residuals.png",
       units = "in", height = 5, width = 6)

# Accelerated Failure Time Model ----
clinical_data_aft <- copy(clinical_data)
clinical_data_aft[, survival_time := days_to_death]
clinical_data_aft[, status := 1]
clinical_data_aft[is.na(survival_time), status := 0]
clinical_data_aft[is.na(survival_time), survival_time := days_to_last_follow_up]

# clean up?
clinical_data_aft <- clinical_data_aft[survival_time > 0]
clinical_data_aft[, age_at_index_orig := age_at_index]
clinical_data_aft[, age_at_index := scale(age_at_index)]
clinical_data_aft[, year_of_diagnosis := scale(year_of_diagnosis)]

## Backward feature selection ----
### Model 1: Full model (age, race, stage and year of diagnosis) ----
model_aft_1 <- survreg(Surv(time = survival_time, event = status) ~
                       age_at_index + stage + year_of_diagnosis + race,
                     data = clinical_data_aft)

summary(model_aft_1)

### Model 2: Remove race ----
model_aft_2 <- survreg(Surv(time = survival_time, event = status) ~
                         age_at_index + stage + year_of_diagnosis,
                       data = clinical_data_aft)

summary(model_aft_2)

### Model 3: Remove year of diagnosis ----
model_aft_3 <- survreg(Surv(time = survival_time, event = status) ~
                       age_at_index + stage,
                     data = clinical_data_aft)

summary(model_aft_3)

AIC(model_aft_1, model_aft_2, model_aft_3)
BIC(model_aft_1, model_aft_2, model_aft_3)

## Inference ----
probs <- ciTools::add_probs(clinical_data_aft[, .(age_at_index, stage,
                                                  survival_time, status)],
                            model_aft_3, q = 365 * 5,
                            name = c("prob", "lcb", "ucb"),
                            comparison = ">")
probs <- as.data.table(probs)
probs[, age := clinical_data_aft$age_at_index_orig]

ggplot(probs, aes(x = age, y = prob)) +
  geom_point() +
  facet_wrap(~stage, ncol = 4) +
  geom_ribbon(aes(ymin = lcb, ymax = ucb), alpha = 0.5) +
  labs(y = "Probability", x = "Age") +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))

ggsave("figs/aft_survival_by_age_and_stage.png",
       units = "in", height = 3, width = 12)

# Next steps ----
# 1. look deepter into interactions between covariates
# 2. look at the correlation between therapy type and survival time
# 3. look into the Bayesian version of the model to allow for
# uncertainty quantification of the posterior parameters
