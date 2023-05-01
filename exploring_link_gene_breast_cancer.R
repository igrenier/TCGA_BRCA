# Title: Understanding which genes are most associated with patient survival
# Author: @igrenier
# Date: 04/28/2023

# Libraries ----
library(data.table)
library(ggplot2)

# Helper functions ----
# Build linear formulas given a dependent variable and a vector of
# independent variables
build_cor_formula <- function(dep_var, ind_var_vec) {

  N <- length(ind_var_vec)
  formula <- paste0(dep_var, " ~ ")
  for(k in 1:(N - 1)){

    formula <- paste0(formula, ind_var_vec[k], " + ")

  }
  formula <- paste0(formula, ind_var_vec[N])
}

# Read the data ----
## Clinical Data ----
## read file
clinical_data <- fread("data/clinical.tsv", sep ="\t", header = TRUE)
clinical_data <- as.data.table(clinical_data)

## reduce to needed variables for project and remove duplicates
clinical_data <- clinical_data[, .(case_id, vital_status,
                                   days_to_last_follow_up,
                                   days_to_death)]
clinical_data <- unique(clinical_data)

## clean data format
clinical_data[days_to_last_follow_up == "'--", days_to_last_follow_up := NA]
clinical_data[, days_to_last_follow_up := as.numeric(days_to_last_follow_up)]
clinical_data[days_to_death == "'--", days_to_death := NA]
clinical_data[, days_to_death := as.numeric(days_to_death)]

## remove patients with less than 3 month history
clinical_data <- clinical_data[days_to_last_follow_up > 90 |
                                 days_to_death > 90]
clinical_data[, survival_time := days_to_last_follow_up]
clinical_data[vital_status == "Dead", survival_time := days_to_death]

## create subset to look at 5-year survival rate only
clinical_data_5_year <- clinical_data[days_to_last_follow_up > 5 * 365 |
                                        is.na(days_to_last_follow_up)]
clinical_data_5_year[, status := 1]
clinical_data_5_year[days_to_death < 5 * 365, status := 0]

### how many cases for each status?
table(clinical_data_5_year$vital_status) # 249/68

## Gene Data ----
## read and extract information from file dictionary
gene_file_dictionary <- as.data.table(read.csv("data/gene_file_dictionary.csv",
                                               header = TRUE))

## read each folder and append to common data table with added case id
## Note 1: we only extract cases with 3 months of history
## Note 2: we run into the issue of cases IDs having 2 or more gene folders,
##         we will ignore and keep the first instance only. To be discussed.
## Note 3: some case ids do not have the gene expression folder, we skip them.
n_cases <- clinical_data[, .N]
for(cs in 1:n_cases) {

  data_case_id <- clinical_data$case_id[cs]
  data_dir <- paste0("data/genes/",
                     gene_file_dictionary[case_id == data_case_id, file_id],
                     "/")[1]
  data_file <- gene_file_dictionary[case_id == data_case_id, file_name][1]

  # check that the case-ID has gene expression folder
  if(is.na(data_file)){
    next
  }

  data_case <- as.data.table(read.table(paste0(data_dir, data_file),
                                        sep = "\t",
                                        header = TRUE))
  # add column with case_id
  data_case <- data_case[, .(gene_name, gene_type,
                             tpm_unstranded)]

  # filter our first 4 rows and keep protein-coding genes only
  data_case <- data_case[!is.na(tpm_unstranded) & gene_type == "protein_coding"]

  # remove duplicates and keep first row (based on gene_name)
  data_case <- data_case[!duplicated(data_case, by = "gene_name")]

  if(cs == 1){

    gene_data_row <- data_case$tpm_unstranded
    names(gene_data_row) <- data_case$gene_name

    gene_data <- data.table("case_id" = cs)
    gene_data <- cbind(gene_data, t(gene_data_row))

  } else {

    gene_data_row <- data_case$tpm_unstranded
    names(gene_data_row) <- data_case$gene_name
    gene_data <- rbind(gene_data, t(c("case_id" = cs, gene_data_row)))

  }

}

## attach/merge survival time to gene data table
gene_data[, case_id := clinical_data$case_id[case_id]]
gene_data <- merge(gene_data,
                   clinical_data[, .(case_id, survival_time)],
                   by = "case_id")

gene_data_5_year <- merge(gene_data,
                          clinical_data_5_year[, .(case_id, status)],
                          by = "case_id")

# Analyse the 5 year data ----
## remove the case id
gene_data_5_year[, case_id := NULL]

## make a copy of the gene without status
gene_data_5_year_0 <- copy(gene_data_5_year)
gene_data_5_year_0[, status := NULL]
gene_data_5_year_0[, survival_time := NULL]

## remove the columns where the marker is 0 for all cases
gene_data_all_0 <- colnames(gene_data_5_year_0)[colSums(gene_data_5_year_0) == 0]
gene_data_5_year_0 <- gene_data_5_year_0[, (gene_data_all_0) := NULL]
gene_data_5_year <- gene_data_5_year[, (gene_data_all_0) := NULL]

## Correlation ----
corr_test_5 <- apply(gene_data_5_year_0, 2,
                   function(x) cor(x, y = gene_data_5_year$status))
corr_test_5 <- sort(corr_test_5)
corr_5_dt <- data.table("correlation" = corr_test_5)

## largest negative correlation
neg_cor_names <- names(corr_test_5[1:5])
## largest positive correlation
pos_cor_names <- names(corr_test_5[length(corr_test_5) - 4:0])

ggplot(corr_5_dt, aes(x = correlation)) +
  geom_histogram(fill = "grey", color = "black") +
  theme_bw() +
  labs(x = "Correlation") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/correlation.png",
       units = "in", height = 5, width = 6)

ggplot(gene_data_5_year[, .(status, get(pos_cor_names[5]))],
       aes(x = factor(status), group = factor(status), y = V2)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(labels = c('Dead','Alive')) +
  labs(y = pos_cor_names[5]) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))

ggsave("figs/pos_cor_example.png",
       units = "in", height = 5, width = 6)

ggplot(gene_data_5_year[, .(status, get(neg_cor_names[1]))],
       aes(x = factor(status), group = factor(status), y = V2)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(labels = c('Dead','Alive')) +
  labs(y = neg_cor_names[5]) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))

ggsave("figs/neg_cor_example.png",
       units = "in", height = 5, width = 6)

### Logistic regression

# Full Model
formula_1 <- build_cor_formula("status", c(neg_cor_names, pos_cor_names))
mymodel_1 <- glm(formula_1, data = gene_data_5_year)
summary(mymodel_1)

# Model 2: remove the 2nd and 3rd positive correlation
formula_2 <- build_cor_formula("status",
                               c(neg_cor_names,
                                 pos_cor_names[c(1, 4:5)]))
mymodel_2 <- glm(formula_2, data = gene_data_5_year)
summary(mymodel_2)

# Analyse the continuous survival time ----
## remove the case id
gene_data[, case_id := NULL]

## make a copy of the gene without status
gene_data_0 <- copy(gene_data)
gene_data_0[, survival_time := NULL]

## remove the columns where the marker is 0 for all cases
gene_data_all_0 <- colnames(gene_data_0)[colSums(gene_data_0) == 0]
gene_data_0 <- gene_data_0[, (gene_data_all_0) := NULL]
gene_data <- gene_data[, (gene_data_all_0) := NULL]

## Correlation ----
corr_test <- apply(gene_data_0, 2,
                   function(x) cor(x, y = log(gene_data$survival_time)))
corr_dt <- data.table("gene"= names(corr_test),
                      "correlation" = corr_test)
corr_test <- sort(corr_test)

## largest negative correlation
neg_cor_names <- names(corr_test[1:5])
## largest positive correlation
pos_cor_names <- names(corr_test[length(corr_test) - 4:0])

ggplot(corr_dt, aes(x = correlation)) +
  geom_histogram(fill = "grey", color = "black") +
  theme_bw() +
  labs(x = "Correlation") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/correlation_continuous.png",
       units = "in", height = 5, width = 6)

ggplot(gene_data[, .(survival_time, get(pos_cor_names[5]))],
       aes(x = log(survival_time), y = V2)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method='lm', formula= y~x, color = "cornflowerblue") +
  labs(y = pos_cor_names[5], x = "log(Survival Time)") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))

ggsave("figs/pos_cor_example_continuous.png",
       units = "in", height = 5, width = 6)

ggplot(gene_data[, .(survival_time, get(neg_cor_names[1]))],
       aes(x = log(survival_time), y = V2)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method='lm', formula= y~x, color = "cornflowerblue") +
  labs(y = neg_cor_names[1], x = "log(Survival Time)") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))

ggsave("figs/neg_cor_example_continuous.png",
       units = "in", height = 5, width = 6)
## Rank Correlation ----
corr_test <- apply(gene_data_0, 2,
                   function(x) cor.test(x, y = log(gene_data$survival_time),
                                        method = "spearman")$estimate)
corr_dt[, "rank_rho" := corr_test]
corr_test <- sort(corr_test)

ggplot(corr_dt, aes(x = rank_rho)) +
  geom_histogram(fill = "grey", color = "black") +
  theme_bw() +
  labs(x = "Spearman's Rho") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))
ggsave("figs/rank_rho_correlation_continuous.png",
       units = "in", height = 5, width = 6)

## largest 10 comparison
setorder(corr_dt, cols = "correlation")
corr_dt[1:10]
setorder(corr_dt, cols = "rank_rho")
corr_dt[1:10]
setorder(corr_dt, cols = -"correlation")
corr_dt[1:10]
setorder(corr_dt, cols = -"rank_rho")
corr_dt[1:10]

## Principal Component Analysis ----
pc <- prcomp(gene_data_0,
             center = TRUE,
             scale. = TRUE)

### Extract the variance information
pca_dt <- data.table("sdev" = pc$sdev)
pca_dt[, component := 1:.N]
pca_dt[, proportion_var := sdev^2 / sum((sdev^2))]
pca_dt[, cum_proportion := cumsum(proportion_var)]

### Plot the diminishing return of components
ggplot(pca_dt[component < 50], aes(x = component, y = cum_proportion)) +
  geom_point() +
  theme_bw() +
  labs(y = "Cumulative Proportion of Variance",
       x = "Principal Component Index") +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "cornflowerblue",
                                    fill = NA, size  = 1))

ggsave("figs/pca.png",
       units = "in", height = 5, width = 6)

# Could be done in the future
# 1. Archetypal analysis: instead of grouping by looking for the mean of a
#   cluster, we look at finding the boundaries/extreme cases
