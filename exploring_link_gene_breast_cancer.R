# Title: Understanding which genes are most associated with patient survival
# Author: @igrenier
# Date: 04/28/2023

# Libraries ----
library(data.table)
library(ggplot2)

# Read the data ----
## Clinical Data ----
## read file
clinical_data <- fread("data/clinical.tsv", sep ="\t", header = TRUE)
clinical_data <- as.data.table(clinical_data)

## reduce to needed variables for project and remove duplicates
clinical_data <- clinical_data[, .(case_id, vital_status,
                                   days_to_last_follow_up)]
clinical_data <- unique(clinical_data)

## clean data format
clinical_data[days_to_last_follow_up == "'--", days_to_last_follow_up := NA]
clinical_data[, days_to_last_follow_up := as.numeric(days_to_last_follow_up)]

## clean data to remove patients with less than 5 year survival?
clinical_data <- clinical_data[days_to_last_follow_up > 5 * 365 |
                                 is.na(days_to_last_follow_up)]

### how many cases for each status?
table(clinical_data$vital_status) # 201/116 (good balance for analysis)

## Gene Data ----
## read and extract information from file dictionary
gene_file_dictionary <- as.data.table(read.csv("data/gene_file_dictionary.csv",
                                               header = TRUE))

## read each folder and append to common data table with added case id
## Note 1: we only extract cases with available vital status
## Note 2: we run into the issue of 124 cases IDs having 2 or more gene folders,
##         we will ignore and keep the first instance only. To be discussed.
## Note 3: one case id does not have the gene expression folder, we skip it.
n_cases <- clinical_data[, .N]
for(cs in 1:n_cases) {

  data_case <- clinical_data$case_id[cs]
  data_dir <- paste0("data/genes/",
                     gene_file_dictionary[case_id == data_case, file_id],
                     "/")[1]
  data_file <- gene_file_dictionary[case_id == data_case, file_name][1]

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

## attach/merge vital status to gene data table
gene_data[, case_id := clinical_data$case_id[case_id]]
gene_data <- merge(gene_data,
                   clinical_data[, .(case_id, vital_status)],
                   by = "case_id")

# Analyse the data ----
## remove the case id
gene_data[, case_id := NULL]
gene_data[vital_status == "Dead", vital_status := 0]
gene_data[vital_status == "Alive", vital_status := 1]
gene_data[, vital_status := as.numeric(vital_status)]
gene_data_0 <- copy(gene_data)
gene_data_0[, vital_status := NULL]
gene_data_all_0 <- colnames(gene_data_0)[colSums(gene_data_0) == 0]
gene_data_0 <- gene_data_0[, (gene_data_all_0) := NULL]
gene_data <- gene_data[, (gene_data_all_0) := NULL]

# correlation
corr_test <- apply(gene_data_0, 2, function(x) cor(x, y = gene_data$vital_status))
corr_test <- sort(corr_test)

# five smallest
neg_cor_names <- corr_test[1:5]
# five largest
pos_cor_names <- corr_test[length(corr_test) - 4:0]

# logistic regression
names(neg_cor_names)
names(pos_cor_names)
mymodel_1 <- glm(vital_status ~ EDA2R  + DEFB132 + ZNF705G + PTCHD4 + TYRO3 +
                   IL12B + TGFB1 + NFKBIA + CD74 + RGS1,
               data = trg)
summary(mymodel_1)

# pca
pc <- prcomp(gene_data_0,
             center = TRUE,
             scale. = TRUE)
summary(pc)
trg <- predict(pc, gene_data)
trg <- data.frame(trg, gene_data)

mymodel <- glm(vital_status ~ , data = gene_data)
summary(mymodel)

# Could be done in the future
# 1. Archetypal analysis: instead of grouping by looking for the mean of a
#   cluster, we look at finding the boundaries/extreme cases
