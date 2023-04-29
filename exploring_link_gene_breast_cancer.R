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
for(doc in clinical_data$case_id) {

  data_dir <- paste0("data/genes/",
                     gene_file_dictionary[case_id == doc, file_id],
                     "/")[1]
  data_file <- gene_file_dictionary[case_id == doc, file_name][1]

  # check that the case-ID has gene expression folder
  if(is.na(data_file)){
    next
  }

  data_case <- as.data.table(read.table(paste0(data_dir, data_file),
                                        sep = "\t",
                                        header = TRUE))
  # add column with case_id
  data_case[, case_id := doc]
  data_case <- data_case[, .(case_id, gene_name, gene_type,
                             tpm_unstranded, fpkm_unstranded)]

  # filter our first 4 rows and keep protein-coding genes only
  data_case <- data_case[!is.na(tpm_unstranded) & gene_type == "protein_coding"]

  # remove duplicates and keep first row (based on gene_name)
  data_case <- data_case[!duplicated(data_case, by = "gene_name")]

  if(doc == clinical_data$case_id[1]){

    gene_data <- data_case

  } else {

    gene_data <- rbind(gene_data, data_case)
  }

}

## attach/merge vital status to gene data table
gene_data <- merge(gene_data,
                   clinical_data[, .(case_id, vital_status)],
                   by = "case_id")

# Analyse the data ----

