# Title and author information --------------------------------------------
#!/usr/bin/R

###################
#                 #
# colon_cancer.R #
#                 #
###################

#Title: Colon cancer public data analysis
#
#Copyright (C) 2021-2022  Christopher A. Gaulke
#author contact: chris.gaulke@gmail.com
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#For a copy of the GNU General Public License see
#<http://www.gnu.org/licenses/>.

# Purpose: Data wrangling and viz


# SET ENVIRONMENT ---------------------------------------------------------

options(stringsAsFactors = F)
library(ggplot2)
library(vegan)
#library(randomForest)

# DATA: IMPORT DATA -------------------------------------------------------

colon_cancer_mb_cpm.df <- read.table(file = "../coadread_tcga_pan_can_atlas_2018/data_microbiome.txt",
                                     row.names = NULL,
                                     sep = "\t",
                                     header = T,
                                     quote ="",
                                     comment.char = "#",
                                     as.is = T)

colon_cancer_patient.metadata <- read.table("../coadread_tcga_pan_can_atlas_2018/data_clinical_patient.txt",
                                       row.names = 1,
                                       skip = 4,
                                       sep = "\t",
                                       header = T,
                                       quote ="",
                                       comment.char = "#"
                                       )

colon_cancer_sample.metadata <- read.table("../coadread_tcga_pan_can_atlas_2018/data_clinical_sample.txt",
                                           row.names = 1,
                                           skip = 4,
                                           sep = "\t",
                                           header = T,
                                           quote ="",
                                           comment.char = "#")


# ANALYSIS: DATA WRANGLE ------------------------------------------------------

# Are all samples in both files?
all(rownames(colon_cancer_patient.metadata) %in%
      rownames(colon_cancer_sample.metadata))


# Are all samples in the same order in both files?
all(rownames(colon_cancer_patient.metadata) ==
      rownames(colon_cancer_sample.metadata))

#df[rows,cols]
#must include both rows and columns to reorder

colon_cancer_patient.metadata <-
  colon_cancer_patient.metadata[rownames(colon_cancer_sample.metadata),]

# Are all samples in the same order in both files?
all(rownames(colon_cancer_patient.metadata) ==
      rownames(colon_cancer_sample.metadata))


#fix colon_cancer_mb_cpm.df colnames to match metadata rownames
# don't forget to escape wildcard char "." !

colnames(colon_cancer_mb_cpm.df) <-
  gsub(x = colnames(colon_cancer_mb_cpm.df), pattern = "\\.", "-")

#see if any sample names match
colnames(colon_cancer_mb_cpm.df)[5:587] %in% rownames(colon_cancer_patient.metadata)

#replace trailing -01 and see if any sample names match

all(gsub(x = colnames(colon_cancer_mb_cpm.df), pattern = "-01", "")[5:587] %in%  rownames(colon_cancer_patient.metadata))

#fix colnames again
colnames(colon_cancer_mb_cpm.df) <-
  gsub(x = colnames(colon_cancer_mb_cpm.df), pattern = "-01", "")

#filter metadata

colon_cancer_patient.metadata <-
  colon_cancer_patient.metadata[which(rownames(colon_cancer_patient.metadata) %in%
                                      colnames(colon_cancer_mb_cpm.df)),]


colon_cancer_sample.metadata <-
  colon_cancer_sample.metadata[which(rownames(colon_cancer_sample.metadata) %in%
                                        colnames(colon_cancer_mb_cpm.df)),]

#order metadata
colon_cancer_patient.metadata <-
  colon_cancer_patient.metadata[colnames(colon_cancer_mb_cpm.df)[5:587],]


colon_cancer_sample.metadata <-
  colon_cancer_sample.metadata[colnames(colon_cancer_mb_cpm.df)[5:587],]

all(rownames(colon_cancer_sample.metadata) == colnames(colon_cancer_mb_cpm.df)[5:587])
all(rownames(colon_cancer_patient.metadata) == colnames(colon_cancer_mb_cpm.df)[5:587])


# ANALYSIS: Viz -----------------------------------------------------------

hist(colon_cancer_patient.metadata$AGE)
table(colon_cancer_patient.metadata$RACE)

