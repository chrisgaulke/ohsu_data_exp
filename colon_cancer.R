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
library(randomForest)

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

#now lets make a data frame without those pesky extra columns

rf.cancer <- colon_cancer_mb_cpm.df
#make sensible row names
rownames(rf.cancer) <- rf.cancer$ENTITY_STABLE_ID

#remove non numeric tax column
rf.cancer <- rf.cancer[,5:587]

#transpose table so columns are now rows to facilitate analysis downstream
rf.cancer <- t(rf.cancer)

#always double check!
all(rownames(colon_cancer_patient.metadata) == rownames(rf.cancer))


# ANALYSIS: AGE -----------------------------------------------------------

#Lets looks to see how age of diagnosis is distributed in these data
#recall from the full metadata file that "AGE" is actually age of diagnosis
#For now we will assume all of these samples have cancer

hist(colon_cancer_patient.metadata$AGE)

#lets also sweep out some summary stats here

mean(colon_cancer_patient.metadata$AGE)

#Get tukeys five number summary of data. This will include (in this order)
#the minimum, lower-hinge, median, upper-hinge, and the maximum
fivenum(colon_cancer_patient.metadata$AGE)

#Now we can set up tentative age bins
#early < 60 ; normal 60 - 70; late > 70

colon_cancer_patient.metadata$GROUP <-

  factor(

    cut(
      colon_cancer_patient.metadata$AGE, c(0,60,70,100)
      )

    )

#now make these numbers more readable
levels(colon_cancer_patient.metadata$GROUP) <- c("0-60", "61-70", ">70")

#cut can be tricky so always check to make sure things are as you want them

View(colon_cancer_patient.metadata[,c("AGE", "GROUP")])

#now lets take a look at how many samples are in each group
table(colon_cancer_patient.metadata$GROUP)

#finally lets looks to see how well a classifier does at binning colon cancer
#patients

randomForest(y = colon_cancer_patient.metadata$GROUP, x = rf.cancer)
