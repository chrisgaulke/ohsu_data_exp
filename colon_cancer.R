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

#Questions


#1. What are the differences between early onset (<50 years) and late onset (>65)
#   microbiomes
#2. Are there microbiome differences in patients who have high CFD,
#   C7 complement protein expression- this can be divided into two groups
#   by median expression or quartiles and compare top 25% vs bottom 25%
#   (regardless of age, other variables)


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

#gene expression data for C7 and CFD
cancer_gene.df <- read.table(file = "../coadread_tcga_pan_can_atlas_2018/modified_gene_data/gene_quantification_C7_CFD.txt",
                             row.names = 1,
                             sep = "\t",
                             header = T,
                             quote ="",
                             comment.char = "#",
                             as.is = T)

cancer_gene.df$Entrez_Gene_Id <- NULL #do not want


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


# ANALYSIS: DATA WRANGLE GENES --------------------------------------------

#start with genes and fix names
colnames(cancer_gene.df) <-
  gsub(x = colnames(cancer_gene.df), pattern = "\\.", "-")

#fix colnames again
colnames(cancer_gene.df) <-
  gsub(x = colnames(cancer_gene.df), pattern = "-01", "")

cancer_gene.df <- cancer_gene.df[,
                    which(
                      colnames(cancer_gene.df) %in%
                        rownames(colon_cancer_patient.metadata)
                      )
                    ]

cancer_gene.df <- cancer_gene.df[,rownames(colon_cancer_patient.metadata)]

cancer_gene.df <- t(cancer_gene.df)
cancer_gene.df <- as.data.frame(cancer_gene.df)

#add to metadata

colon_cancer_patient.metadata$CFD <- cancer_gene.df$CFD
colon_cancer_patient.metadata$C7 <- cancer_gene.df$C7

# ANALYSIS: BETA DIVERSITY  ------------------------------------------------------

#Lets looks to see how age of diagnosis is distributed in these data
#recall from the full metadata file that "AGE" is actually age of diagnosis
#For now we will assume all of these samples have cancer

hist(colon_cancer_patient.metadata$AGE)

#note that this is equiv to:

hist(colon_cancer_patient.metadata[,"AGE"])

#lets also sweep out some summary stats here

mean(colon_cancer_patient.metadata$AGE)

#Get tukeys five number summary of data. This will include (in this order)
#the minimum, lower-hinge, median, upper-hinge, and the maximum
fivenum(colon_cancer_patient.metadata$AGE)

#Now we can set up tentative age bins
#early < 50 ; normal 51 - 79; late > 79

colon_cancer_patient.metadata$GROUP <-
    cut(
      colon_cancer_patient.metadata$AGE, c(0,50,65,100)
      )


#now make these numbers more readable
levels(colon_cancer_patient.metadata$GROUP) <- c("0-50", "51-65", ">65")

#cut can be tricky so always check to make sure things are as you want them

#View(colon_cancer_patient.metadata[,c("AGE", "GROUP")])

#now lets take a look at how many samples are in each group
table(colon_cancer_patient.metadata$GROUP)

#lets have a look at how

cancer.prcomp <- prcomp(rf.cancer,scale =F, center = T )

#make a skree plot to identify important dimensions
plot(cancer.prcomp)

cancer_prcomp.df <- cancer.prcomp$x
cancer_prcomp.df <- as.data.frame(cancer_prcomp.df)
cancer_prcomp.df$GROUP <- colon_cancer_patient.metadata$GROUP
cancer_prcomp.df$C7 <- colon_cancer_patient.metadata$C7
cancer_prcomp.df$CFD <- colon_cancer_patient.metadata$CFD
cancer_prcomp.df$CENTER <- colon_cancer_patient.metadata$CENTER

cancer_prcomp_age.plot <- ggplot(cancer_prcomp.df, aes(x = PC1,
                                                       y = PC2,
                                                       color = GROUP))

cancer_prcomp_age.plot +
  geom_point()

# we can color the points by gene expression to get a high level view of
cancer_prcomp_age.plot <- ggplot(cancer_prcomp.df, aes(x = PC1,
                                                       y = PC2,
                                                       color = C7))

cancer_prcomp_age.plot +
  geom_point()

cancer_prcomp_age.plot <- ggplot(cancer_prcomp.df, aes(x = PC1,
                                                       y = PC2,
                                                       color = CFD))

cancer_prcomp_age.plot +
  geom_point()

#Center
cancer_prcomp_age.plot <- ggplot(cancer_prcomp.df, aes(x = PC1,
                                                       y = PC2,
                                                       color = CENTER))

cancer_prcomp_age.plot +
  geom_point()

#adonis note because these data are prenormalized we will get a warning here
#we should approach this result with skepticism.
cancer_age.adonis <- adonis(rf.cancer ~ colon_cancer_patient.metadata$GROUP, permutations = 1000)
cancer_age_C7.adonis <- adonis(rf.cancer ~ colon_cancer_patient.metadata$C7, permutations = 1000)
cancer_age_CFD.adonis <- adonis(rf.cancer ~ colon_cancer_patient.metadata$CFD, permutations = 1000)

cancer_age.adonis #not sig
cancer_age_C7.adonis #sig
cancer_age_CFD.adonis #sig

#Since we can't trust adonis but together with the PCA it seems like there may
# be an association between C7 or CFD and the microbiome we can try to determine
# if the vector of gene expression values correlate with the dimensions of the
# ordination. We can do this with a simple cor.test or with envfit. We will try
# both here

cancer_metadata.envfit <- envfit(cancer.prcomp,
                        colon_cancer_patient.metadata[,c("GROUP",
                                                         "C7",
                                                         "CFD",
                                                         "CENTER"
                                                         )],
                        na.rm = T,
                        permutations = 5000,choices = c(1,2))

cancer_metadata.envfit

#we can also just use independent correlations with the different PCs we will
# limit this to the first three dims

cor.test(cancer_prcomp.df$PC1,cancer_prcomp.df$C7, method = "spearman")# n.s
cor.test(cancer_prcomp.df$PC2,cancer_prcomp.df$C7, method = "spearman")# **
cor.test(cancer_prcomp.df$PC3,cancer_prcomp.df$C7, method = "spearman")# n.s

cor.test(cancer_prcomp.df$PC1,cancer_prcomp.df$CFD, method = "spearman")# n.s
cor.test(cancer_prcomp.df$PC2,cancer_prcomp.df$CFD, method = "spearman")# **
cor.test(cancer_prcomp.df$PC3,cancer_prcomp.df$CFD, method = "spearman")# n.s


# ANALYSIS: RANDOM FORESTS -----------------------------------------------


#finally lets looks to see how well a classifier does at binning colon cancer
#patients

age.rf <- randomForest(y = colon_cancer_patient.metadata$GROUP, x = rf.cancer)

#take a look at performance
age.rf

#classifier does pretty awful. The OOB error is high and the class error is
#~1 for two groups. This would suggest that there aren't strong differences
# between age groups with regard to microbiome composition.
