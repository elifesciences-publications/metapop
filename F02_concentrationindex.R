#' #############################################################
#' Metapopulation ecology links antibiotic resistance, consumption
#' and patient transfers in a network of hospital wards
#' 
#' Shapiro et al. 2019
#' (c) Jean-Philippe Rasigade, Julie T Shapiro
#' Université Claude Bernard Lyon 1
#' CIRI Inserm U1111
#' 
#' MIT LICENSE

#' #############################################################
#' 
#' Computes the concentration
#' index for infection episodes with each variant;
#' and for the consumption of each class of antibiotics.
#' 
#' Analysis begins using the output dataframe "abundance" from F01_dataprep.R script
#' 

library(iNEXT) # Tested with version 2.0.19  for R version 3.6.0
library(dplyr) # Tested with version 0.8.3
library(data.table)

load("abundance.Rdata") #Output dataframe from F01_dataprep.R script

###############################################################
# Concentration index for infection episodes - TABLE 1.

#Select infection episodes from abundance
conc.bact <- select(abundance, starts_with("N_", ignore.case = FALSE))

#Transpose the data so that rows are variant
t.conc.bact <- transpose(conc.bact)

#Set rownames to variant
t.conc.bact2 <- data.frame(Sites=factor(colnames(conc.bact), levels=colnames(conc.bact)),
                           t.conc.bact)

#Make rownames the variant
rownames(t.conc.bact2) <- t.conc.bact2$Sites

#Eliminate "Site" column
t.conc.bact3 <- t.conc.bact2[,-1]

#Split the matrix into a series of lists (vector)
#Each bacterial taxa becomes its own list containing a vector of abundances in each ward
t.conc.bact4 <- as.list(setNames(split(t.conc.bact3, seq(nrow(t.conc.bact3))), rownames(t.conc.bact3)))

#Use the ChaoSimpson function in the iNEXT package to calculate the 
#asymptotic dispersion index and 95% confidence intervals for each bacterial taxon + resistance profile
# Parameter B controls no. of bootstrap replicates. Set B = 200 for fast estimation but
# use B = 10,000 to obtain stable estimates for publication.
bact.disp.index <- ChaoSimpson(t.conc.bact4, datatype="abundance", B = 200)

# Note that this is a bootstrap procedure, confidence intervals can vary slightly
# from one run to another.

# Take complement of dispersion index to obtain concentration index in percents
bact.conc.index <- data.table(bact.disp.index)[
  , lapply(.SD, function(x) (1 - x)*100), .SDcols = names(bact.disp.index)[-3]
  ][ , variant := row.names(bact.disp.index)]
print(bact.conc.index[])

###################################################################################################
# Concentration index for the consumption of antibiotic classes - TABLE 2

#Select consumption of each antibiotic from "abundance" dataframe
conc.atb <- abundance %>%
  select(starts_with("ddd_", ignore.case = FALSE)) %>%
  select(-ddd_total)

#Transpose the data so that rows are antibiotic classes
t.conc.atb <- transpose(conc.atb)

#Set rownames to the antibiotic classes
t.conc.atb2 <- data.frame(Sites=factor(colnames(conc.atb), levels=colnames(conc.atb)),
                          t.conc.atb)

#Make rownames the antibiotic classes
rownames(t.conc.atb2) <- t.conc.atb2$Sites

#Eliminate "Site" column
t.conc.atb3 <- t.conc.atb2[,-1]

#Split the matrix into a series of lists (vector)
#Each antibiotic becomes its own list containing a vector of use in each ward
t.conc.atb4 <- as.list(setNames(split(t.conc.atb3, seq(nrow(t.conc.atb3))), rownames(t.conc.atb3)))

#Use the ChaoSimpson function in the iNext package to calculate the 
#asymptotic dispersion index and 95% confidence intervals for each class of antibiotics
# Parameter B controls no. of bootstrap replicates. Set B = 200 for fast estimation but
# use B = 10,000 to obtain stable estimates for publication.
atb.disp.index <- ChaoSimpson(t.conc.atb4, datatype="abundance", B = 200)

# Take complement of dispersion index to obtain concentration index in percents
atb.conc.index <- data.table(atb.disp.index)[
  , lapply(.SD, function(x) (1 - x)*100), .SDcols = names(atb.disp.index)[-3]
  ][ , group := row.names(atb.disp.index)]
print(atb.conc.index[])
