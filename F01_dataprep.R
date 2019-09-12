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
#' 
#' #############################################################
#' 
#' Performs computations for connectivity and incidence control.
#' 
#' The output dataframe "abundance" is used in downstream analyses.
#' 

library(data.table) # Tested with version 1.12.2 for R version 3.6.0

#############################################################################################
#' Load data including:
#' 
#' ward_incidence -- the number of patients tested positive
#'   for each variant during the study period. Rows are wards (n = 574), columns are variants.
#'   
#' ward_prevalence -- the proportion of patients tested positive 
#'   at least once for each variant during the study period. Rows are wards (n = 574), columns are variants.
#'   
#' ward_metadata -- for each ward, the ward size in column "n_beds" and ward type in column "PatStat" 
#'   as a proxy to patient fragility, where 2 denotes ICUs and blood cancer units,
#'   1 denotes progressive care units and 0 denotes other wards. (n = 448 wards)
#'   
#' ward_atb_use -- for each ward, the consumption of each antibiotic group in defined daily doses. Column "ddd_total"
#'   is the consumption of all antibiotics (ATC class J01).
#'   (n = 342 wards. HCL wards not in this dataframe did not prescribe antibiotics)
#'   
#' sample_type -- the probability that each sample type (location) is positive for a variant,
#'   aggregated over all wards. Rows are sample types, columns are variants.
#'   
#' sample_locations -- distribution of sample types/locations in wards. For each ward and sample location,
#'   column "n_patients" is the number of patients in which "n_samples" samples were taken.
#'   
#' transfers -- the no. of patients transferred from the upstream ward in column "ward_from"
#'   to the downstream ward in column "ward_to".
#'
#'###########################################################
#' Ward identifiers are specific to the study.

print(load("rawdata.Rdata"))



#' #########################################
# Incidence control

# Keep track of all variant names for convenience
groupNames <- names(ward_incidence)[-1]

# Append per-variant probability to each sample location, written P(variant|location) 
sample_locations <- merge(sample_locations, sample_type, by = "sample_type")

# For each group of patients, compute the probability that all samples remained negative
sample_locations <- sample_locations[ , lapply(.SD, function(x) prod((1 - x)^n_samples) ), by = .(ward, n_patients), .SDcols = groupNames]

# Take the complement of this probability and aggregate over patients in each ward to obtain the expected
# number of patients tested positive at least once.
incidence_control <- sample_locations[ , lapply(.SD, function(x) sum( n_patients * (1 - x)) ), by = .(ward), .SDcols = groupNames ]


#' ########################################
#' Connectivity
#' 

# Append per-variant proportion of positive patients in each upstream ward (i.e., "ward_from").
transfers <- merge(transfers, ward_prevalence, by.x = "ward_from", by.y = "ward")

# Multiply this probability with the number of transferred patients to obtain the expected
# number of positive transferred patients. For each downstream ward (i.e., "ward_to"), add the
# expected number of positive transferred patients from all upstream wards to obtain the total
# expected number of positive admitted patients, i.e. the connectivity.
connectivity <- transfers[ , lapply(.SD, function(x) sum(x * n_transferred)), by = .(ward_to), .SDcols = groupNames]
setnames(connectivity, "ward_to", "ward")

#' ########################################
# Combined dataset

# Prefix each variant name with "N_" for incidence, "C_" for incidence control, "S_" for connectivity
N_ <- data.table(ward_incidence)
C_ <- data.table(incidence_control)
S_ <- data.table(connectivity)

setnames(N_, groupNames, paste("N_", groupNames, sep = ""))
setnames(C_, groupNames, paste("C_", groupNames, sep = ""))
setnames(S_, groupNames, paste("S_", groupNames, sep = ""))

# Additional data related to ward characteristics and antibiotic consumption
W_ <- data.table(ward_metadata)
D_ <- data.table(ward_atb_use)

abundance <- N_[]
abundance <- merge(abundance, C_, all.x = T)
abundance <- merge(abundance, S_, all.x = T)
abundance <- merge(abundance, D_, all.x = T)
abundance <- merge(abundance, W_, all.x = F) # Restrict inclusion to wards within HCL with known activity (n = 357)

# Wards with zero connectivity and/or consumption
# are excluded from the S_ and D_ dataframes by design. Their respective
# values should be zero but are replaced with NA during left full joint.

# Replace missing values with zeroes for connectivity and consumption
zeroCols  <- c(names(S_)[-1], names(D_)[-1])
abundance <- abundance[ , (zeroCols) := lapply(.SD, function(x) {x[is.na(x)] <- 0.0; return(x)}), .SDcols = zeroCols]

# Check that no NA remains
tmp <- lapply(abundance, function(x) stopifnot(all(!is.na(x))))

save(abundance, file = "abundance.Rdata")



