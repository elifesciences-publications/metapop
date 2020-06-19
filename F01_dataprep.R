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
#' sample_sites -- distribution of sample sites taken in wards. Columns SST, GUT, DEEP, RT, DEV and OTHERS
#'  indicate the number of samples taken at each site in a patient. Column N indicates the no. of patients
#'  with the same distribution of repeated samples.#'   
#'   
#' transfers -- the no. of patients transferred from the upstream ward in column "ward_from"
#'   to the downstream ward in column "ward_to".
#'
#'###########################################################
#' Ward identifiers are specific to the study.

print(load("rawdata.Rdata"))



#' #########################################
# Incidence control

# Melt sample probabilities
sample_type_melt <- melt(sample_type[!is.na(sample_type)],
          id.vars = "sample_type",
          measure.vars = setdiff(names(sample_type), "sample_type"),
          variable.name = "variant",
          value.name = "prob")

# Melt sample site distribution and add a dummy index 'row' per row
sample_sites_melt <- melt(sample_sites[, row := 1:nrow(sample_sites)],
                          id.vars = c("ward", "N", "row"),
                          measure.vars = setdiff(names(sample_sites), c("ward", "N", "row")),
                          variable.name = "sample_type",
                          value.name = "count"
                          )[count > 0][order(row)]


# Merge with per-sample probability of being positive for each variant
sample_probs <- merge(sample_sites_melt, sample_type_melt,
                      by = "sample_type", allow.cartesian = TRUE, sort = FALSE)

# For each variant, compute the probability that each series of samples remained negative
sample_probs[ , neg_prob := (1 - prob)^count]

# For each variant and row, compute the probability that all samples remained negative
sample_probs <- sample_probs[ ,  .(neg_prob = prod(neg_prob)), by = .(ward, N, row, variant)]

# For each variant and ward, compute the expected incidence as the weighted sum of probability
# that >1 sample was positive (complement of the probability of negativity)
incidence_control <- sample_probs[, .(incidence = sum(N*(1 - neg_prob))), by = .(ward, variant)]

# CHECKPOINT Distribution per variant
incidence_control[ , sum(incidence), by = .(variant)]

# dcast per ward and variant
incidence_control <- dcast(incidence_control, ward ~ variant, value.var = "incidence")

#' ########################################
#' Connectivity
#' 

# Keep track of all variant names for convenience
groupNames <- names(ward_incidence)[-1]

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
abundance <- merge(abundance, C_, all.x = T, by = "ward")
abundance <- merge(abundance, S_, all.x = T, by = "ward")
abundance <- merge(abundance, D_, all.x = T, by = "ward")
abundance <- merge(abundance, W_, all.x = F, by = "ward") # Restrict inclusion to wards within HCL with known activity (n = 357)

# Wards with zero connectivity and/or consumption
# are excluded from the S_ and D_ dataframes by design. Their respective
# values should be zero but are replaced with NA during left full joint.

# Replace missing values with zeroes for connectivity and consumption
zeroCols  <- c(names(S_)[-1], names(D_)[-1])
abundance <- abundance[ , (zeroCols) := lapply(.SD, function(x) {x[is.na(x)] <- 0.0; return(x)}), .SDcols = zeroCols]

# Check that no NA remains
tmp <- lapply(abundance, function(x) stopifnot(all(!is.na(x))))

save(abundance, file = "abundance.Rdata")



