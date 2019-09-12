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
#' Prepares data for the downstream generalized linear models.
#' 
#' Analysis begins using the output dataframe "abundance" from F01_dataprep.R script
#' 
#' Output dataframe "transf.dat" will be used in downstream models.
#' 

library(data.table)
library(dplyr) # Tested with version 0.8.0.1 for R version 3.6.0
library(tidyr) #Tested with version 0.8.3 for R version 3.6.0

load("abundance.Rdata")

#Select and gather variables from abundance dataframe
#that are specific to each bacteria type: N_patients, Control, Connectivity

#Gather infection episode (N_) columns to place all bacteria in one column
#Key is "BacType" followed by "N_patients" = number of patients
N.dat <- abundance %>%
  select(ward, starts_with("N_", ignore.case = FALSE)) %>%
  gather(starts_with("N_") ,key="BacType", value="N_patients")

N.dat$BacType <- sub(pattern="N_", replacement = '', N.dat$BacType)

#Gather Control (C_) columns
C.dat <- abundance %>%
  select(ward, starts_with("C_", ignore.case = FALSE))%>%
  gather(starts_with("C_") ,key="BacType", value="C_control")

C.dat$BacType <- sub(pattern="C_", replacement = '', C.dat$BacType)

#Gather Connectivity (S_) columns
S.dat <- abundance %>%
  select(ward, starts_with("S_", ignore.case = FALSE))%>%
  gather(starts_with("S_") ,key="BacType", value="S_connectivity")

S.dat$BacType <- sub(pattern="S_", replacement = '', S.dat$BacType)
  

#### Select columns common to all bacteria types for each
## ward (antibiotic consumption, n_beds, PatStat)

ward.dat <- select(abundance, ward, starts_with("ddd_", ignore.case = FALSE), n_beds, PatStat)

# Join N patients, control, and connectivity dataframes by BOTH BacType AND ward
join1 <- full_join(N.dat, C.dat, by = c("BacType"="BacType","ward"="ward"))
join2 <- full_join(join1, S.dat, by= c("BacType"="BacType","ward"="ward"))

# Join this data frame to the ward data
mod.dat.raw <- full_join(join2, ward.dat, by="ward")


#Before running models, transform the response variables using a log-2 transformation

#Change data.frame to a data.table
mod.dat.raw.dt <- data.table(mod.dat.raw)

#Select the columns that will be transformed
logVars <- c("C_control", "S_connectivity","ddd_total","ddd_carba", "ddd_c1g_c2g",
             "ddd_c3g_classic","ddd_c3g_pyo","ddd_glyco","ddd_oxa","ddd_fq","ddd_bsp","ddd_nsp","n_beds")

#Transform the selected variables using data.table package
#To avoid infinity values from log-transformation, 
#first convert all 0 values to 1/2 the non-zero minimum value
#Then apply a log-2 transformation
transf.dat <- mod.dat.raw.dt [ , (logVars) := lapply(.SD, function(x) {
  xmin <- min(x[x > 0])
  x[x < xmin] <- xmin / 2
  return(log2(x))
}) , .SDcol = logVars]

save(transf.dat, mod.dat.raw, file = "modeldata.Rdata")
