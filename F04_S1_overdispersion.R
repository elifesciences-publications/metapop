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
#' Performs computations for global generalized linear models
#' and accompanying figures.
#' 
#' Tests for overdispersion in count models

library(data.table)
library(MASS)

load(file = "modeldata.Rdata") # Output from F03_modeldataprep.R

logliklist <- list()
taxonName <- unique(transf.dat.factor$BacType)

#' For each variant, fit a regular Poisson model, a quasi-Poisson model (which estimates under/overdispersion
#' through its theta parameter) and a negative binomial model. We keep track of the log-likelihoods of
#' the Poisson and negative binomial models. Remark that quasi models don't have a proper log-likelihood so we
#' cannot directly compare the quasi-Poisson and negative binomial models.
for(i in 1:length(taxonName)){
  print(taxonName[i])
  selTaxon<-which(transf.dat.factor$BacType == taxonName[i]) #Divide the data,
  # this variable will be called to identify each subset of the data
  
  mod_poisson <- glm(
    N_patients ~ C_control + S_connectivity + n_beds + ddd_total + PatStat, data = transf.dat.factor,subset=selTaxon,
    family=poisson)

  mod_quasi <- glm(
    N_patients ~ C_control + S_connectivity + n_beds + ddd_total + PatStat, data = transf.dat.factor,subset=selTaxon,
    family=quasipoisson)
  
  mod_negbin <- glm.nb(
    N_patients ~ C_control + S_connectivity + n_beds + ddd_total + PatStat, data = transf.dat.factor,subset=selTaxon)
  
  logLik(mod_negbin)
  
  summary_poisson <- summary(mod_poisson)
  summary_quasi <- summary(mod_quasi)
  summary_negbin <- summary(mod_negbin)
  
  logliklist[[i]] <- data.table(variant = taxonName[i],
                                poisson_loglik = logLik(mod_poisson),
                                negbin_loglik = logLik(mod_negbin),
                                quasipoisson_theta = summary_quasi$dispersion
                                )
}

logliktable <- rbindlist(logliklist)

print(logliktable)
