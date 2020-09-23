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
#' Computations for global generalized linear models and specific antibiotic
#' models for all 3rd generation cephalosporin (3GCR) resistant taxa pooled together.
#' 

library(data.table) # Tested with version 1.12.2 for R version 3.6.0
library(dplyr)      # Tested with version 0.8.0.1 for R version 3.6.0
library(visreg)     #Tested with version 2.5-0 for R version 3.6.0

#Begin with dataframe "mod.dat.raw.factor" from F03 script
load("modeldata.Rdata")
head(mod.dat.raw.factor)

# List of 3GCR taxa
c3gR.list <- c("ESCCOL_C3G_R", "KLEPNE_C3G_R","KLEPNE_CARBA_R",
               "ENTCLO_C3G_R", "ENTCLO_CARBA_R","PSEAER_S","PSEAER_CARBA_R",
               "ACIBAU_CARBA_S", "ACIBAU_CARBA_R","ENCFAC_VANCO_S","ENCFAC_VANCO_R",
               "STAAUR_OXA_R")

# Subset the data to include only taxa in the c3gR.list 
c3gR.dat.raw <- subset(mod.dat.raw.factor, BacType %in% c3gR.list)

# Sum N_patients, C_control, and S_connectivity of different variants by ward
c3gR.dat.raw2 <- c3gR.dat.raw %>%
  group_by(ward) %>%
  summarize(C_controlC3G = sum(C_control),
            N_patsC3G = sum(N_patients),
            S_connC3G = sum(S_connectivity))

# Select ward-level variables
#######################################################################
# IMPORTANT NOTE: Detach package MASS (conflicts with dplyr) if code 
# below produces error
#######################################################################
ddd.dat.c3g <- select(c3gR.dat.raw, ward, n_beds, PatStat,  starts_with("ddd_"))
ddd.dat.c3g <- distinct(ddd.dat.c3g)

# Join the summed N, S, C data by taxa to the ward-level variables
c3g.R.dat.raw3 <- left_join(c3gR.dat.raw2, ddd.dat.c3g, by="ward")

# Change data.frame to a data.table
c3g.R.dat.raw3 <- data.table(c3g.R.dat.raw3)

# Select the columns that will be log2 transformed
logVars.c3g <- c("C_controlC3G", "S_connC3G","ddd_total","ddd_carba", "ddd_c1g_c2g",
                 "ddd_c3g_classic","ddd_c3g_pyo","ddd_glyco","ddd_oxa","ddd_fq","ddd_bsp","ddd_nsp", "ddd_amin", "ddd_amox", "n_beds")

# Transform the selected variables using data.table package
# To avoid infinity values from log-transformation, 
# first convert all 0 values to 1/2 the non-zero minimum value
# Then apply a log-2 transformation
c3g.dat <- c3g.R.dat.raw3[ , (logVars.c3g) := lapply(.SD, function(x) {
  xmin <- min(x[x > 0])
  x[x < xmin] <- xmin / 2
  return(log2(x))
}) , .SDcol = logVars.c3g]


# Global 3GCR model
c3gR.mod1 <- glm(N_patsC3G ~ C_controlC3G + n_beds + S_connC3G + ddd_total + PatStat, 
                 family = quasipoisson, data = c3g.dat)

# Run the 3GCR model using specific antibiotics
c3gR.mod2 <- glm(N_patsC3G ~ C_controlC3G + S_connC3G +
                   ddd_carba + ddd_c1g_c2g
                 + ddd_c3g_classic + ddd_c3g_pyo + ddd_glyco + ddd_oxa
                 + ddd_fq + ddd_bsp + ddd_nsp + ddd_amin + ddd_amox, 
                 family = quasipoisson, data = c3g.dat)

#' #############################################################
#' 
#' Computations for global generalized linear models and specific antibiotic
#' models for all Carbapenem resistant (CR) taxa pooled together.
#' 

#Begin with dataframe "mod.dat.raw.factor" from F03 script
head(mod.dat.raw.factor)

# List of CR taxa
carbaR.list <-c("ESCCOL_CARBA_R","KLEPNE_CARBA_R","ENTCLO_CARBA_R","PSEAER_CARBA_R","ACIBAU_CARBA_R",
                "ENCFAC_VANCO_S","ENCFAC_VANCO_R","STAAUR_OXA_R")

# Subset the data to include only taxa in the carbaR.list 
carbaR.dat.raw <- subset(mod.dat.raw.factor, BacType %in% carbaR.list)

# Sum N_patients, C_control, and S_connectivity of different variants by ward
carbaR.dat.raw2 <- carbaR.dat.raw %>%
  group_by(ward) %>%
  summarize(C_controlCarba = sum(C_control),
            N_patsCarba = sum(N_patients),
            S_connCarba = sum(S_connectivity))

# Select ward-level variables
ddd.dat.carba <- select(carbaR.dat.raw, ward, n_beds,PatStat,  starts_with("ddd_"))
ddd.dat.carba <- distinct(ddd.dat.carba)

# Join the summed N, S, C data by taxa to the ward-level variables
carbaR.dat.raw3 <- left_join(carbaR.dat.raw2, ddd.dat.carba, by="ward")

#Change data.frame to a data.table
carbaR.dat.raw3 <- data.table(carbaR.dat.raw3)

#Select the columns that will be transformed
logVars.carba <- c("C_controlCarba", "S_connCarba","ddd_total","ddd_carba", "ddd_c1g_c2g",
                   "ddd_c3g_classic","ddd_c3g_pyo","ddd_glyco","ddd_oxa","ddd_fq","ddd_bsp","ddd_nsp", "ddd_amin", "ddd_amox","n_beds")

# Transform the selected variables using data.table package
# Note: To avoid infinity values from log-transformation, 
# first convert all 0 values to 1/2 the non-zero minimum value
# Then apply a log-2 transformation
carbaR.dat <- carbaR.dat.raw3[ , (logVars.carba) := lapply(.SD, function(x) {
  xmin <- min(x[x > 0])
  x[x < xmin] <- xmin / 2
  return(log2(x))
}) , .SDcol = logVars.carba]


# Global CR model
carbaR.mod1 <- glm(N_patsCarba ~ C_controlCarba + n_beds + S_connCarba + ddd_total + PatStat, 
                   family = quasipoisson, data = carbaR.dat)

# CR model using specific antibiotics
carbaR.mod2 <- glm(N_patsCarba ~ C_controlCarba + S_connCarba +
                     ddd_carba + ddd_c1g_c2g
                   + ddd_c3g_classic + ddd_c3g_pyo + ddd_glyco + ddd_oxa
                   + ddd_fq + ddd_bsp + ddd_nsp + ddd_amin + ddd_amox , 
                   family = quasipoisson, data = carbaR.dat)

#' #############################################################
#' 
#' Creates Figure 3, showing the percent change represented by the 
#' beta coefficients for the global
#' model (Figure 3a) and graphs the relationship between the consumption of specific
#' antibiotics and infection incidence (Figure 3b)
#' 

# Combine 3GCR and CR models into a single list
# Global model list
combined.simplist <- list(c3gR.mod1,carbaR.mod1)

# Specific antibiotic model list
combined.atbslist <- list(c3gR.mod2,carbaR.mod2)

# Calculate 95% confidence interval for each model in each list
coef.combined.simp <- lapply(combined.simplist, function(x) {
  cis.comb <- confint(x)
  betas <- (exp(coefficients(x)) - 1)*100
  ci95 <- (exp(cis.comb) - 1)*100
  return(cbind(ci95, betas))
}) 

coef.combined.atbs <- lapply(combined.atbslist, function(x) {
  cis.comb <- confint(x)
  betas <- (exp(coefficients(x)) - 1)*100
  ci95 <- (exp(cis.comb) - 1)*100
  return(cbind(ci95, betas))
}) 

# Rename each object in the list
names(coef.combined.simp) <- c("C3GR","CarbaR")
names(coef.combined.atbs) <- c("C3GR","CarbaR")

###Coefficients #############################

# Prepare Figure 3a:
# Extract from each list (3GCR and CR) the beta coefficient and lower and upper confidence intervals
# for each variable
{
  # Patient fragility (ward type) Category 1
  coef_ps1 <- unlist(lapply(coef.combined.simp, function(x) x[grep("PatStat1", rownames(x)),3]))
  coef_ps_li1 <- unlist(lapply(coef.combined.simp, function(x) x[grep("PatStat1", rownames(x)),1]))
  coef_ps_ui1 <- unlist(lapply(coef.combined.simp, function(x) x[grep("PatStat1", rownames(x)),2]))
  
  # Patient fragility (ward type) Category 2
  coef_ps2 <- unlist(lapply(coef.combined.simp, function(x) x[grep("PatStat2", rownames(x)),3]))
  coef_ps_li2 <- unlist(lapply(coef.combined.simp, function(x) x[grep("PatStat2", rownames(x)),1]))
  coef_ps_ui2 <- unlist(lapply(coef.combined.simp, function(x) x[grep("PatStat2", rownames(x)),2]))
  
  # Number of beds (Ward size)
  coef_n <- unlist(lapply(coef.combined.simp, function(x) x[grep("n_beds", rownames(x)),3]))
  coef_n_li <- unlist(lapply(coef.combined.simp, function(x) x[grep("n_beds", rownames(x)),1]))
  coef_n_ui <- unlist(lapply(coef.combined.simp, function(x) x[grep("n_beds", rownames(x)),2]))
  
  # Connectivity
  coef_s <- unlist(lapply(coef.combined.simp, function(x) x[grep("S_", rownames(x)),3]))
  coef_s_li <- unlist(lapply(coef.combined.simp, function(x) x[grep("S_", rownames(x)),1]))
  coef_s_ui <- unlist(lapply(coef.combined.simp, function(x) x[grep("S_", rownames(x)),2]))
  
  # Antibiotic consumption
  coef_atb <- unlist(lapply(coef.combined.simp, function(x) x[grep("ddd_total", rownames(x)),3]))
  coef_atb_li <- unlist(lapply(coef.combined.simp, function(x) x[grep("ddd_total", rownames(x)),1]))
  coef_atb_ui <- unlist(lapply(coef.combined.simp, function(x) x[grep("ddd_total", rownames(x)),2])) 
}

# Assign labels and colors to 3GCR and CR models
buglabs <- c("All 3GCR", "All CR")
bugcols <- c( "deepskyblue", "darkmagenta")

ord <- T

errorbar_width <- 0.04

# Plot the beta coefficients and 95% confidence intervals for each variable in each model
svg(file = "glm_pooled_pane1.svg", 1.5, 6)
{
  p <- length(coef.combined.simp)
  
  par(mfrow = c(6,1))
  par(mar = c(1,4,1,4))
  xl <- c(0.75, 2.25)
  marker.cex <- 1.25

  yl <- c(-40,50)
  plot(coef_ps2[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "Intensive care", bty = "n", type = "n")
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_ps_li2[ord], 1:p, coef_ps_ui2[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_ps2[ord], pch = 19, col = bugcols, cex = marker.cex)
  
  yl <- c(-20,100)
  plot(coef_ps1[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "Progressive care", bty = "n", type = "n")
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_ps_li1[ord], 1:p, coef_ps_ui1[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_ps1[ord], pch = 19, col = bugcols, cex = marker.cex)
  
  yl <- c(-15,15)
  plot(coef_n[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "Ward size", bty = "n", type = "n")
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_n_li[ord], 1:p, coef_n_ui[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_n[ord], pch = 19, col = bugcols, cex = marker.cex)
  
  yl <- c(-10,15)
  plot(coef_s[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "Connectivity", bty = "n", type = "n")
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_s_li[ord], 1:p, coef_s_ui[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_s[ord], pch = 19, col = bugcols, cex = marker.cex)
  
  yl <- c(-5,15)
  plot(coef_atb[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "Antibiotic use", bty = "n", type = "n")
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_atb_li[ord], 1:p, coef_atb_ui[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_atb[ord], pch = 19, col = bugcols, cex = marker.cex)
  
  axis(1, at = 1:p, labels = buglabs[ord], las = 2)  
}
dev.off()


# Use visreg package to show the response curve of number of infection episodes 
# to consumption of specific antibiotic classes 
# First line = 3GCR, second line = CR, columns = 3GC, Carba, TZP (Figure 3a and 3b)
# Third line = histogram of 3GC, Carba, TZP use (Figure 3c)

#Print the percentage change represented by beta coefficients and confidence interval for each
prettyConfint <- function(response, drug) {
  beta <- sprintf("β = %.1f (%.1f, %.1f)", coef.combined.atbs[[response]][drug, 3], coef.combined.atbs[[response]][drug, 1], coef.combined.atbs[[response]][drug, 2])
  
} 
# NOTE - due to encoding the beta symbol above may be changed
# to another symbol when saving / reopening file. Copy-paste beta symbol to correct.

svg(file = "visreg_pooled.svg", 6, 6)
{
  par(mfrow = c(3, 3))
  par(mar = c(4, 4, 1, 1))
  coef_legend_cex = 0.9
  
  visreg(c3gR.mod2, "ddd_c3g_classic", scale=c("response"),xlab="CTX/CRO use (ddd/bed/y)",ylab="No. 3GCR episodes/ward/y",xlim = c(-6,6)*1.1, ylim=c(0,12), rug = F, xaxt = "n")
  xseq <- seq(-2,2,1)
  axis(1, at = xseq/log(2)*log(10), labels = 10^xseq)
  legend("topleft", title = prettyConfint("C3GR", "ddd_c3g_classic"), legend = "", bty = "n", cex = coef_legend_cex)
  
  visreg(c3gR.mod2, "ddd_carba", scale=c("response"),xlab="IPM/MEM use (ddd/bed/y)",ylab="No. 3GCR episodes/ward/y",xlim = c(-6,6)*1.1,ylim=c(0,12), rug = F, xaxt = "n")
  xseq <- seq(-2,2,1)
  axis(1, at = xseq/log(2)*log(10), labels = 10^xseq)
  legend("topleft", title = prettyConfint("C3GR", "ddd_carba"), legend = "", bty = "n", cex = coef_legend_cex)
  
  visreg(c3gR.mod2, "ddd_bsp", scale=c("response"),xlab="TZP use (ddd/bed/y)",ylab="No. 3GCR episodes/ward/y",xlim = c(-6,6)*1.1,ylim=c(0,12), rug = F, xaxt = "n")
  xseq <- seq(-2,2,1)
  axis(1, at = xseq/log(2)*log(10), labels = 10^xseq)
  legend("topleft", title = prettyConfint("C3GR", "ddd_bsp"), legend = "", bty = "n", cex = coef_legend_cex)
  
  visreg(carbaR.mod2, "ddd_c3g_classic", scale=c("response"),xlab="CTX/CRO use (ddd/bed/y)",ylab="No. CR episodes/ward/y",xlim = c(-6,6)*1.1, ylim=c(0,5), rug = F, xaxt = "n")
  xseq <- seq(-2,2,1)
  axis(1, at = xseq/log(2)*log(10), labels = 10^xseq)
  legend("topleft", title = prettyConfint("CarbaR", "ddd_c3g_classic"), legend = "", bty = "n", cex = coef_legend_cex)
  
  
  visreg(carbaR.mod2, "ddd_carba", scale=c("response"),xlab="IPM/MEM use (ddd/bed/y)",ylab="No. CR episodes/ward/y", xlim = c(-6,6)*1.1, ylim=c(0,5), rug = F, xaxt = "n")
  xseq <- seq(-2,2,1)
  axis(1, at = xseq/log(2)*log(10), labels = 10^xseq)
  legend("topleft", title = prettyConfint("CarbaR", "ddd_carba"), legend = "", bty = "n", cex = coef_legend_cex)
  
  visreg(carbaR.mod2, "ddd_bsp", scale=c("response"),xlab="TZP use (ddd/bed/y)",ylab="No. CR episodes/ward/y",xlim = c(-6,6)*1.1,ylim=c(0,5), rug = F, xaxt = "n")
  xseq <- seq(-2,2,1)
  axis(1, at = xseq/log(2)*log(10), labels = 10^xseq)
  legend("topleft", title = prettyConfint("CarbaR", "ddd_bsp"), legend = "", bty = "n", cex = coef_legend_cex)
  

  hist_xlim <- c(log10(0.01), log10(100))
  hist_col <- "azure2"
  
  hist_c3g_classic <- c3g.dat$ddd_c3g_classic / log(10) * log(2)
  hist_c3g_classic[hist_c3g_classic == min(hist_c3g_classic)] <- hist_xlim[1]
  
  hist_carba <- c3g.dat$ddd_carba / log(10) * log(2)
  hist_carba[hist_carba == min(hist_carba)] <- hist_xlim[1]
  
  hist_bsp <- c3g.dat$ddd_bsp / log(10) * log(2)
  hist_bsp[hist_bsp == min(hist_bsp)] <- hist_xlim[1]
  
  hist(hist_c3g_classic,xlab="CTX/CRO use (ddd/bed/y)",main="",
       freq=TRUE,xlim = hist_xlim, ylim=c(0,180),breaks=12,col=hist_col, xaxt="n")
  xseq <- seq(-3,2,1)
  axis(1, at = xseq, labels = 10^xseq)
  
  hist(hist_carba,xlab="IPM/MEM use (ddd/bed/y)",main="",
       freq=TRUE,xlim = hist_xlim, ylim=c(0,180),breaks=12,col=hist_col, xaxt="n")
  xseq <- seq(-3,2,1)
  axis(1, at = xseq, labels = 10^xseq)
  
  hist(hist_bsp,xlab="TZP use (ddd/bed/y)",main="",
       freq=TRUE,xlim = hist_xlim, ylim=c(0,180),breaks=12,col=hist_col, xaxt="n")
  xseq <- seq(-3,2,1)
  axis(1, at = xseq, labels = 10^xseq)
  
  
}
dev.off()

# Find resulting figures in files glm_pooled_pane1.svg and visreg_pooled.svg

#####################################################################
# Calculate the model coefficients for pooled 3GCR and Carba-R 
# infections:

#Begin with dataframe "mod.dat.raw.factor" from F03 script
load("modeldata.Rdata")
head(mod.dat.raw.factor)

# Make a list of all taxa that are resistant to either 3GC or carbapenems
c3g.carb.list <- c("ESCCOL_C3G_R","ESCCOL_CARBA_R", "KLEPNE_C3G_R",
                   "KLEPNE_CARBA_R","ENTCLO_C3G_R", "ENTCLO_CARBA_R",
                   "PSEAER_S","PSEAER_CARBA_R","ACIBAU_CARBA_S", 
                   "ACIBAU_CARBA_R","ENCFAC_VANCO_S","ENCFAC_VANCO_R",
                "STAAUR_OXA_R")

# Subset the data to include only taxa in the c3g.carb.list 
c3gR.carb.dat.raw <- subset(mod.dat.raw.factor, BacType %in% c3g.carb.list)

# Sum N_patients, C_control, and S_connectivity of different variants by ward
c3gR.carb.dat.raw2 <- c3gR.carb.dat.raw %>%
  group_by(ward) %>%
  summarize(C_controlC3G.Car = sum(C_control),
            N_patsC3G.Car = sum(N_patients),
            S_connC3G.Car = sum(S_connectivity))

# Select ward-level variables
ddd.dat.c3g.carba <- select(c3gR.dat.raw, ward, n_beds, PatStat,  starts_with("ddd_"))
ddd.dat.c3g.carba <- distinct(ddd.dat.c3g.carba)

# Join the summed N, S, C data by taxa to the ward-level variables
c3gR.carb.dat.raw3 <- left_join(c3gR.carb.dat.raw2, ddd.dat.c3g.carba, by="ward")

# Change data.frame to a data.table
c3gR.carb.dat.raw3 <- data.table(c3gR.carb.dat.raw3)

# Select the columns that will be log2 transformed
logVars.c3g.carba <- c("C_controlC3G.Car", "S_connC3G.Car",
                       "ddd_total","ddd_carba", "ddd_c1g_c2g",
                 "ddd_c3g_classic","ddd_c3g_pyo","ddd_glyco",
                 "ddd_oxa","ddd_fq","ddd_bsp","ddd_nsp", "ddd_amin", "ddd_amox","n_beds")

# Transform the selected variables using data.table package
# To avoid infinity values from log-transformation, 
# first convert all 0 values to 1/2 the non-zero minimum value
# Then apply a log-2 transformation
c3g.carba.dat <- c3gR.carb.dat.raw3[ , (logVars.c3g.carba) := lapply(.SD, function(x) {
  xmin <- min(x[x > 0])
  x[x < xmin] <- xmin / 2
  return(log2(x))
}) , .SDcol = logVars.c3g.carba]


# Run the pooled 3GCR - Carba-R model using specific antibiotics
c3gR.carba.mod2 <- glm(N_patsC3G.Car ~ C_controlC3G.Car + S_connC3G.Car +
                   ddd_carba + ddd_c1g_c2g
                 + ddd_c3g_classic + ddd_c3g_pyo + ddd_glyco + ddd_oxa
                 + ddd_fq + ddd_bsp + ddd_nsp + ddd_amin + ddd_amox, 
                 family = quasipoisson, data = c3g.carba.dat)

summary(c3gR.carba.mod2)

# Calculate the confidence intervals for each variable in the model:
confint(c3gR.carba.mod2)

#Convert raw coefficients and confidence intervals to percentages
c3gR.carba.betas <- (exp(coefficients(c3gR.carba.mod2)) - 1)*100
c3gR.carba.mod2.ci95 <- (exp(confint(c3gR.carba.mod2)) - 1)*100

# Bind confidence interval to beta coefficient from the model 
# (to improve readability)
cbind(c3gR.carba.betas, c3gR.carba.mod2.ci95)

#summary(step(c3gR.carba.mod2))
