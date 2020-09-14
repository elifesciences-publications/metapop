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
#' Computations for generalized linear models using specific
#' classes of antibiotics and accompanying figures.
#' 

library(tidyverse) #Tested with version 1.2.1 for R version 3.6.0
library(beeswarm)  #Tested with version 0.2.3 for R version 3.6.0
library(data.table)
library(MASS)

load(file = "modeldata.Rdata") # Output from F03_modeldataprep.R

head(transf.dat.factor)

# Identify variant names
taxonName<-unique(transf.dat.factor$BacType)

# Store per-variant model results in a list
ResultList2<-list()

# Loop over variants and populate ResultList2
for(i in 1:length(taxonName)){
  selTaxon<-which(transf.dat.factor$BacType == taxonName[i])
  
  # Fit quasi-Poisson GLM
  metab = glm(N_patients ~ C_control + S_connectivity + ddd_carba + ddd_c1g_c2g
              + ddd_c3g_classic + ddd_c3g_pyo + ddd_glyco + ddd_oxa
              + ddd_fq + ddd_bsp + ddd_nsp + ddd_amin + ddd_amox, 
              data = transf.dat.factor,subset=selTaxon, family = quasipoisson)

  # Model summary
  summary <- summary(metab)
  
  # Access coefficients, standard error, etc
  info <- summary$coefficients
  
  # Select beta coefficients for all variables in the model
  betas1 <- info[, 'Estimate']
  
  # Calculate the profile confidence intervals
  ci95.1 <- confint(metab)
  
  
  #Convert raw coefficients and confidence intervals to percentages
  betas <- (exp(betas1) - 1)*100
  ci95 <- (exp(ci95.1) - 1)*100
  
  # Concatenate beta estimates and 95% CIs
  betas.cis <- cbind(betas,ci95)
  
  # Gather results
  resList<-list(
    name = paste(taxonName[i]),
    metab = metab,
    summary = summary,
    betas=betas,
    ci95=ci95,
    betas.cis=betas.cis
  )
  
  # Append results to the model list
  ResultList2[[length(ResultList2)+1]]<-resList
  
}

# Use the [[]] operator to access individual results
ResultList2[[17]]$name
ResultList2[[17]]$betas.cis

##################################################################################################
# Figure 2 of percent change represented by beta coefficients 
# and confidence intervals for each antibiotic class in each model

# Extract the betas + CI's table for each element of the list 
beta.cisList2 <- lapply(ResultList2, function (x){unlist(x[['betas.cis']])})

# Set names for each list
names.betacisList2 <- unlist(lapply(ResultList2, function(x) x[c('name')]))

# Apply those names to the list
names(beta.cisList2) <- names.betacisList2

# Make a separate list for Abau and Efaec due to much larger confidence intervals
beta.cisList2.a <- beta.cisList2[c(1,2,3,4,5,6,7,8,9,10,11,16,17)]

beta.cisList2.b <- beta.cisList2[c(12,13,14,15)]

# Identify the beta coefficient, lower, and upper confidence intervals for each variable
# and each variant
{
  coef_carba.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_carba", rownames(x)),1]))
  coef_carba_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_carba", rownames(x)),2]))
  coef_carba_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_carba", rownames(x)),3]))
  
  coef_c1gc2g.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_c1g_c2g", rownames(x)),1]))
  coef_c1gc2g_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_c1g_c2g", rownames(x)),2]))
  coef_c1gc2g_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_c1g_c2g", rownames(x)),3]))
  
  coef_c3g_classic.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_c3g_classic", rownames(x)),1]))
  coef_c3g_classic_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_c3g_classic", rownames(x)),2]))
  coef_c3g_classic_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_c3g_classic", rownames(x)),3]))
  
  coef_c3g_pyo.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_c3g_pyo", rownames(x)),1]))
  coef_c3g_pyo_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_c3g_pyo", rownames(x)),2]))
  coef_c3g_pyo_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_c3g_pyo", rownames(x)),3]))
  
  coef_glyco.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_glyco", rownames(x)),1]))
  coef_glyco_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_glyco", rownames(x)),2]))
  coef_glyco_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_glyco", rownames(x)),3]))
  
  coef_oxa.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_oxa", rownames(x)),1]))
  coef_oxa_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_oxa", rownames(x)),2]))
  coef_oxa_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_oxa", rownames(x)),3]))
  
  coef_fq.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_fq", rownames(x)),1]))
  coef_fq_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_fq", rownames(x)),2]))
  coef_fq_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_fq", rownames(x)),3]))
  
  coef_bsp.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_bsp", rownames(x)),1]))
  coef_bsp_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_bsp", rownames(x)),2]))
  coef_bsp_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_bsp", rownames(x)),3]))
  
  coef_nsp.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_nsp", rownames(x)),1]))
  coef_nsp_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_nsp", rownames(x)),2]))
  coef_nsp_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_nsp", rownames(x)),3]))  
  
  coef_amin.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_amin", rownames(x)),1]))
  coef_amin_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_amin", rownames(x)),2]))
  coef_amin_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_amin", rownames(x)),3]))  
  
  coef_amox.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_amox", rownames(x)),1]))
  coef_amox_li.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_amox", rownames(x)),2]))
  coef_amox_ui.a <- unlist(lapply(beta.cisList2.a, function(x) x[grep("ddd_amox", rownames(x)),3]))  
}

# Repeat for Abau and Efae
{
  coef_carba.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_carba", rownames(x)),1]))
  coef_carba_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_carba", rownames(x)),2]))
  coef_carba_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_carba", rownames(x)),3]))
  
  coef_c1gc2g.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_c1g_c2g", rownames(x)),1]))
  coef_c1gc2g_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_c1g_c2g", rownames(x)),2]))
  coef_c1gc2g_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_c1g_c2g", rownames(x)),3]))
  
  coef_c3g_classic.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_c3g_classic", rownames(x)),1]))
  coef_c3g_classic_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_c3g_classic", rownames(x)),2]))
  coef_c3g_classic_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_c3g_classic", rownames(x)),3]))
  
  coef_c3g_pyo.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_c3g_pyo", rownames(x)),1]))
  coef_c3g_pyo_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_c3g_pyo", rownames(x)),2]))
  coef_c3g_pyo_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_c3g_pyo", rownames(x)),3]))
  
  coef_glyco.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_glyco", rownames(x)),1]))
  coef_glyco_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_glyco", rownames(x)),2]))
  coef_glyco_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_glyco", rownames(x)),3]))
  
  coef_oxa.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_oxa", rownames(x)),1]))
  coef_oxa_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_oxa", rownames(x)),2]))
  coef_oxa_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_oxa", rownames(x)),3]))
  
  coef_fq.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_fq", rownames(x)),1]))
  coef_fq_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_fq", rownames(x)),2]))
  coef_fq_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_fq", rownames(x)),3]))
  
  coef_bsp.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_bsp", rownames(x)),1]))
  coef_bsp_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_bsp", rownames(x)),2]))
  coef_bsp_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_bsp", rownames(x)),3]))
  
  coef_nsp.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_nsp", rownames(x)),1]))
  coef_nsp_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_nsp", rownames(x)),2]))
  coef_nsp_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_nsp", rownames(x)),3]))  
  
  coef_amin.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_amin", rownames(x)),1]))
  coef_amin_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_amin", rownames(x)),2]))
  coef_amin_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_amin", rownames(x)),3]))  
  
  coef_amox.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_amox", rownames(x)),1]))
  coef_amox_li.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_amox", rownames(x)),2]))
  coef_amox_ui.b <- unlist(lapply(beta.cisList2.b, function(x) x[grep("ddd_amox", rownames(x)),3]))  
}

###########################################################################################
# Prepare the features of the figure

# Set variant labels 
buglabels.a <- c("EC", "3GCREC", "CREC", "KP", "3GCRKP", "CRKP", "EB", "3GCREB",
                 "CREB", "PA", "CRPA", "SA", "MRSA")

# Labels for Abau and Efae
buglabels.b <-c("AB", "CRAB", "EF", "VREF")

# Set colors for the beta coefficient points for each taxon
bugcolors.a <- c("violetred1", "violetred3", "violetred4",
                 "aquamarine", "aquamarine3", "aquamarine4",
                 "goldenrod1", "goldenrod3", "goldenrod4",
                 "turquoise2", "turquoise4"
                 ,"lightpink", "lightpink3")

# Set colors for the beta coefficient points for Abau and Efae
bugcolors.b <- c("red2", "red4",
                 "olivedrab3", "olivedrab4")

ord <- T

errorbar_width <- 0.04

###################################################################
# Generate Figure 2 

# Panel 1 (Figure 2a) - all taxa except Abau and Efae
svg(file = "atbmod_pane1.svg", 4, 10)
{
  showpanes <- function(yl) { 
    rect(0.75,  yl[1], 3.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(3.75,  yl[1], 6.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(6.75,  yl[1], 9.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(9.75,  yl[1], 11.25, yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(11.75, yl[1], 13.25, yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
  }
  
  par(mfrow = c(12,1))
  par(mar = c(0.4,4,0,4))
  
  p <- length(coef_carba.a)
  marker.cex <- 1.2
  yl <- c(-0.5, 0.5)  
  
  yl <- c(-5, 5)
  plot(coef_amin.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "AMIN", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_amin_li.a[ord], 1:p, coef_amin_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_amin.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-30, 45)
  plot(coef_fq.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "FQ", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_fq_li.a[ord], 1:p, coef_fq_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_fq.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-30, 35)
  plot(coef_glyco.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "VAN/TEC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_glyco_li.a[ord], 1:p, coef_glyco_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_glyco.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex) 
  
  yl <- c(-20, 15)
  plot(coef_oxa.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "OXA", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_oxa_li.a[ord], 1:p, coef_oxa_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_oxa.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-25, 90)
  plot(coef_carba.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "IPM/MEM", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_carba_li.a[ord], 1:p, coef_carba_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_carba.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-10, 110)
  plot(coef_bsp.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "TZP", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_bsp_li.a[ord], 1:p, coef_bsp_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_bsp.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-40, 50)
  plot(coef_c3g_pyo.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "CTZ/FEP", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c3g_pyo_li.a[ord], 1:p, coef_c3g_pyo_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c3g_pyo.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-50, 90)
  plot(coef_c3g_classic.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "CTX/CRO", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c3g_classic_li.a[ord], 1:p, coef_c3g_classic_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c3g_classic.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex) 
  
  yl <- c(-20, 25)
  plot(coef_c1gc2g.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "1GC/2GC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c1gc2g_li.a[ord], 1:p, coef_c1gc2g_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c1gc2g.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex) 
  
  yl <- c(-50, 30)
  plot(coef_nsp.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "AMC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_nsp_li.a[ord], 1:p, coef_nsp_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_nsp.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-1.5, 1.5)
  plot(coef_amox.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "AMX", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_amox_li.a[ord], 1:p, coef_amox_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_amox.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  axis(1, at = 1:p, labels = buglabels.a[ord], las = 2)  
}
dev.off()

# Panel 2 (Figure 2b) - Abau and Efae
svg(file = "atbmod_pane2.svg", 1.9, 10)
{
  showpanes <- function(yl) {
    rect(0.75,  yl[1], 2.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(2.75,  yl[1], 4.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
  }
  
  par(mfrow = c(12,1))
  par(mar = c(0.4,4,0,4))
  
  p <- length(coef_carba.b)
  marker.cex <- 1.2
  yl <- c(-1, 1) * 1.25
  xl <- c(0.75, 4.25)
  
  yl <- c(-30, 10)
  plot(coef_amin.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "AMIN", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_amin_li.b[ord], 1:p, coef_amin_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_amin.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  yl <- c(-65, 25)
  plot(coef_fq.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "FQ", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_fq_li.b[ord], 1:p, coef_fq_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_fq.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)

  yl <- c(-55, 150)
  plot(coef_glyco.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "VAN/TEC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_glyco_li.b[ord], 1:p, coef_glyco_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_glyco.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex) 
  
  yl <- c(-25, 70)
  plot(coef_oxa.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "OXA", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_oxa_li.b[ord], 1:p, coef_oxa_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_oxa.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  yl <- c(-20, 150)
  plot(coef_carba.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "IPM/MEM", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_carba_li.b[ord], 1:p, coef_carba_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_carba.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  yl <- c(-40, 200)
  plot(coef_bsp.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "TZP", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_bsp_li.b[ord], 1:p, coef_bsp_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_bsp.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  yl <- c(-30, 150)
  plot(coef_c3g_pyo.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "CTZ/FEP", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c3g_pyo_li.b[ord], 1:p, coef_c3g_pyo_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c3g_pyo.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  yl <- c(-70, 160)
  plot(coef_c3g_classic.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "CTX/CRO", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c3g_classic_li.b[ord], 1:p, coef_c3g_classic_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c3g_classic.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex) 
  
  yl <- c(-30, 90)
  plot(coef_c1gc2g.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "1GC/2GC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c1gc2g_li.b[ord], 1:p, coef_c1gc2g_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c1gc2g.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex) 
  
  yl <- c(-40, 90)
  plot(coef_nsp.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "AMC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_nsp_li.b[ord], 1:p, coef_nsp_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_nsp.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  yl <- c(-6, 2)
  plot(coef_amox.b[ord], ylim = yl, xlim = xl, xaxt = "n", xlab = "", ylab = "AMX", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_amox_li.b[ord], 1:p, coef_amox_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_amox.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  axis(1, at = 1:p, labels = buglabels.b[ord], las = 2)  
}
dev.off()

#' #############################################################
#' 
# Extract GLM results
models <- lapply(ResultList2, function(x) x$betas.cis[-(1:3),])
names(models) <- unlist(lapply(ResultList2, function(x) x$name))

print(models)

stop()

# Define possibly selective associations between infection incidence
# and classes of antibiotics
possibly_selective <- list(
  ESCCOL_C3G_R = "ddd_c3g_classic",
  ESCCOL_CARBA_R = "ddd_carba",
  KLEPNE_S = "ddd_amox",
  KLEPNE_C3G_R = "ddd_c3g_classic",
  KLEPNE_CARBA_R = "ddd_carba",
  ENTCLO_S = "ddd_nsp",
  ENTCLO_C3G_R = c("ddd_c3g_classic", "ddd_bsp"),
  ENTCLO_CARBA_R = "ddd_carba",
  PSEAER_S = "ddd_c3g_classic",
  PSEAER_CARBA_R = "ddd_carba",
  ACIBAU_CARBA_S = "ddd_c3g_classic",
  ACIBAU_CARBA_R = "ddd_carba",
  ENCFAC_VANCO_S = c("ddd_c3g_classic"),
  ENCFAC_VANCO_R = "ddd_glyco",
  STAAUR_OXA_S = "ddd_amox",
  STAAUR_OXA_R = c("ddd_c1g_c2g", "ddd_oxa", "ddd_nsp")
)

# Categorize antibiotic/variant associations as possibly selective or
# not, with its beta coefficient and lower and upper confidence intervals from
# the antibiotic models in ResultsList2. 
pcmodel <- sapply(names(models), function(bug) {
  mat <- models[[bug]]
  dt <- data.table(mat)[
    , drug := rownames(mat)
    ][
      , bug := bug
      ][
        , pc := FALSE
        ]
  pcfilt <- which(dt$drug %in% possibly_selective[[bug]])
  dt[ pcfilt, pc := TRUE]
  return(dt)
}, simplify = F) %>% rbindlist


#Set up display labels for the boxplot (below)
#Possibly selective are labelled as "Poss. selective", all others labelled "Other"
pcmodel[, pc_display := c("Other", "Poss. selective")[pc + 1]][
  , pc_display := factor(pc_display, levels = c("Poss. selective", "Other"))
  ]

# Boxplot showing the difference in the magnitude of the beta coefficients
# for possibly selective associations vs. Other
{
  svg("possiblyselective_boxplot_2.svg", 4, 6)
  boxplot(betas ~ pc_display, pcmodel, ylim = c(-20, 40), outline = F, xlab = "", ylab = "Percent change in incidence")
  beeswarm(betas ~ pc_display, pcmodel, method = "hex", add = T, pch = 19, col = c(rgb(0.8,0,0,.2), rgb(0,0,.8,.2)))
  dev.off()  
}

# Means and difference of the means of coefficients in possibly selective vs other
# associations
t.test(betas ~ pc_display, pcmodel)

# Wilcoxon test for boxplot figure significance
wilcox.test(betas ~ pc_display, pcmodel)

# Find resulting figures in files atbmod_pane1.svg, atbmod_pane2.svg and possiblyselective_boxplot_2.svg

#####################################################################################################

