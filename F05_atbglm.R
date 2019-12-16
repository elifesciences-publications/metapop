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

load(file = "modeldata.Rdata") # Output from F03_modeldataprep.R

head(transf.dat)

# Identify variant names
taxonName<-unique(transf.dat$BacType)

# Store per-variant model results in a list
ResultList2<-list()

# Loop over variants and populate ResultList2
for(i in 1:length(taxonName)){
  selTaxon<-which(transf.dat$BacType == taxonName[i])
  
  # Run Poisson GLM
  metab = glm(N_patients ~ C_control + S_connectivity + ddd_carba + ddd_c1g_c2g
              + ddd_c3g_classic + ddd_c3g_pyo + ddd_glyco + ddd_oxa
              + ddd_fq + ddd_bsp + ddd_nsp, data = transf.dat,subset=selTaxon, family="poisson")
  
  # Model summary
  summary <- summary(metab)
  
  # Access coefficients, standard error, etc
  info <- summary$coefficients
  
  # Select beta coefficients for all variables in the model
  betas <- info[, 'Estimate']
  
  # Calculate the profile confidence intervals
  ci95 <- confint(metab)
  
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
ResultList2[[17]]$ci95

##################################################################################################
# Figure 3 of beta coefficients and confidence intervals for each antibiotic class in each model

# Extract the betas + CI's table for each element of the list 
beta.cisList2 <- lapply(ResultList2, function (x){unlist(x[['betas.cis']])})

# Set names for each list
names.betacisList2 <- unlist(lapply(ResultList2, function(x) x[c('name')]))

# Apply those names to the list
names(beta.cisList2) <- names.betacisList2

# Reorder the elements of the list
beta.cisList2.order <- beta.cisList2[c(1,2,3,4,5,6,7,8,9,10,11,16,17,12,13,14,15)]

# Make one list of all bacteria except Efaec and Abau, which have wider confidence intervals
beta.cisList2.a <- beta.cisList2.order[c(1:13)]

# Separate list for Abau and Efaec due to much larger confidence intervals
beta.cisList2.b <- beta.cisList2.order[c(14:17)]

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
# Generate figure

# Panel 1 - all taxa except Abau and Efae
svg(file = "atbmod_pane1.svg", 4, 8)
{
  showpanes <- function(yl) { 
    rect(0.75,  yl[1], 3.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(3.75,  yl[1], 6.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(6.75,  yl[1], 9.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(9.75,  yl[1], 11.25, yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(11.75, yl[1], 13.25, yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
  }
  
  par(mfrow = c(10,1))
  par(mar = c(0.2,4,0,4))
  
  p <- length(coef_carba.a)
  marker.cex <- 1.2
  yl <- c(-0.5, 0.5)  
  
  yl <- c(-0.65, 0.5)
  plot(coef_nsp.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "AMC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_nsp_li.a[ord], 1:p, coef_nsp_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_nsp.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-0.5, 0.8)
  plot(coef_bsp.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "TZP", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_bsp_li.a[ord], 1:p, coef_bsp_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_bsp.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-0.5, 0.5)
  plot(coef_fq.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "FQ", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_fq_li.a[ord], 1:p, coef_fq_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_fq.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-0.5, 0.5)
  plot(coef_oxa.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "OXA", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_oxa_li.a[ord], 1:p, coef_oxa_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_oxa.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-0.5, 0.5)
  plot(coef_glyco.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "VAN/TEC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_glyco_li.a[ord], 1:p, coef_glyco_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_glyco.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex) 
  
  yl <- c(-0.5, 0.5)
  plot(coef_c3g_pyo.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "CTZ/FEP", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c3g_pyo_li.a[ord], 1:p, coef_c3g_pyo_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c3g_pyo.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  yl <- c(-0.7, 0.8)
  plot(coef_c3g_classic.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "CTX/CRO", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c3g_classic_li.a[ord], 1:p, coef_c3g_classic_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c3g_classic.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex) 
  
  yl <- c(-0.3, 0.4)
  plot(coef_c1gc2g.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "1GC/2GC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c1gc2g_li.a[ord], 1:p, coef_c1gc2g_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c1gc2g.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex) 
  
  yl <- c(-0.4, 0.75)
  plot(coef_carba.a[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "IPM/MEM", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_carba_li.a[ord], 1:p, coef_carba_ui.a[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_carba.a[ord], pch = 19, col = bugcolors.a, cex = marker.cex)
  
  axis(1, at = 1:p, labels = buglabels.a[ord], las = 2)  
}
dev.off()

# Panel 2 - Abau and Efae
svg(file = "atbmod_pane2.svg", 1.9, 8)
{
  showpanes <- function(yl) {
    rect(0.75,  yl[1], 2.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(2.75,  yl[1], 4.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
  }
  
  par(mfrow = c(10,1))
  par(mar = c(0.2,4,0,4))
  
  p <- length(coef_carba.b)
  marker.cex <- 1.2
  yl <- c(-1, 1) * 1.25
  xl <- c(0.75, 4.25)
  
  plot(coef_nsp.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "AMC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_nsp_li.b[ord], 1:p, coef_nsp_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_nsp.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  plot(coef_bsp.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "TZP", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_bsp_li.b[ord], 1:p, coef_bsp_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_bsp.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  plot(coef_fq.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "FQ", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_fq_li.b[ord], 1:p, coef_fq_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_fq.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  plot(coef_oxa.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "OXA", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_oxa_li.b[ord], 1:p, coef_oxa_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_oxa.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  plot(coef_glyco.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "VAN/TEC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_glyco_li.b[ord], 1:p, coef_glyco_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_glyco.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex) 
  
  plot(coef_c3g_pyo.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "CTZ/FEP", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c3g_pyo_li.b[ord], 1:p, coef_c3g_pyo_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c3g_pyo.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  plot(coef_c3g_classic.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "CTX/CRO", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c3g_classic_li.b[ord], 1:p, coef_c3g_classic_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c3g_classic.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex) 
  
  plot(coef_c1gc2g.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "1GC/2GC", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_c1gc2g_li.b[ord], 1:p, coef_c1gc2g_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_c1gc2g.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex) 
  
  plot(coef_carba.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "IPM/MEM", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_carba_li.b[ord], 1:p, coef_carba_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_carba.b[ord], pch = 19, col = bugcolors.b, cex = marker.cex)
  
  axis(1, at = 1:p, labels = buglabels.b[ord], las = 2)  
}
dev.off()

#' #############################################################
#' 
#' Performs computations to examine possibly causal associations between use of 
#' each class of antibiotics and the number of infection episodes. Creates
#' Figure 3b.
#' 

# Extract GLM results
models <- lapply(ResultList2, function(x) x$betas.cis[-(1:3),])
names(models) <- unlist(lapply(ResultList2, function(x) x$name))

print(models)

# Define possibly causal associations between infection incidence
# and classes of antibiotics
possibly_causal <- list(
  ESCCOL_C3G_R = "ddd_c3g_classic",
  ESCCOL_CARBA_R = "ddd_carba",
  KLEPNE_C3G_R = "ddd_c3g_classic",
  KLEPNE_CARBA_R = "ddd_carba",
  ENTCLO_S = "ddd_nsp",
  ENTCLO_C3G_R = "ddd_c3g_classic",
  ENTCLO_CARBA_R = "ddd_carba",
  PSEAER_S = "ddd_c3g_classic",
  PSEAER_CARBA_R = "ddd_carba",
  ACIBAU_CARBA_S = "ddd_c3g_classic",
  ACIBAU_CARBA_R = "ddd_carba",
  ENCFAC_VANCO_S = c("ddd_c3g_classic", "ddd_c3g_pyo"),
  ENCFAC_VANCO_R = "ddd_glyco",
  STAAUR_OXA_R = c("ddd_c1g_c2g", "ddd_oxa", "ddd_nsp")
)

# Categorize antibiotic/variant associations as possibly causal or
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
  pcfilt <- which(dt$drug %in% possibly_causal[[bug]])
  dt[ pcfilt, pc := TRUE]
  return(dt)
}, simplify = F) %>% rbindlist


#Set up display labels for the boxplot (below)
#Possibly causal are labelled as "Poss. causal", all others labelled "Other"
pcmodel[, pc_display := c("Other", "Poss. causal")[pc + 1]][
  , pc_display := factor(pc_display, levels = c("Poss. causal", "Other"))
  ]

# Boxplot showing the difference in the magnitude of the beta coefficients
# for possibly causal associations vs. Other
{
  svg("possiblycausal_boxplot_2.svg", 4, 6)
  boxplot(betas ~ pc_display, pcmodel, ylim = c(-0.1, 0.4), outline = F, xlab = "", ylab = "Regression coefficient")
  beeswarm(betas ~ pc_display, pcmodel, method = "hex", add = T, pch = 19, col = c(rgb(0.8,0,0,.2), rgb(0,0,.8,.2)))
  dev.off()  
}


# Run wilcox.test on the beta coefficients to compare Possibly causal associations
# to others/
wilcox.test(betas ~ pc_display, pcmodel)

# Calculate the median value for beta coefficients for Possible Causal associations and Other
pcmodel[, .(med = median(betas)), by = pc]

# Compare coefficient significance for possibly causal associations and others
fisher.test(table(pcmodel$pc_display, pcmodel$`2.5 %` <= 0))

# Find resulting figures in files atbmod_pane1.svg, atbmod_pane2.svg and possiblycausal_boxplot_2.svg

#####################################################################################################

