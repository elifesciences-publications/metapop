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
#' 

library(data.table)

load(file = "modeldata.Rdata") # Output from F03_modeldataprep.R

#Begin with "transf.dat.factor" from F03 script.
head(transf.dat.factor)

#First identify each unit for the for-loop
taxonName <- unique(transf.dat.factor$BacType)

# Store per-variant model results in a list
ResultList<-list()

# Loop over variants and populate ResultList
for(i in 1:length(taxonName)){
  selTaxon<-which(transf.dat.factor$BacType == taxonName[i]) #Divide the data,
  # this variable will be called to identify each subset of the data
  
  ###### FACTOR CASE: not enough observations, exclude
  

  if(taxonName[i] %in% c("ACIBAU_CARBA_R", "ENCFAC_VANCO_R", "ESCCOL_CARBA_R")) {
    metab = glm(N_patients ~ C_control + S_connectivity + n_beds + ddd_total, data = transf.dat.factor,subset=selTaxon, family=quasipoisson)
    
  } else {# Run Poisson GLM
    metab = glm(N_patients ~ C_control + S_connectivity + n_beds + ddd_total + PatStat, data = transf.dat.factor,subset=selTaxon, family=quasipoisson)
    
  }
  
  # Model summary
  summary <- summary(metab)
  
  # Access coefficients, standard error, etc
  info <- summary$coefficients
  
  # Select beta coefficients for all variables in the model
  betas <- info[, 'Estimate']
  
  # Calculate the profile confidence intervals
  ci95 <- confint(metab)
  
  #Convert raw coefficients and confidence intervals to percentages
  percents.beta <- (exp(betas) - 1)*100
  percents.ci <- (exp(ci95) - 1)*100
  
  # Concatenate beta estimates and 95% CIs
  betas.cis <- cbind(percents.beta,percents.ci)
  
  ######## FACTOR CASE: CREC, CRAB, VREF; preserve order of coefficients
  if(taxonName[i] %in% c("ENCFAC_VANCO_R", "ESCCOL_CARBA_R", "ACIBAU_CARBA_R")) {
    betas <- c(head(betas, 5), PatStat1 = NA, PatStat2 = NA)
    ci95 <- rbind(head(ci95, 5), PatStat1 = c(NA, NA), PatStat2 = c(NA, NA))
    betas.cis <- rbind(head(betas.cis, 5), PatStat1 = c(NA, NA, NA), PatStat2 = c(NA, NA, NA))
  }
  
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
  ResultList[[length(ResultList)+1]]<-resList  
}

# Use the [[]] operator to access individual results
ResultList[[1]]$name
ResultList[[1]]$summary


##################################################################################################
# Figure 1 of percent change represented by beta coefficients 
# and confidence intervals for each variable in each model

# Extract the betas + CI's table for each element of the list
beta.cisList <- lapply(ResultList, function (x){unlist(x[['betas.cis']])})
# Set names for each list
names.betacisList <- unlist(lapply(ResultList, function(x) x[c('name')]))
# Apply those names to the list
names(beta.cisList) <- names.betacisList

# Make a separate list for Abau and Efaec due to much larger confidence intervals
beta.cisList.a <- beta.cisList[c(1,2,3,4,5,6,7,8,9,10,11,16,17)]

beta.cisList.b <- beta.cisList[c(12,13,14,15)]

# Identify the beta coefficient, lower, and upper confidence intervals for each variable
# and each variant
{
  #Antibiotic consumption
  coef_atb    <- unlist(lapply(beta.cisList.a, function(x) x[grep("ddd_total", rownames(x)),1]))
  coef_atb_li <- unlist(lapply(beta.cisList.a, function(x) x[grep("ddd_total", rownames(x)),2]))
  coef_atb_ui <- unlist(lapply(beta.cisList.a, function(x) x[grep("ddd_total", rownames(x)),3]))
  
  #Connectivity
  coef_s    <- unlist(lapply(beta.cisList.a, function(x) x[grep("S_connectivity", rownames(x)),1]))
  coef_s_li <-  unlist(lapply(beta.cisList.a, function(x) x[grep("S_connectivity", rownames(x)),2]))
  coef_s_ui <-  unlist(lapply(beta.cisList.a, function(x) x[grep("S_connectivity", rownames(x)),3]))
  
  #Number of beds 
  coef_n    <- unlist(lapply(beta.cisList.a, function(x) x[grep("n_beds", rownames(x)),1]))
  coef_n_li <- unlist(lapply(beta.cisList.a, function(x) x[grep("n_beds", rownames(x)),2]))
  coef_n_ui <- unlist(lapply(beta.cisList.a, function(x) x[grep("n_beds", rownames(x)),3]))
  
  #Patient fragility (ward type) Category 1
  coef_ps1    <- unlist(lapply(beta.cisList.a, function(x) x[grep("PatStat1", rownames(x)),1]))
  coef_ps_li1 <-  unlist(lapply(beta.cisList.a, function(x) x[grep("PatStat1", rownames(x)),2]))
  coef_ps_ui1 <-  unlist(lapply(beta.cisList.a, function(x) x[grep("PatStat1", rownames(x)),3]))
  
  #Patient fragility (ward type) Category 2
  coef_ps2    <- unlist(lapply(beta.cisList.a, function(x) x[grep("PatStat2", rownames(x)),1]))
  coef_ps_li2 <-  unlist(lapply(beta.cisList.a, function(x) x[grep("PatStat2", rownames(x)),2]))
  coef_ps_ui2 <-  unlist(lapply(beta.cisList.a, function(x) x[grep("PatStat2", rownames(x)),3]))
  
}

# Repeat for Abau and Efae
{
  #Antibiotic consumption
  coef_atb.b    <- unlist(lapply(beta.cisList.b, function(x) x[grep("ddd_total", rownames(x)),1]))
  coef_atb_li.b <- unlist(lapply(beta.cisList.b, function(x) x[grep("ddd_total", rownames(x)),2]))
  coef_atb_ui.b <- unlist(lapply(beta.cisList.b, function(x) x[grep("ddd_total", rownames(x)),3]))
  
  #Connectivity
  coef_s.b    <- unlist(lapply(beta.cisList.b, function(x) x[grep("S_connectivity", rownames(x)),1]))
  coef_s_li.b <-  unlist(lapply(beta.cisList.b, function(x) x[grep("S_connectivity", rownames(x)),2]))
  coef_s_ui.b <-  unlist(lapply(beta.cisList.b, function(x) x[grep("S_connectivity", rownames(x)),3]))
  
  #Number of beds
  coef_n.b    <- unlist(lapply(beta.cisList.b, function(x) x[grep("n_beds", rownames(x)),1]))
  coef_n_li.b <- unlist(lapply(beta.cisList.b, function(x) x[grep("n_beds", rownames(x)),2]))
  coef_n_ui.b <- unlist(lapply(beta.cisList.b, function(x) x[grep("n_beds", rownames(x)),3]))
  
  #Patient fragility (ward type) Category 1
  coef_ps1.b    <- unlist(lapply(beta.cisList.b, function(x) x[grep("PatStat1", rownames(x)),1]))
  coef_ps_li1.b <-  unlist(lapply(beta.cisList.b, function(x) x[grep("PatStat1", rownames(x)),2]))
  coef_ps_ui1.b <-  unlist(lapply(beta.cisList.b, function(x) x[grep("PatStat1", rownames(x)),3]))
  
  #Patient fragility (ward type) Category 2
  coef_ps2.b    <- unlist(lapply(beta.cisList.b, function(x) x[grep("PatStat2", rownames(x)),1]))
  coef_ps_li2.b <-  unlist(lapply(beta.cisList.b, function(x) x[grep("PatStat2", rownames(x)),2]))
  coef_ps_ui2.b <-  unlist(lapply(beta.cisList.b, function(x) x[grep("PatStat2", rownames(x)),3]))
  
}


###########################################################################################
# Prepare the features of the figure

# Set variant labels 
buglabels <- c("EC", "3GCREC", "CREC", "KP", "3GCRKP", "CRKP", "EB", "3GCREB",
               "CREB", "PA", "CRPA", "SA", "MRSA")

# Labels for Abau and Efae
buglabels.b <-c("AB", "CRAB", "EF", "VREF")

# Set colors for the beta coefficient points for each taxon
bugcolors <- c("violetred1", "violetred3", "violetred4",
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
# Generate Figure 1

# Panel 1 (Figure 1a) - all taxa except Abau and Efae
svg(file = "glm_global_pane1.svg", 3.5, 6)
{
  showpanes <- function(yl) {
    rect(0.75,  yl[1], 3.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(3.75,  yl[1], 6.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(6.75,  yl[1], 9.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(9.75,  yl[1], 11.25, yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(11.75, yl[1], 13.25, yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
  }
  
  par(mfrow = c(6,1))
  par(mar = c(1,4,1,4))
  
  p <- length(coef_atb)
  
  yl <- c(-100,200)
  plot(coef_ps2[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "Intensive care", bty = "n", type = "n", yaxt = "n")
  axis(2, at = c(-100, 0, 100, 200), labels = c(-100, 0, 100, 200))
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_ps_li2[ord], 1:p, coef_ps_ui2[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_ps2[ord], pch = 19, col = bugcolors, cex = 1.5)
  
  yl <- c(-200,600)
  plot(coef_ps1[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "Progressive care", bty = "n", type = "n", yaxt = "n")
  axis(2, at = c(-200, 200, 600))
  showpanes(yl)
  
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_ps_li1[ord], 1:p, coef_ps_ui1[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_ps1[ord], pch = 19, col = bugcolors, cex = 1.5)
  
  yl <- c(-60,50)
  plot(coef_n[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "Ward size", bty = "n", type = "n", yaxt = "n")
  axis(2, at = c(-40, 0, 40))
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_n_li[ord], 1:p, coef_n_ui[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_n[ord], pch = 19, col = bugcolors, cex = 1.5)
  
  yl <- c(-20,40)
  plot(coef_s[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "Connectivity", bty = "n", type = "n", yaxt = "n")
  axis(2, at = c(-20, 0, 20, 40))
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_s_li[ord], 1:p, coef_s_ui[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_s[ord], pch = 19, col = bugcolors, cex = 1.5)
  
  yl <- c(-20,100)
  plot(coef_atb[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "Antibiotic use", bty = "n", type = "n", yaxt = "n")
  axis(2, at = c(0, 50, 100))
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_atb_li[ord], 1:p, coef_atb_ui[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_atb[ord], pch = 19, col = bugcolors, cex = 1.5)
  
  axis(1, at = 1:p, labels = buglabels[ord], las = 2)  
}
dev.off()


# Panel 2 (Figure 1b) - Abau and Efae
svg(file = "glm_global_pane2.svg", 1.66, 6)
{
  par(mfrow = c(6,1))
  par(mar = c(1,4,1,4))
  
  showpanes <- function(yl) {
    rect(0.75,  yl[1], 2.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(2.75,  yl[1], 4.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
  }
  
  p <- length(coef_ps1.b)
  xl <- c(0.75, 4.25)
  
  yl <- c(-100, 120)
  plot(coef_ps2.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "", bty = "n", type = "n", yaxt = "n")
  axis(2, at = c(-100, 0, 100))
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_ps_li2.b[ord], 1:p, coef_ps_ui2.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_ps2.b[ord], pch = 19, col = bugcolors.b, cex = 1.5)
  
  
  yl <- c(-100,200)
  plot(coef_ps1.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "", bty = "n", type = "n", yaxt = "n")
  axis(2, at = c(-100, 0, 100, 200))
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_ps_li1.b[ord], 1:p, coef_ps_ui1.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_ps1.b[ord], pch = 19, col = bugcolors.b, cex = 1.5)
  
  yl <- c(-80,150)
  plot(coef_n.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_n_li.b[ord], 1:p, coef_n_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_n.b[ord], pch = 19, col = bugcolors.b, cex = 1.5)
  
  yl <- c(-30, 90)
  plot(coef_s.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "", bty = "n", type = "n", yaxt = "n")
  axis(2, at = c(-20, 0, 20, 40, 60, 80))
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_s_li.b[ord], 1:p, coef_s_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_s.b[ord], pch = 19, col = bugcolors.b, cex = 1.5)
  
  yl <- c(-20,150)
  plot(coef_atb.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "", bty = "n", type = "n", yaxt = "n")
  axis(2, at = c(0, 50, 100, 150))
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_atb_li.b[ord], 1:p, coef_atb_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_atb.b[ord], pch = 19, col = bugcolors.b, cex = 1.5)
  
  axis(1, at = 1:p, labels = buglabels.b[ord], las = 2)  
}
dev.off()

# Find resulting figures in files glm_global_pane1.svg and glm_global_pane2.svg


