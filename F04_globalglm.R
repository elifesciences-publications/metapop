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

load(file = "modeldata.Rdata") # Output from F03_modeldataprep.R

#Begin with "transf.dat" from F03 script.
head(transf.dat)

#First identify each unit for the for-loop
taxonName <- unique(transf.dat$BacType)

# Store per-variant model results in a list
ResultList<-list()

# Loop over variants and populate ResultList
for(i in 1:length(taxonName)){
  selTaxon<-which(transf.dat$BacType == taxonName[i]) #Divide the data,
  # this variable will be called to identify each subset of the data
  
  # Run Poisson GLM
  metab = glm(N_patients ~ C_control + S_connectivity + n_beds + ddd_total + PatStat, data = transf.dat,subset=selTaxon, family="poisson")
  
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
  ResultList[[length(ResultList)+1]]<-resList  
}

# Use the [[]] operator to access individual results
ResultList[[1]]$name
ResultList[[1]]$betas.cis

##################################################################################################
# Figure 2 of beta coefficients and confidence intervals for each variable in each model

# Extract the betas + CI's table for each element of the list
beta.cisList <- lapply(ResultList, function (x){unlist(x[['betas.cis']])})
# Set names for each list
names.betacisList <- unlist(lapply(ResultList, function(x) x[c('name')]))
# Apply those names to the list
names(beta.cisList) <- names.betacisList

# Reorder the elements of the list, needed for consistent variant ordering
beta.cisList.order <- beta.cisList[c(1,2,3,4,5,6,7,8,9,10,11,16,17,12,13,14,15)]

# Make one list of all bacteria except Efaec and Abau, which have wider confidence intervals
beta.cisList.a <- beta.cisList.order[c(1:13)]

# Separate list for Abau and Efaec due to much larger confidence intervals
beta.cisList.b <- beta.cisList.order[c(14:17)]


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
  
  #Patient fragility (ward type)
  coef_ps    <- unlist(lapply(beta.cisList.a, function(x) x[grep("PatStat", rownames(x)),1]))
  coef_ps_li <-  unlist(lapply(beta.cisList.a, function(x) x[grep("PatStat", rownames(x)),2]))
  coef_ps_ui <-  unlist(lapply(beta.cisList.a, function(x) x[grep("PatStat", rownames(x)),3]))
  
}

# Repeat for Abau and Efae
{
  coef_atb.b    <- unlist(lapply(beta.cisList.b, function(x) x[grep("ddd_total", rownames(x)),1]))
  coef_atb_li.b <- unlist(lapply(beta.cisList.b, function(x) x[grep("ddd_total", rownames(x)),2]))
  coef_atb_ui.b <- unlist(lapply(beta.cisList.b, function(x) x[grep("ddd_total", rownames(x)),3]))
  
  coef_s.b    <- unlist(lapply(beta.cisList.b, function(x) x[grep("S_connectivity", rownames(x)),1]))
  coef_s_li.b <-  unlist(lapply(beta.cisList.b, function(x) x[grep("S_connectivity", rownames(x)),2]))
  coef_s_ui.b <-  unlist(lapply(beta.cisList.b, function(x) x[grep("S_connectivity", rownames(x)),3]))
  
  coef_n.b    <- unlist(lapply(beta.cisList.b, function(x) x[grep("n_beds", rownames(x)),1]))
  coef_n_li.b <- unlist(lapply(beta.cisList.b, function(x) x[grep("n_beds", rownames(x)),2]))
  coef_n_ui.b <- unlist(lapply(beta.cisList.b, function(x) x[grep("n_beds", rownames(x)),3]))
  
  coef_ps.b    <- unlist(lapply(beta.cisList.b, function(x) x[grep("PatStat", rownames(x)),1]))
  coef_ps_li.b <-  unlist(lapply(beta.cisList.b, function(x) x[grep("PatStat", rownames(x)),2]))
  coef_ps_ui.b <-  unlist(lapply(beta.cisList.b, function(x) x[grep("PatStat", rownames(x)),3]))
  
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
# Generate figure

# Panel 1 - all taxa except Abau and Efae
svg(file = "glm_global_pane1.svg", 4, 6)
{
  showpanes <- function(yl) {
    rect(0.75,  yl[1], 3.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(3.75,  yl[1], 6.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(6.75,  yl[1], 9.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(9.75,  yl[1], 11.25, yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(11.75, yl[1], 13.25, yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
  }
  
  par(mfrow = c(5,1))
  par(mar = c(1,4,1,4))
  
  p <- length(coef_ps)
  
  yl <- c(-0.75,0.75)
  plot(coef_ps[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "Ward type", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_ps_li[ord], 1:p, coef_ps_ui[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_ps[ord], pch = 19, col = bugcolors, cex = 1.5)
  
  yl <- c(-0.5,0.5)
  plot(coef_n[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "Ward size", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_n_li[ord], 1:p, coef_n_ui[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_n[ord], pch = 19, col = bugcolors, cex = 1.5)
  
  yl <- c(-0.2,0.2)
  plot(coef_s[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "Connectivity", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_s_li[ord], 1:p, coef_s_ui[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_s[ord], pch = 19, col = bugcolors, cex = 1.5)
  
  yl <- c(-0.45,0.45)
  plot(coef_atb[ord], ylim = yl, xaxt = "n", xlab = "", ylab = "Antibiotic use", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_atb_li[ord], 1:p, coef_atb_ui[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_atb[ord], pch = 19, col = bugcolors, cex = 1.5)
  
  axis(1, at = 1:p, labels = buglabels[ord], las = 2)  
}
dev.off()


# Panel 2 - Abau and Efae
svg(file = "glm_global_pane2.svg", 1.9, 6)
{
  par(mfrow = c(5,1))
  par(mar = c(1,4,1,4))
  
  showpanes <- function(yl) {
    rect(0.75,  yl[1], 2.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
    rect(2.75,  yl[1], 4.25,  yl[2], col = rgb(0.9,0.9,0.9,0.3), border = NA)
  }
  
  p <- length(coef_ps.b)
  xl <- c(0.75, 4.25)
  
  yl <- c(-1.5,1.5)
  plot(coef_ps.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_ps_li.b[ord], 1:p, coef_ps_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_ps.b[ord], pch = 19, col = bugcolors.b, cex = 1.5)
  
  yl <- c(-1.5,1.5)
  plot(coef_n.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_n_li.b[ord], 1:p, coef_n_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_n.b[ord], pch = 19, col = bugcolors.b, cex = 1.5)
  
  yl <- c(-0.75, 0.75)
  plot(coef_s.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_s_li.b[ord], 1:p, coef_s_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_s.b[ord], pch = 19, col = bugcolors.b, cex = 1.5)
  
  yl <- c(-1,1)
  plot(coef_atb.b[ord], xlim = xl, ylim = yl, xaxt = "n", xlab = "", ylab = "", bty = "n", type = "n")
  showpanes(yl)
  abline(0,0, lty = 2, col = "lightgrey")
  arrows(1:p, coef_atb_li.b[ord], 1:p, coef_atb_ui.b[ord], length = errorbar_width, angle = 90, code = 3, col = "darkgrey")
  points(coef_atb.b[ord], pch = 19, col = bugcolors.b, cex = 1.5)
  
  axis(1, at = 1:p, labels = buglabels.b[ord], las = 2)  
}
dev.off()

# Find resulting figures in files glm_global_pane1.svg and glm_global_pane2.svg


