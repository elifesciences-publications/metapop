#######Correlation figures ########################
### Need to make correlation figure with Control on X axis vs. on Y axis:
### 1. Occurence of each variant
### 2. Use of each antibiotic
### 3. Use of total antibiotics
### 4. Connectivity
### -> Include coefficient of correlation

## Labels for each plot
buglabels <- list(ESCCOL_S = "EC",
                  ESCCOL_C3G_R = "3GCREC",
                  ESCCOL_CARBA_R = "CREC",
                  KLEPNE_S = "KP",KLEPNE_C3G_R="3GCRKP",
                  KLEPNE_CARBA_R = "CRKP",
                  ENTCLO_S = "EB",ENTCLO_C3G_R = "3GCREB",
                  ENTCLO_CARBA_R="CREB",PSEAER_S="PA",
                  PSEAER_CARBA_R="CRPA",
                  ACIBAU_CARBA_S="AB",ACIBAU_CARBA_R="CRAB",
                  ENCFAC_VANCO_S="EF",ENCFAC_VANCO_R="VREF",
  STAAUR_OXA_S = "SA",STAAUR_OXA_R="MRSA"
)

########################################
### Begin with raw data ######

head(mod.dat.raw) # From F01_dataprep

raw.dat <- mod.dat.raw
head(raw.dat)

# Total ATB use (is currently in ddd per bed)
raw.dat$TotATB <- raw.dat$ddd_total*raw.dat$n_beds
head(raw.dat)

# Create totals for N_patients and Connectivity per ward
raw.dat2 <- raw.dat %>% 
  group_by(ward)%>%
  summarise(Tot_NPats = sum(N_patients), 
            TotConn = sum(S_connectivity),TotControl = sum(C_control))
raw.dat3 <- right_join(raw.dat, raw.dat2, by="ward")
head(raw.dat3)


##################################################################
###### Figures for Log-transformed data, with axes "libre" #######

## Fig S2 - Control Vs. Incidence 
##Labels
xlabplot <- "log2 Incidence control"
ylabplot <- "log2 No. episodes"

{
  svg("IncidenceVsControl.CI.svg",4.8,12) #4X10 was square
  par(mfrow=c(6,3))
  for (i in unique(raw.dat3$BacType)){
    pointcol <- rgb(0,0,0.8,0.25)
    d_subset <- raw.dat3[raw.dat3$BacType==i,] #keep only rows for that file number
    xx <- log2(d_subset$C_control)
    yy <- log2(d_subset$N_patients)
    yy[!is.finite(yy)] <- min(yy[is.finite(yy)]) - 1
    xx[!is.finite(xx)] <- min(xx[is.finite(xx)]) - 1
    #xylim <- c(-3, max(c(xx, yy)))
    
    plot(xx, yy, ylab = ylabplot, 
         xlab = xlabplot, main = sprintf("%s", buglabels[[i]])
         
         #, xlim = xylim, ylim = xylim
         , pch = 19, col = pointcol
    )
    #abline(lm(yy ~ xx), col = "blue") #### NOTE: Uncomment for raw data figs
    # cortext <- paste('Cor =',round(cor(d_subset$C_control, d_subset$N_patients),2))
    #p.value.cor <- cor.test(xx, yy)$p.value
    #cortext <- sprintf("R² = %.2f, p = %.1e", cor(xx, yy)^2, p.value.cor)
    corr.ci <- cor.test(xx, yy)$conf.int
    cortext <- sprintf("R² = %.2f (%.2f, %.2f)", 
                       cor(xx, yy)^2, corr.ci[1]^2, corr.ci[2]^2)
    legend("top",inset = c(-0.2,-0.2),xpd=TRUE, legend = "", 
           title=cortext,cex=0.7,bty="n",
           box.lwd=0.5)
    
  }
  xx2 <- log2(d_subset$TotControl)
  yy2 <- log2(d_subset$Tot_NPats)
  yy2[!is.finite(yy2)] <- min(yy2[is.finite(yy2)]) - 1
  xx2[!is.finite(xx2)] <- min(xx2[is.finite(xx2)]) - 1
  #xylim2 <- c(-3, max(c(xx2, yy2)))
  plot(xx2, yy2,ylab = ylabplot, 
       xlab = xlabplot, main = "All Bacteria",
       #xlim = xylim2, ylim = xylim2, 
       pch = 19, col = pointcol)
  #p.value.cor <- cor.test(xx2, yy2)$p.value
  #cortext <- sprintf("R² = %.2f, p = %.1e", cor(xx2, yy2)^2, p.value.cor)
  corr.ci2 <- cor.test(xx2, yy2)$conf.int
  cortext2 <- sprintf("R² = %.2f (%.2f, %.2f)", 
                     cor(xx2, yy2)^2, corr.ci2[1]^2, corr.ci2[2]^2)
  legend("top",inset = c(-0.2,-0.2),xpd=TRUE
         , legend = "", title=cortext2,cex=0.7,bty="n",
         box.lwd=0.5)
  dev.off()
}##

## Fig S3 - Control Vs. Antibiotic use 
##Labels
xlabplot <- "log2 Incidence control"
ylabplot <- "log2 Antibiotic use"

{
  svg("ControlVsATB.CI.svg",4.8,12) #4X10 was square
  par(mfrow=c(6,3))
  for (i in unique(raw.dat3$BacType)){
    pointcol <- rgb(0,0,0.8,0.25)
    d_subset <- raw.dat3[raw.dat3$BacType==i,] #keep only rows for that file number
    xx <- log2(d_subset$C_control)
    yy <- log2(d_subset$ddd_total)
    yy[!is.finite(yy)] <- min(yy[is.finite(yy)]) - 1
    xx[!is.finite(xx)] <- min(xx[is.finite(xx)]) - 1
    #xylim <- c(-3, max(c(xx, yy)))
    
    plot(xx, yy, ylab = ylabplot, 
         xlab = xlabplot, main = sprintf("%s", buglabels[[i]])
         
         #, xlim = xylim, ylim = xylim
         , pch = 19, col = pointcol
    )
    #abline(lm(yy ~ xx), col = "blue") #### NOTE: Uncomment for raw data figs
    # cortext <- paste('Cor =',round(cor(d_subset$C_control, d_subset$N_patients),2))
    #p.value.cor <- cor.test(xx, yy)$p.value
    #cortext <- sprintf("R² = %.2f, p = %.1e", cor(xx, yy)^2, p.value.cor)
    corr.ci <- cor.test(xx, yy)$conf.int
    cortext <- sprintf("R² = %.2f (%.2f, %.2f)", 
                       cor(xx, yy)^2, corr.ci[1]^2, corr.ci[2]^2)
    legend("top",inset = c(-0.2,-0.2),xpd=TRUE, legend = "", 
           title=cortext,cex=0.7,bty="n",
           box.lwd=0.5)
    
  }
  xx2 <- log2(d_subset$TotControl)
  yy2 <- log2(d_subset$ddd_total)
  yy2[!is.finite(yy2)] <- min(yy2[is.finite(yy2)]) - 1
  xx2[!is.finite(xx2)] <- min(xx2[is.finite(xx2)]) - 1
  #xylim2 <- c(-3, max(c(xx2, yy2)))
  plot(xx2, yy2,ylab = ylabplot, 
       xlab = xlabplot, main = "All Bacteria",
       #xlim = xylim2, ylim = xylim2, 
       pch = 19, col = pointcol)
  #p.value.cor <- cor.test(xx2, yy2)$p.value
  #cortext <- sprintf("R² = %.2f, p = %.1e", cor(xx2, yy2)^2, p.value.cor)
  corr.ci2 <- cor.test(xx2, yy2)$conf.int
  cortext2 <- sprintf("R² = %.2f (%.2f, %.2f)", 
                     cor(xx2, yy2)^2, corr.ci[1]^2, corr.ci[2]^2)
  legend("top",inset = c(-0.2,-0.2),xpd=TRUE
         , legend = "", title=cortext2,cex=0.7,bty="n",
         box.lwd=0.5)
  dev.off()
}##

## Fig S4 - Antibiotic use Vs. Incidence 
##Labels
xlabplot <- "log2 Antibiotic use"
ylabplot <- "log2 No. episodes"

{
  svg("ATBvsIncidence.CI.svg",4.8,12) #4X10 was square
  par(mfrow=c(6,3))
  for (i in unique(raw.dat3$BacType)){
    pointcol <- rgb(0,0,0.8,0.25)
    d_subset <- raw.dat3[raw.dat3$BacType==i,] #keep only rows for that file number
    xx <- log2(d_subset$ddd_total)
    yy <- log2(d_subset$N_patients)
    yy[!is.finite(yy)] <- min(yy[is.finite(yy)]) - 1
    xx[!is.finite(xx)] <- min(xx[is.finite(xx)]) - 1
    #xylim <- c(-3, max(c(xx, yy)))
    
    plot(xx, yy, ylab = ylabplot, 
         xlab = xlabplot, main = sprintf("%s", buglabels[[i]])
         
         #, xlim = xylim, ylim = xylim
         , pch = 19, col = pointcol
    )
    #abline(lm(yy ~ xx), col = "blue") #### NOTE: Uncomment for raw data figs
    # cortext <- paste('Cor =',round(cor(d_subset$C_control, d_subset$N_patients),2))
    #p.value.cor <- cor.test(xx, yy)$p.value
    #cortext <- sprintf("R² = %.2f, p = %.1e", cor(xx, yy)^2, p.value.cor)
    corr.ci <- cor.test(xx, yy)$conf.int
    cortext <- sprintf("R² = %.2f (%.2f, %.2f)", 
                       cor(xx, yy)^2, corr.ci[1]^2, corr.ci[2]^2)
    legend("top",inset = c(-0.2,-0.2),xpd=TRUE, legend = "", 
           title=cortext,cex=0.7,bty="n",
           box.lwd=0.5)
    
  }
  xx2 <- log2(d_subset$ddd_total)
  yy2 <- log2(d_subset$Tot_NPats)
  yy2[!is.finite(yy2)] <- min(yy2[is.finite(yy2)]) - 1
  xx2[!is.finite(xx2)] <- min(xx2[is.finite(xx2)]) - 1
  #xylim2 <- c(-3, max(c(xx2, yy2)))
  plot(xx2, yy2,ylab = ylabplot, 
       xlab = xlabplot, main = "All Bacteria",
       #xlim = xylim2, ylim = xylim2, 
       pch = 19, col = pointcol)
  #p.value.cor <- cor.test(xx2, yy2)$p.value
  #cortext <- sprintf("R² = %.2f, p = %.1e", cor(xx2, yy2)^2, p.value.cor)
  corr.ci2 <- cor.test(xx2, yy2)$conf.int
  cortext2 <- sprintf("R² = %.2f (%.2f, %.2f)", 
                     cor(xx2, yy2)^2, corr.ci2[1]^2, corr.ci2[2]^2)
  legend("top",inset = c(-0.2,-0.2),xpd=TRUE
         , legend = "", title=cortext2,cex=0.7,bty="n",
         box.lwd=0.5)
  dev.off()
}##

## Fig S5 - Connectivity Vs. Incidence 
##Labels
xlabplot <- "log2 Connectivity"
ylabplot <- "log2 No. episodes"

{
  svg("ConnectivityVsIncidence.CI.svg",4.8,12) #4X10 was square
  par(mfrow=c(6,3))
  for (i in unique(raw.dat3$BacType)){
    pointcol <- rgb(0,0,0.8,0.25)
    d_subset <- raw.dat3[raw.dat3$BacType==i,] #keep only rows for that file number
    xx <- log2(d_subset$S_connectivity)
    yy <- log2(d_subset$N_patients)
    yy[!is.finite(yy)] <- min(yy[is.finite(yy)]) - 1
    xx[!is.finite(xx)] <- min(xx[is.finite(xx)]) - 1
    #xylim <- c(-3, max(c(xx, yy)))
    
    plot(xx, yy, ylab = ylabplot, 
         xlab = xlabplot, main = sprintf("%s", buglabels[[i]])
         
         #, xlim = xylim, ylim = xylim
         , pch = 19, col = pointcol
    )
    #abline(lm(yy ~ xx), col = "blue") #### NOTE: Uncomment for raw data figs
    # cortext <- paste('Cor =',round(cor(d_subset$C_control, d_subset$N_patients),2))
    #p.value.cor <- cor.test(xx, yy)$p.value
    #cortext <- sprintf("R² = %.2f, p = %.1e", cor(xx, yy)^2, p.value.cor)
    corr.ci <- cor.test(xx, yy)$conf.int
    cortext <- sprintf("R² = %.2f (%.2f, %.2f)", 
                       cor(xx, yy)^2, corr.ci[1]^2, corr.ci[2]^2)
    legend("top",inset = c(-0.2,-0.2),xpd=TRUE, legend = "", 
           title=cortext,cex=0.7,bty="n",
           box.lwd=0.5)
    
  }
  xx2 <- log2(d_subset$TotConn)
  yy2 <- log2(d_subset$Tot_NPats)
  yy2[!is.finite(yy2)] <- min(yy2[is.finite(yy2)]) - 1
  xx2[!is.finite(xx2)] <- min(xx2[is.finite(xx2)]) - 1
  #xylim2 <- c(-3, max(c(xx2, yy2)))
  plot(xx2, yy2,ylab = ylabplot, 
       xlab = xlabplot, main = "All Bacteria",
       #xlim = xylim2, ylim = xylim2, 
       pch = 19, col = pointcol)
  #p.value.cor <- cor.test(xx2, yy2)$p.value
  corr.ci2 <- cor.test(xx2, yy2)$conf.int
  cortext2 <- sprintf("R² = %.2f (%.2f, %.2f)", 
                     cor(xx2, yy2)^2, corr.ci2[1]^2, corr.ci2[2]^2)
  legend("top",inset = c(-0.2,-0.2),xpd=TRUE
         , legend = "", title=cortext2,cex=0.7,bty="n",
         box.lwd=0.5)
  dev.off()
}##
