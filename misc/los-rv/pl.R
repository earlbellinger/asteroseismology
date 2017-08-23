library(magicaxis)
library(deming)

source(file.path('..', '..', 'scripts', 'utils.R'))

speed_of_light = 299792 # km/s



shift <- function(P_o, kms, dkms) {
    nu_o <- 1 / P_o
    v_r.over.c_upper = (kms+dkms)/speed_of_light
    v_r.over.c_lower = (kms-dkms)/speed_of_light
    doppler_beta_upper = sqrt((1+v_r.over.c_upper)/(1-v_r.over.c_upper))
    doppler_beta_lower = sqrt((1+v_r.over.c_lower)/(1-v_r.over.c_lower))
    nu_s_upper <- doppler_beta_upper * nu_o
    nu_s_lower <- doppler_beta_lower * nu_o
    P_s_upper <- 1 / nu_s_upper
    P_s_lower <- 1 / nu_s_lower
    P_s <- (P_s_upper + P_s_lower)/2
    sd <- P_s_lower - P_s_upper
    dP <- P_o - P_s
    return(list(sd=sd, dP=dP))
}

radial_velocities <- data.frame(
 galaxy=c("LMC", "SMC", "M101", "N1015", "N1309", "N1365", "N1448", "N2442",
          "N3021", "N3370", "N3447", "N3972", "N3982", "N4038", "N4258", 
          "N4424", "N4536", "N4639", "N5584", "N5917", "N7250"),
    kms=c(262.2, 145.6, 243, 2132.4, 2625.7, 1637.1, 1165.5, 1449.3, 1600, 
          1279, 1064, 838, 1110, 1672.0, 462, 438, 1808, 978, 1652, 1934.1, 
          1168),
   dkms=c(3.4, 0.6, 5, 0.1, 0.1, 0.1, 0.1, 0.1, 3, 3, 0.1, 10, 3, 1, 27, 7, 1, 
          16, 3, 0.1, 1))

all.ceps <- read.fwf('apjaa3610t5_mrt.txt', 
    widths=c(5, 9, 9, 8, 6, 6, 3, 6, 3, 5, 5, 5, 5, 2, 3)+1, skip=40,
    col.names=c("Galaxy", "RAdeg", "DEdeg", "ID", "Period",
                "F555Wmag", "e_F555Wmag", "F814Wmag", "e_F814Wmag",
                "AmpV", "AmpI", "AmpWh", "Z", "Flag", "Ref"),
    stringsAsFactors=F)

mean.period <- do.call(rbind, Map(function(galaxy) {
    data.frame(galaxy=trimws(galaxy), 
               period=mean(all.ceps[all.ceps$Galaxy == galaxy,]$Period))
}, galaxy=unique(all.ceps$Galaxy)))

mean.period <- rbind(mean.period,
    data.frame(galaxy='LMC', 
       period=mean(read.table('cepF-lmc.dat', header=T, stringsAsFactors=F)$P)),
    data.frame(galaxy='SMC', 
       period=mean(read.table('cepF-smc.dat', header=T, stringsAsFactors=F)$P)))

radial_velocities <- merge(radial_velocities, mean.period, by=c('galaxy'))

for (ii in 1:nrow(radial_velocities)) {
    dPsd <- with(radial_velocities[ii,], shift(period, kms, dkms))
    cat(paste(
        round(radial_velocities[ii,]$period,5),
        '&',
        round(dPsd$dP*24,5), 
        "$\\pm$", 
        round(dPsd$sd*24,5), "\\\\\n"))
}



multiphase <- read.table('cepF-lmc-lasso.dat', header=T, stringsAsFactors=F)
phase.cols <- multiphase[,grepl('Phase', names(multiphase))]
LMC.ml <- data.frame(ID=multiphase$ID,
                     I.mean = apply(phase.cols, 1, mean),
                     I.sd =   sqrt(apply(phase.cols, 1, var)))

multiphase <- read.table('cepF-smc-lasso.dat', header=T, stringsAsFactors=F)
phase.cols <- multiphase[,grepl('Phase', names(multiphase))]
SMC.ml <- data.frame(ID=multiphase$ID,
                     I.mean = apply(phase.cols, 1, mean),
                     I.sd =   sqrt(apply(phase.cols, 1, var)))

LMC <- merge(LMC.ml, merge(
    read.table('cepF-lmc.dat', header=T, stringsAsFactors=F),
    read.table('red-lmc.dat', header=T, stringsAsFactors=F),
    by='ID'), by='ID')
SMC <- merge(SMC.ml, merge(
    read.table('cepF-smc.dat', header=T, stringsAsFactors=F),
    read.table('red-smc.dat', header=T, stringsAsFactors=F),
    by='ID'), by='ID')

LMC <- LMC[complete.cases(LMC),]
SMC <- SMC[complete.cases(SMC),]

# Uncorrected LMC and SMC
PL <- function(LMC.I, LMC.dI, LMC.P, LMC.dP,
               SMC.I, SMC.dI, SMC.P, SMC.dP, ...,
               text.cex=1, mgp=utils.mgp, font="Times New Roman") {
    
    #par(mar=c(4, 5, 1, 1))
    par(mgp=mgp)
    
    plot(NA, axes=F, tcl=0, 
         xaxs='i', yaxs='i', 
         xlim=c(0.0001, log10(max(c(SMC$P, LMC$P)))),
         ylim=rev(range(SMC.I, LMC.I)),
         xlab=expression("Period"~P/"days"),
         ylab=expression("I-band Magnitude"))
    #magaxis(1:4, labels=c(1,1,0,0), tcl=-0.25, las=1, family=font, mgp=mgp,
    #        unlog='x', cex.axis=text.cex)
    
    SMC.logP <- log10(SMC.P)
    SMC.dlogP <- SMC.dP / SMC.P
    
    LMC.logP <- log10(LMC.P)
    LMC.dlogP <- LMC.dP / LMC.P
    
    smc <- data.frame(I=SMC.I, dI=SMC.dI, logP=SMC.logP, dlogP=SMC.dlogP)
    lmc <- data.frame(I=LMC.I, dI=LMC.dI, logP=LMC.logP, dlogP=LMC.dlogP)
    
    repeat {
        SMC.fmla <- smc$I ~ smc$logP
        residuals <- resid(lm(SMC.fmla)) / smc$dI
        outliers <- boxplot.stats(residuals, coef=3)$out
        if (length(outliers) == 0) break
        smc <- smc[-which(residuals %in% outliers),]
    }
    
    repeat {
        LMC.fmla <- lmc$I ~ lmc$logP
        residuals <- resid(lm(LMC.fmla)) / lmc$dI
        outliers <- boxplot.stats(residuals, coef=3)$out
        if (length(outliers) == 0) break
        lmc <- lmc[-which(residuals %in% outliers),]
    }
    
    #SMC.fmla <- SMC.I ~ SMC.logP
    #LMC.fmla <- LMC.I ~ LMC.logP
    #arrows(SMC.logP, SMC.I - SMC.dI,
    #       SMC.logP, SMC.I + SMC.dI,
    #       length=0, lwd=1, angle=90, code=3, col="darkred")
    points(SMC.fmla, pch=20, cex=0.75, col='#8b000066') 
    #arrows(LMC.logP, LMC.I - LMC.dI,
    #       LMC.logP, LMC.I + LMC.dI,
    #       length=0, lwd=1, angle=90, code=3)
    points(LMC.fmla, pch=20, cex=0.75, col="#00000066") 
    
    SMC.lm <- lm(SMC.fmla)
    LMC.lm <- lm(LMC.fmla)
    
    #print(summary(SMC.lm))
    #print(summary(LMC.lm))
    
    LMC.coef <- coef(summary(LMC.lm))
    SMC.coef <- coef(summary(SMC.lm))
    Z <- (LMC.coef[2,1] - SMC.coef[2,1]) / 
         sqrt( LMC.coef[2,2]**2 + SMC.coef[2,2]**2 )
    print(paste("Z =", Z))
    print(paste("p =", dnorm(Z)))
    
    abline(SMC.lm, col='white', lwd=3)
    abline(SMC.lm, col='darkred', lwd=2)
    abline(LMC.lm, col='white', lwd=3)
    abline(LMC.lm, lwd=2)
    
    for (galaxy in c("LMC", "SMC")) {
        for (side in c("All", "Short", "Long")) {
            dta <- if (galaxy == "LMC") lmc else smc
            sel <- if (side == "All") T else 
                   if (side == "Short") dta$logP < 1 else
                                        dta$logP >= 1 
            dta <- dta[sel,]

            mdl <- deming(I ~ logP, data=dta, xstd=dlogP, ystd=dI) 
            
            cat(paste0(galaxy, ' & ', side, ' & ', 
                       mdl$coef[[2]], ' $\\pm$ ', sqrt(mdl$var[[4]]), ' & ', 
                       mdl$coef[[1]], ' $\\pm$ ', sqrt(mdl$var[[1]]), ' & \n'))
            
            #ab <- summary(lm(I ~ logP, data=dta))$coefficients
            
            #cat(paste0(galaxy, ' & ', side, ' & ', 
            #           ab[[2,1]], ' $\\pm$ ', ab[[2,2]], ' & ', 
            #           ab[[1,1]], ' $\\pm$ ', ab[[1,2]], ' & \n'))
        }
    }
    
    #SMC.mdl <- deming(SMC.fmla, xstd=SMC.dlogP, ystd=SMC.dI, conf=.99) 
    #LMC.mdl <- deming(LMC.fmla, xstd=LMC.dlogP, ystd=LMC.dI, conf=.99) 
    #abline(SMC.mdl, col='darkred')#, untf=T)
    #abline(LMC.mdl)#, untf=T)
    legend("topleft", col=c("darkred", "black"), pch=20, 
           legend=c("SMC", "LMC"), bty='n', cex=0.8*text.cex)
    magaxis(side=2:4, family=font, tcl=-0.25, labels=c(1,0,0), 
        mgp=mgp+c(0, 0.2, 0), las=1, cex.axis=text.cex, unlog='x')
    magaxis(side=1, family=font, tcl=-0.25, labels=T, mgp=mgp,#-c(0, 0.2, 0), 
        las=1, cex.axis=text.cex, unlog='x')
}

make_plots(PL, 'PL', 
   LMC.I=LMC$I.mean, LMC.dI=LMC$I.sd, LMC.P=LMC$P, LMC.dP=LMC$dP,
   SMC.I=SMC$I.mean, SMC.dI=SMC$I.sd, SMC.P=SMC$P, SMC.dP=SMC$dP)

PL(LMC$I.mean, LMC$I.sd, LMC$P, LMC$dP,
   SMC$I.mean, SMC$I.sd, SMC$P, SMC$dP)


## Correct for reddening
#LMC.red <- read.table('red-lmc.dat', header=T, stringsAsFactors=F)
#SMC.red <- read.table('red-smc.dat', header=T, stringsAsFactors=F)

#LMC.corr <- merge(LMC, LMC.red, by='ID')
#SMC.corr <- merge(SMC, SMC.red, by='ID')

#LMC.corr$I <- LMC.corr$I - 1.41 * LMC.corr$ev_i
LMC$I.red <- LMC$I.mean - 1.41 * LMC$ev_i

mean.red <- mean(SMC$ev_i[SMC$ev_i != 0])
SMC[SMC$ev_i == 0,]$ev_i <- mean.red
SMC$I.red <- SMC$I - 0.48 - 1.41 * SMC$ev_i
PL(LMC$I.red, LMC$I.sd, LMC$P, LMC$dP,
   SMC$I.red, SMC$I.sd, SMC$P, SMC$dP)

## Correct for line-of-sight radial velocity
LMC.rv <- 262.2 # km/s 
LMC.drv <- 3.4 
SMC.rv <- 145.6 
SMC.drv <- 0.6 

#LMC.rvcorr <- LMC.corr
#SMC.rvcorr <- SMC.corr

LMC.doppler_beta <- LMC.rv/speed_of_light
LMC.dbeta <- LMC.drv/speed_of_light
LMC.doppler_shift <- sqrt((1+LMC.doppler_beta) / (1-LMC.doppler_beta))
LMC.dshift <- sqrt((1+LMC.dbeta) / (1-LMC.dbeta)) - 1
LMC$P.corr <- LMC$P / LMC.doppler_shift
LMC$dP.corr <- LMC$P.corr * sqrt( (LMC$dP / LMC$P)**2 +
                                  (LMC.dshift / LMC.doppler_shift)**2 )

SMC.doppler_beta <- SMC.rv/speed_of_light
SMC.dbeta <- SMC.drv/speed_of_light
SMC.doppler_shift <- sqrt((1+SMC.doppler_beta) / (1-SMC.doppler_beta))
SMC.dshift <- sqrt((1+SMC.dbeta) / (1-SMC.dbeta)) - 1
SMC$P.corr <- SMC$P / SMC.doppler_shift
SMC$dP.corr <- SMC$P.corr * sqrt( (SMC$dP / SMC$P)**2 +
                                  (SMC.dshift / SMC.doppler_shift)**2 )



cat('Mean increase in uncertainty (factor):', mean(LMC$dP.corr/LMC$dP))
plot(LMC$P, (LMC$dP.corr)/LMC$dP)

cairo_pdf("PL-doppler.pdf", family='Times New Roman',
          width=6.97522/2, height=4.17309*2/3)
par(mar=c(2.5,3,1,1), mgp=c(1.5, 0.15, 0))
PL(LMC$I.red, LMC$I.sd, LMC$P.corr, LMC$dP.corr,
   SMC$I.red, SMC$I.sd, SMC$P.corr, SMC$dP.corr)
dev.off()
#PL(LMC.rvcorr, SMC.rvcorr)

plot(SMC$P.corr, SMC$dP.corr / SMC$dP, pch=20, cex=0.75, 
     xlab="Period", ylab="New Uncertainty / Old Uncertainty", 
     xaxs='i', yaxs='i', xlim=c(0, 15), ylim=c(0, 16))

plot(LMC$P.corr, LMC$dP.corr / LMC$dP, pch=20, cex=0.75, 
     xlab="Period", ylab="New Uncertainty / Old Uncertainty", 
     xaxs='i', yaxs='i', xlim=c(0, 15), ylim=c(0, 50))
abline(h=mean(LMC$dP.corr / LMC$dP), lty=2)


sigma <- function(...,
                  text.cex=1, mgp=utils.mgp, mar=utils.mar, 
                  font="Times New Roman") {
    par(mar=mar+c(0, 0.1, 0, 0))
    LMC.sigma <- abs(LMC$P - LMC$P.corr) / LMC$dP
    SMC.sigma <- abs(SMC$P - SMC$P.corr) / SMC$dP
    plot(NA, axes=F, pch=20, cex=0.75, 
         xlab=expression("Period"~P/"days"), 
         ylab="", 
         xaxs='i', yaxs='i', 
         xlim=c(0.00001, log10(max(c(SMC$P, LMC$P)))),
         ylim=c(1, 4))
         #ylim=c(0, 4000))
    points(log10(SMC$P), log10(SMC.sigma), pch=4, cex=0.33, col='#b36200')
    points(log10(LMC$P), log10(LMC.sigma), pch=3, cex=0.33, 
        col=adjustcolor('black', alpha=0.4))
    legend("topright", col=c('#b36200', 'black'), pch=c(4,3), 
           inset=c(0.04, 0.04), legend=c("SMC", "LMC"), cex=text.cex)
    magaxis(side=3:4, family=font, tcl=-0.25, labels=F,
        mgp=mgp+c(0, 0.2, 0), las=1, cex.axis=text.cex, unlog='xy')
    
    magaxis(side=1, family=font, tcl=-0.25, labels=T, mgp=mgp, #-c(0, 0.2, 0), 
        las=1, cex.axis=text.cex, unlog='xy')
    
    magaxis(side=2, family=font, tcl=-0.25, labels=T, mgp=mgp+c(0, 0.2, 0), 
        las=1, cex.axis=text.cex, unlog='xy')
        
    mtext(expression("Sigma distance"~"|"*P["o"]-P["s"]*"|"/sigma["o"]),
        side=2, line=2.25, cex=text.cex)
}
make_plots(sigma, 'sigma')




cairo_pdf("sigma.pdf", family='Times New Roman',
          width=6.97522/2, height=4.17309*2/3)
par(mar=c(3,3.5,1,1), mgp=c(2.1, 0.15, 0))
sigma()
dev.off()





sigma <- function(stars, col.pal, ...,
                  text.cex=1, mgp=utils.mgp, mar=utils.mar, 
                  font=utils.font) {
    par(mar=mar+c(0, 0.1, 0, 0))
    MC.sigma <- abs(stars$P - stars$P.corr) / stars$dP
    plot(NA, axes=F, pch=20, cex=0.75, 
         xlab=expression("Period"~P/"days"), 
         ylab="", 
         #xaxs='i', 
         yaxs='i', 
         xlim=c(0.000001, log10(50)),
         ylim=c(1, 4))
    points(log10(stars$P), log10(MC.sigma), pch=20, cex=0.2, col=col.pal)
    #legend("topright", col=c('#b36200', 'black'), pch=c(4,3), 
    #       inset=c(0.04, 0.04), legend=c("SMC", "LMC"), cex=text.cex)
    #magaxis(side=3:4, family=font, tcl=-0.25, labels=F,
    #    mgp=mgp+c(0, 0.2, 0), las=1, cex.axis=text.cex, unlog='xy')
    
    magaxis(side=1, family=font, tcl=-0.25, labels=T, mgp=mgp, #-c(0, 0.2, 0), 
        las=1, cex.axis=text.cex, unlog='xy')
    
    magaxis(side=2, family=font, tcl=-0.25, labels=T, mgp=mgp+c(0, 0.2, 0), 
        las=1, cex.axis=text.cex, unlog='xy')
        
    mtext(expression("Sigma distance"~"|"*P["o"]-P["s"]*"|"/sigma["o"]),
        side=2, line=2.25, cex=text.cex)
}
make_plots(sigma, 'sigma-LMC', stars=LMC, col.pal='black', wide=F, tall=F)
make_plots(sigma, 'sigma-SMC', stars=SMC, col.pal='black', wide=F, tall=F)





plot(LMC$P.corr, 100*(LMC$P - LMC$P.corr)/LMC$P, pch=20, cex=0.75, 
     xlab="Period (days)", ylab="New Period - Old Period (minutes)", 
     xaxs='i', yaxs='i', xlim=c(0, 15))



LMC.diff <- abs(LMC$P - LMC$P.corr) / LMC$P.corr
SMC.diff <- abs(SMC$P - SMC$P.corr) / SMC$P.corr
mean(LMC.diff*100)
mean(SMC.diff*100)
fivenum(c((LMC$P - LMC$P.corr)/LMC$dP, (SMC$P - SMC$P.corr)/SMC$dP))


make_plots(PL, 'PL', 
   LMC.I=LMC$I.red, LMC.dI=LMC$I.sd, LMC.P=LMC$P.corr, LMC.dP=LMC$dP.corr,
   SMC.I=SMC$I.red, SMC.dI=SMC$I.sd, SMC.P=SMC$P.corr, SMC.dP=SMC$dP.corr)

PL(LMC$I.red, LMC$I.sd, LMC$P.corr, LMC$dP.corr,
   SMC$I.red, SMC$I.sd, SMC$P.corr, SMC$dP.corr)








### NGC 5584
DF <- read.table('NGC5584.dat', header=1, stringsAsFactors=F)
shifts <- shift(DF$Period, 1652, 3)

M <- DF$F814Wmag - 1.45*(DF$F555Wmag - DF$F814Wmag) #- 5*log10(22.49 * 10**6) + 5
dM <- with(DF, 
    sqrt(      (e_F814Wmag/1000)**2 + 
          (1.45*e_F555Wmag/1000)**2 + 
          (1.45*e_F814Wmag/1000)**2 ) )
    #+ 
    #(10**6*0.07/(22.49*10**6*log(10)))**2  ) )
logP <- log10(DF$Period)

shifted <- DF$Period - shifts$dP
logP.2 <- log10(shifted)
dlogP <- shifts$sd / shifted / log(10)

#PLR <- M ~ logP

DF.1 <- data.frame(M, logP, logP.2, dM, dlogP)
repeat {
    mdl.1 <- lm(M ~ logP, data=DF.1, weights=1/dM**2)
    resids <- resid(mdl.1) / DF.1$dM
    if (!(any(resids > 3))) break
    DF.1 <- DF.1[-which.max(resids),]
}

a <- summary(mdl.1)$coefficients[[2]]
b <- summary(mdl.1)$coefficients[[1]]
da <- summary(mdl.1)$coefficients[[4]]
db <- summary(mdl.1)$coefficients[[3]]

#DF.2 <- data.frame(M, logP=logP.2, dM, dlogP)
#repeat {
mdl.2 <- deming(M ~ logP.2, xstd=dlogP, ystd=dM, data=DF.1) 
#    resids <- resid(mdl.2) / DF.2$dM
#    if (!(any(resids > 4))) break
#    DF.2 <- DF.2[-which.max(resids),]
#}

dplr.a <- mdl.2$coef[[2]]
dplr.da <- sqrt(mdl.2$var[[4]])
dplr.b <- mdl.2$coef[[1]]
dplr.db <- sqrt(mdl.2$var[[1]])

print(paste((da-dplr.da)/da*100, (db-dplr.db)/db*100))
print(pf(da**2/dplr.da**2, nrow(DF.1)-1, nrow(DF.1)-1))
print(pf(db**2/dplr.db**2, nrow(DF.1)-1, nrow(DF.1)-1))
#print(dchisq((nrow(DF.1)-1)*(da/dplr.da)**2, nrow(DF.1)-1))
#print(dchisq((nrow(DF.1)-1)*(db/dplr.db)**2, nrow(DF.1)-1))

logP <- DF.1$logP
logP.2 <- DF.1$logP.2
M <- DF.1$M
dM <- DF.1$dM
dlogP <- DF.1$dlogP

PL <- function(..., text.cex=1, mgp=utils.mgp, font="Times New Roman",
        mar=utils.mar) {
    par(mgp=mgp, mar+c(0, 0.1, 0, 0))
    plot(NA, axes=F, tcl=0, 
         xlim=c(log10(20), 2),
         ylim=rev(range(M, M+dM, M-dM)),
         xlab=expression("Period"~P/"days"),
         ylab="")
    
    arrows(logP, M-dM, logP, M+dM, 
        length=0, lwd=1.5, angle=90, code=3, col=red)
    
    arrows(logP.2, M-dM, logP.2, M+dM, 
        length=0, lwd=1.5, angle=90, code=3, col='#3E8EDE')
    arrows(logP.2-dlogP, M, logP.2+dlogP, M, 
        length=0, lwd=1.5, angle=90, code=3, col='#3E8EDE')
    
    abline(mdl.1, col=red, lty=2, lwd=3)
    abline(mdl.2, col='#3E8EDE', lty=4, lwd=2)
    
    #newx <- c(10, 200)
    #conf.int <- predict(mdl.1, newdata=data.frame(logP=newx), interval=c("confidence"))
    #lines(newx, conf.int[,2], col=red,     lty=2)
    #lines(newx, conf.int[,3], col=red,     lty=2)
    
    legend("topleft", col=c(red, '#3E8EDE'), pch=c(3,3), cex=text.cex,
        inset=c(0.02, 0.04),
        legend=c("Uncorrected Cepheids", "Doppler-corrected Cepheids"))
    legend("bottomright", col=c(red, '#3E8EDE'), lty=c(2,4), cex=text.cex,
        inset=c(0.02, 0.04),
        legend=as.expression(c(bquote( 
            (.(sprintf("%.3f", round(a,3)))%+-%.(round(da,3)))~log[10]~P+
            (.(sprintf("%.2f", round(b,2)))%+-%.(round(db,2)))),
                 bquote(
            (.(sprintf("%.2f", round(dplr.a,2)))%+-%.(round(dplr.da,2)))~log[10]~P+
            (.(round(dplr.b,2))%+-%.(sprintf("%.2f", round(dplr.db,2))))))))
    
    magaxis(side=2:4, family=font, tcl=-0.25, labels=c(1,0,0), 
        mgp=mgp+c(0, 0.2, 0), las=1, cex.axis=text.cex, unlog='x')
    magaxis(side=1, family=font, tcl=-0.25, labels=T, mgp=mgp,#-c(0, 0.2, 0), 
        las=1, cex.axis=text.cex, unlog='x')
    
    mtext(expression("Wesenheit Index"~(I-1.45%*%(V-I))/"mag"),
        side=2, line=2, cex=text.cex)
}

make_plots(PL, 'N5584-PL', paper_pdf_height=3)



all.ceps <- read.fwf('apjaa3610t5_mrt.txt', 
    widths=c(5, 9, 9, 8, 6, 6, 3, 6, 3, 5, 5, 5, 5, 2, 3)+1, skip=40,
    col.names=c("Galaxy", "RAdeg", "DEdeg", "ID", "Period",
                "F555Wmag", "e_F555Wmag", "F814Wmag", "e_F814Wmag",
                "AmpV", "AmpI", "AmpWh", "Z", "Flag", "Ref"),
    stringsAsFactors=F)

for (galaxy in unique(all.ceps$Galaxy)) {
    DF <- all.ceps[all.ceps$Galaxy == galaxy,]
	print(1/mean((DF$e_F814Wmag/1000) / DF$F814Wmag))
}

for (galaxy in unique(all.ceps$Galaxy)) {
    DF <- all.ceps[all.ceps$Galaxy == galaxy,]
    
    if (! trimws(galaxy) %in% radial_velocities$galaxy) next 
    
    rv <- radial_velocities[radial_velocities$galaxy == trimws(galaxy),]
    shifts <- shift(DF$Period, rv$kms, rv$dkms)
    
    M <- DF$F814Wmag - 1.45*(DF$F555Wmag - DF$F814Wmag) 
    dM <- with(DF, 
        sqrt(      (e_F814Wmag/1000)**2 + 
              (1.45*e_F555Wmag/1000)**2 + 
              (1.45*e_F814Wmag/1000)**2 ) )
        #+ 
        #(10**6*0.07/(22.49*10**6*log(10)))**2  ) )
    logP <- log10(DF$Period)
    
    shifted <- DF$Period - shifts$dP
    logP.2 <- log10(shifted)
    dlogP <- shifts$sd / shifted / log(10)

    #PLR <- M ~ logP

    DF.1 <- data.frame(M, logP, logP.2, dM, dlogP)
    should.quit <- FALSE
    repeat {
        if (nrow(DF.1) < 20) {
            should.quit <- TRUE
            break
        }
        mdl.1 <- lm(M ~ logP, data=DF.1, weights=1/dM**2)
        resids <- resid(mdl.1) / DF.1$dM
        if (!(any(resids > 3))) break
        DF.1 <- DF.1[-which.max(resids),]
    }
    if (should.quit) next

    a <- summary(mdl.1)$coefficients[[2]]
    b <- summary(mdl.1)$coefficients[[1]]
    da <- summary(mdl.1)$coefficients[[4]]
    db <- summary(mdl.1)$coefficients[[3]]

    #DF.2 <- data.frame(M, logP=logP.2, dM, dlogP)
    #repeat {
    mdl.2 <- deming(M ~ logP.2, xstd=dlogP, ystd=dM, data=DF.1) 
    #    resids <- resid(mdl.2) / DF.2$dM
    #    if (!(any(resids > 4))) break
    #    DF.2 <- DF.2[-which.max(resids),]
    #}

    dplr.a <- mdl.2$coef[[2]]
    dplr.da <- sqrt(mdl.2$var[[4]])
    dplr.b <- mdl.2$coef[[1]]
    dplr.db <- sqrt(mdl.2$var[[1]])

    p.value.da <- pf(da**2/dplr.da**2, nrow(DF.1)-1, nrow(DF.1)-1)
    p.value.db <- pf(db**2/dplr.db**2, nrow(DF.1)-1, nrow(DF.1)-1)
    cat(galaxy, signif(c(a, da, b, db, dplr.a, dplr.da, dplr.b, dplr.db, 
              p.value.da, p.value.db), 5), '\\\\', '\n')
}







