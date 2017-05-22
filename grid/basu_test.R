library(magicaxis)

basu.cols <- c('Model_no', 'Shells', 'Age', 'X_center', 'Y_center', 'Z_center', 
    'logL', 'logR', 'log_g', 'logT', 'm_core/M', 'm_envp/M', 'Gr_Energy', 
    'X_env', 'Z_env', 'He_core_mass', 'M_T_max', 'eta')

basu_1 <- read.table(file.path('basu_test', 
    'm1.0_fehm0.5_al1.690_zams_tams.track'), col.names=basu.cols)

basu_1.2 <- read.table(file.path('basu_test', 
    'm1.2_feh0.0_al1.690_zams_tams.track'), col.names=basu.cols)

scaling_nu_max <- function(R, M, Teff, Teff_sun=5777, nu_max_sun=3090) {
    M * nu_max_sun / ( R**2 * sqrt(Teff/Teff_sun) )
}

scaling_Delta_nu <- function(R, M, Delta_nu_Sun=135) {
    Delta_nu_Sun * M / R**3
}

filename <- file.path('test', 'M=1.2_Y=0.269838_Z=0.0165316')
if (dir.exists(filename)) {
    dir_files <- list.files(filename)
    DF <- data.frame()
    for (log_dir in file.path(filename, dir_files[grep('LOGS_', dir_files)])) {
        DF <- plyr:::rbind.fill(DF, 
            read.table(file.path(log_dir, 'history.data'), 
                header=TRUE, skip=5))
    }
    DF <- DF[order(DF$star_age),]
} else if (!file.exists(filename)) {
    exit(paste(filename), 'not found') 
} else {
    DF <- read.table(filename, header=1, skip=5)
}

# clip PMS
decreasing_L <- which(diff(DF$log_L) < 0 & DF$center_h1[-1] > 0.6)
if (any(decreasing_L)) {
    pms <- max(decreasing_L)
    print(paste("Clipping", pms, "points"))
    DF <- DF[-1:-pms,]
}

png('HR_1.2.png', width=800, height=600, type='cairo')#, family='Palatino')
par(mar=c(5,5,3,1), family='Palatino', cex=1.25, cex.lab=1.25, cex.axis=1.25)
plot(DF$log_Teff, DF$log_L, 
    type='l', axes=F, tcl=0, 
    xlim=rev(range(DF$log_Teff, basu_1.2$logT)), 
    ylim=range(DF$log_L, basu_1.2$logL), 
    xlab=expression("Temperature"~log(T["eff"]/K)), 
    ylab=expression("Luminosity"~log(L/L["Sun"])))
magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0), las=1,
    logpretty=F)#, unlog='xy')
axis(3, at=log10(c(41000, 31000, 9500, 7240, 5920, 5300, 3850)),
    labels=c("O→", "B→", "A→", "F→", "G→", "K→", "M→"))
points(log10(5777), 0)
points(log10(5777), 0, pch=20, cex=0.1)

legend("topleft", lty=1:2, col=c("black", "red"), 
    legend=c("Earl M=1.2, Y0=0.269838, Z0=0.0165316",
             "Basu M=1.2, Y0=0.275359, Z0=0.0164060"))

with(basu_1.2, lines(logT, logL, lty=2, col='red'))

dev.off()




filename <- file.path('test', 'M=1_Y=0.254246_Z=0.00558098')
if (dir.exists(filename)) {
    dir_files <- list.files(filename)
    DF <- data.frame()
    for (log_dir in file.path(filename, dir_files[grep('LOGS_', dir_files)])) {
        DF <- plyr:::rbind.fill(DF, 
            read.table(file.path(log_dir, 'history.data'), 
                header=TRUE, skip=5))
    }
    DF <- DF[order(DF$star_age),]
} else if (!file.exists(filename)) {
    exit(paste(filename), 'not found') 
} else {
    DF <- read.table(filename, header=1, skip=5)
}

# clip PMS
decreasing_L <- which(diff(DF$log_L) < 0 & DF$center_h1[-1] > 0.6)
if (any(decreasing_L)) {
    pms <- max(decreasing_L)
    print(paste("Clipping", pms, "points"))
    DF <- DF[-1:-pms,]
}

png('HR.png', width=800, height=600, type='cairo')#, family='Palatino')
par(mar=c(5,5,3,1), family='Palatino', cex=1.25, cex.lab=1.25, cex.axis=1.25)
plot(DF$log_Teff, DF$log_L, 
    type='l', axes=F, tcl=0, 
    xlim=rev(range(DF$log_Teff, basu_1$logT)), 
    ylim=range(DF$log_L, basu_1$logL),  
    xlab=expression("Temperature"~log(T["eff"]/K)), 
    ylab=expression("Luminosity"~log(L/L["Sun"])))
magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0), las=1,
    logpretty=F)#, unlog='xy')
axis(3, at=log10(c(41000, 31000, 9500, 7240, 5920, 5300, 3850)),
    labels=c("O→", "B→", "A→", "F→", "G→", "K→", "M→"))
points(log10(5777), 0)
points(log10(5777), 0, pch=20, cex=0.1)

legend("topleft", lty=1:2, col=c("black", "red"), 
    legend=c("Earl M=1, Y0=0.254246, Z0=0.00558098",
             "Basu M=1, Y0=0.277579, Z0=0.00540600"))

with(basu_1, lines(logT, logL, lty=2, col='red'))

dev.off()











filename <- file.path('test', 'M=1.2_Y=0.269838_Z=0.0165316')
if (dir.exists(filename)) {
    dir_files <- list.files(filename)
    DF <- data.frame()
    for (log_dir in file.path(filename, dir_files[grep('LOGS_', dir_files)])) {
        DF <- plyr:::rbind.fill(DF, 
            read.table(file.path(log_dir, 'history.data'), 
                header=TRUE, skip=5))
    }
    DF <- DF[order(DF$star_age),]
} else if (!file.exists(filename)) {
    exit(paste(filename), 'not found') 
} else {
    DF <- read.table(filename, header=1, skip=5)
}

# clip PMS
decreasing_L <- which(diff(DF$log_L) < 0 & DF$center_h1[-1] > 0.6)
if (any(decreasing_L)) {
    pms <- max(decreasing_L)
    print(paste("Clipping", pms, "points"))
    DF <- DF[-1:-pms,]
}
DF_1.2 <- DF


filename <- file.path('test', 'M=1_Y=0.254246_Z=0.00558098')
if (dir.exists(filename)) {
    dir_files <- list.files(filename)
    DF <- data.frame()
    for (log_dir in file.path(filename, dir_files[grep('LOGS_', dir_files)])) {
        DF <- plyr:::rbind.fill(DF, 
            read.table(file.path(log_dir, 'history.data'), 
                header=TRUE, skip=5))
    }
    DF <- DF[order(DF$star_age),]
} else if (!file.exists(filename)) {
    exit(paste(filename), 'not found') 
} else {
    DF <- read.table(filename, header=1, skip=5)
}

# clip PMS
decreasing_L <- which(diff(DF$log_L) < 0 & DF$center_h1[-1] > 0.6)
if (any(decreasing_L)) {
    pms <- max(decreasing_L)
    print(paste("Clipping", pms, "points"))
    DF <- DF[-1:-pms,]
}
DF_1 <- DF


DF_1.Dnu <- scaling_Delta_nu(10**DF_1$log_R, 1)
DF_1.2.Dnu <- scaling_Delta_nu(10**DF_1.2$log_R, 1.2)
basu_1.Dnu <- scaling_Delta_nu(10**basu_1$logR, 1)
basu_1.2.Dnu <- scaling_Delta_nu(10**basu_1.2$logR, 1.2)


png('Teff_Dnu.png', width=800, height=600, type='cairo')
par(mar=c(5,5,3,1), family='Palatino', cex=1.25, cex.lab=1.25, cex.axis=1.25)
plot(NA, 
    type='l', axes=F, tcl=0, 
    xlim=rev(range(DF_1$log_Teff, DF_1.2$log_Teff, basu_1$logT, basu_1.2$logT)),
    ylim=range(DF_1.Dnu, DF_1.2.Dnu, basu_1.Dnu, basu_1.2.Dnu), 
    xlab=expression("Temperature"~log(T["eff"]/K)), 
    ylab=expression("Scaling Large Frequency Separation"~Delta*nu/mu*Hz),
    log='y')
magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0), las=1,
    logpretty=F, unlog='y')#, unlog='xy')
axis(3, at=log10(c(41000, 31000, 9500, 7240, 5920, 5300, 3850)),
    labels=c("O→", "B→", "A→", "F→", "G→", "K→", "M→"))
points(log10(5777), 135)
points(log10(5777), 135, pch=20, cex=0.1)

legend("bottomleft", lty=1:4, lwd=2, col=c("black", "red", "black", "blue"), 
    legend=c("Earl 1 M", 
             "Basu 1 M", 
             "Earl 1.2 M", 
             "Basu 1.2 M"))

lines(basu_1.2$logT, basu_1.2.Dnu, lty=4, col='blue', lwd=2)
lines(basu_1$logT, basu_1.Dnu, lty=3, col='red', lwd=2)
lines(DF_1.2$log_Teff, DF_1.2.Dnu, lty=2, lwd=2)
lines(DF_1$log_Teff, DF_1.Dnu, lty=1, lwd=2)
#with(basu_1.2, lines(logT, logL, lty=2, col='red'))

dev.off()



DF_1.numax <- scaling_nu_max(10**DF_1$log_R, 1, 10**DF_1$log_Teff)
DF_1.2.numax <- scaling_nu_max(10**DF_1.2$log_R, 1.2, 10**DF_1.2$log_Teff)
basu_1.numax <- scaling_nu_max(10**basu_1$logR, 1, 10**basu_1$logT)
basu_1.2.numax <- scaling_nu_max(10**basu_1.2$logR, 1.2, 10**basu_1.2$logT)


png('Teff_numax.png', width=800, height=600, type='cairo')
par(mar=c(5,5,3,1), family='Palatino', cex=1.25, cex.lab=1.25, cex.axis=1.25)
plot(NA, 
    type='l', axes=F, tcl=0, 
    xlim=rev(range(DF_1$log_Teff, DF_1.2$log_Teff, basu_1$logT, basu_1.2$logT)),
    ylim=range(DF_1.numax, DF_1.2.numax, basu_1.numax, basu_1.2.numax), 
    xlab=expression("Temperature"~log(T["eff"]/K)), 
    ylab=expression("Scaling"~nu["max"]/mu*Hz),
    log='y')
magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0), las=1,
    logpretty=F, unlog='y')#, unlog='xy')
axis(3, at=log10(c(41000, 31000, 9500, 7240, 5920, 5300, 3850)),
    labels=c("O→", "B→", "A→", "F→", "G→", "K→", "M→"))
points(log10(5777), 3090)
points(log10(5777), 3090, pch=20, cex=0.1)

legend("bottomleft", lty=1:4, lwd=2, col=c("black", "red", "black", "blue"), 
    legend=c("Earl 1 M", 
             "Basu 1 M", 
             "Earl 1.2 M", 
             "Basu 1.2 M"))

lines(basu_1.2$logT, basu_1.2.numax, lty=4, col='blue', lwd=2)
lines(basu_1$logT, basu_1.numax, lty=3, col='red', lwd=2)
lines(DF_1.2$log_Teff, DF_1.2.numax, lty=2, lwd=2)
lines(DF_1$log_Teff, DF_1.numax, lty=1, lwd=2)
#with(basu_1.2, lines(logT, logL, lty=2, col='red'))

dev.off()
