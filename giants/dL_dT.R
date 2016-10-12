library(magicaxis)

png('knee_finder.png', width=800, height=600)

directories <- list.dirs('simulations-old', recursive=F)
par(mar=c(5, 5, 1, 1), cex=1.25, cex.lab=1.25, cex.axis=1.25, 
    family='Palatino')
plot(NA, axes=F, tcl=0,
    xlab=expression("Temperature"~log(T["eff"]/K)),
    ylab=expression("Luminosity"~log(L/L["☉"])),
    xlim=c(4.3, log10(5000)), 
    ylim=c(-0.1, 3), family='Palatino',
    cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)
magaxis(1:4, tcl=0.25, labels=c(1,1,0,0), family='Palatino', las=1, 
    cex=1.25, cex.lab=1.25, cex.axis=1.25)
for (directory in directories) {
    DF <- read.table(file.path(directory, 'history-ms.data'),
        skip=5, header=1)
    lines(DF$log_Teff, DF$log_L, lty=2)
    
    DF <- read.table(file.path(directory, 'LOGS', 'history.data'),
        skip=5, header=1)
    tip <- which(DF$log_L == max(DF$log_L))
    DF <- DF[1:tip,]
    
    lines(DF$log_Teff, DF$log_L)
    
    DF <- DF[DF$he_core_mass>0 & DF$star_age <= 13.8*10**9,]
    DF <- DF[DF$log_Teff >= 3.675 & DF$log_Teff <= 3.74,]
    if (nrow(DF) == 0) next
    
    new_Ts <- seq(min(DF$log_Teff), max(DF$log_Teff), 0.0001)
    new_Ls <- splinefun(DF$log_Teff, DF$log_L)(new_Ts)
    dlogL.dlogT <- splinefun(DF$log_Teff, DF$log_L)(new_Ts, 1)
    separator <- which(dlogL.dlogT < -3)
    separator <- separator[length(separator)]
    
    points(new_Ts[separator], new_Ls[separator], 
        col='red', pch=20, cex=1)
}

dev.off()


png('knee_cropper.png', width=800, height=600)

directories <- list.dirs('simulations', recursive=F)
par(mar=c(5, 5, 1, 1), cex=1.25, cex.lab=1.25, cex.axis=1.25, 
    family='Palatino')
plot(NA, axes=F, tcl=0,
    xlab=expression("Temperature"~log(T["eff"]/K)),
    ylab=expression("Luminosity"~log(L/L["☉"])),
    xlim=c(4.3, log10(5000)), 
    ylim=c(-0.1, 3), family='Palatino',
    cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25)
magaxis(1:4, tcl=0.25, labels=c(1,1,0,0), family='Palatino', las=1, 
    cex=1.25, cex.lab=1.25, cex.axis=1.25)
for (directory in directories) {
    DF <- read.table(file.path(directory, 'LOGS', 'history.data'),
        skip=5, header=1)
    tip <- which(DF$log_L == max(DF$log_L))
    DF <- DF[1:tip,]
    
    lines(DF$log_Teff, DF$log_L)
    
    DF <- DF[DF$he_core_mass>0 & DF$star_age <= 13.8*10**9,]
    DF <- DF[DF$log_Teff >= 3.675 & DF$log_Teff <= 3.74,]
    if (nrow(DF) == 0) next
    
    new_Ts <- seq(min(DF$log_Teff), max(DF$log_Teff), 0.0001)
    new_Ls <- splinefun(DF$log_Teff, DF$log_L)(new_Ts)
    dlogL.dlogT <- splinefun(DF$log_Teff, DF$log_L)(new_Ts, 1)
    separator <- which(dlogL.dlogT < -3)
    separator <- separator[length(separator)]
    
    #points(new_Ts[separator], new_Ls[separator], 
    #    col='red', pch=20, cex=1)
}

dev.off()

