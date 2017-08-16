library(magicaxis)

options(scipen=5)
options(bitmapType='cairo')

CygA.fname <- file.path('CygA', 'LOGS_MS', 'history.data')
CygB.fname <- file.path('CygB', 'LOGS_MS', 'history.data')

CygA.DF <- read.table(CygA.fname, header=1, skip=5)
CygB.DF <- read.table(CygB.fname, header=1, skip=5)

CygA.profiles <- read.table(file.path('CygA', 'LOGS_MS', 'profiles.index'),
    header=F, skip=1, col.names=c('model_number', 'priority', 'profile_number'))
CygB.profiles <- read.table(file.path('CygB', 'LOGS_MS', 'profiles.index'),
    header=F, skip=1, col.names=c('model_number', 'priority', 'profile_number'))

# clip PMS
#if (F) {
decreasing_L <- which(diff(CygA.DF$log_L) < 0 & CygA.DF$center_h1[-1] > 0.6)
if (any(decreasing_L)) {
    pms <- max(decreasing_L)
    CygA.DF <- CygA.DF[-1:-pms,]
    CygA.profiles <- CygA.profiles[-1:-pms,]
    CygA.DF$star_age <- CygA.DF$star_age - min(CygA.DF$star_age)
}

decreasing_L <- which(diff(CygB.DF$log_L) < 0 & CygB.DF$center_h1[-1] > 0.6)
if (any(decreasing_L)) {
    pms <- max(decreasing_L)
    CygB.DF <- CygB.DF[-1:-pms,]
    CygB.profiles <- CygB.profiles[-1:-pms,]
    CygB.DF$star_age <- CygB.DF$star_age - min(CygB.DF$star_age)
}
#}

CygA.all <- cbind(CygA.DF, CygA.profiles)
CygB.all <- cbind(CygB.DF, CygB.profiles)

i.max <- 30

for (ii in 0:i.max) {

CygA.rows <- 1 : (ii * nrow(CygA.all)/i.max)
CygB.rows <- 1 : (ii * nrow(CygB.all)/i.max)

CygA.rows <- c(1, which(CygA.all$star_age <= ii*max(CygA.all$star_age)/i.max))
CygB.rows <- c(1, which(CygB.all$star_age <= ii*max(CygB.all$star_age)/i.max))

CygA.DF <- CygA.all[CygA.rows,]
CygB.DF <- CygB.all[CygB.rows,]



#png(file.path('HR', paste0(ifelse(ii<10, '0', ''), ii, 'HR.png')), 
#    width=800, height=600, type='cairo', family='Times')#, res=300)
#par(mar=c(5,5,3,1), family='Times', cex=1.25, cex.lab=1.25, cex.axis=1.25)
png(file.path('HR', paste0(ifelse(ii<10, '0', ''), ii, 'HR.png')), 
    width=1260.248*3/3, height=753.973*3/2, type='cairo', family='Times', 
    res=400)
par(mar=c(3, 4, 1, 1), family='Times', #cex=1.25, cex.lab=1.25, cex.axis=1.25,
    mgp=c(2, 0.25, 0))

plot(NA, axes=F, tcl=0, 
    #type='l', lwd=2,
    xlim=rev(range(10**CygA.all$log_Teff, 10**CygB.all$log_Teff, 
        (5825-50), (5825+50))), 
    ylim=range(10**(CygA.all$log_L), 10**(CygB.all$log_L), 
        1.56-0.05, 1.56+0.05), 
    xlab=expression("Temperature"~T["eff"]/K), 
    ylab=expression("Luminosity"~L/L["Sun"]))
points(5777, 1)
points(5777, 1, pch=20, cex=0.1)

#arrows(log10(5825-50), log10(1.56), log10(5825+50), log10(1.56),
#    code=3, angle=90, length=0.02, lwd=2, lty=1)
#arrows(log10(5825), log10(1.56-0.05), log10(5825), log10(1.56+0.05),
#    code=3, angle=90, length=0.02, lwd=2, lty=1)

#rect(log10(5825-2*50), log10(1.56-2*0.05), log10(5825+2*50), log10(1.56+2*0.05),
#    col="#efefef", border=NA)
#rect(log10(5750-2*50), log10(1.27-2*0.04), log10(5750+2*50), log10(1.24+2*0.04),
#    col="#efefef", border=NA)
rect(5825-50, 1.56-0.05, 5825+50, 1.56+0.05,
    col="#d0d0d0", border=NA)
rect(5750-50, 1.27-0.04, 5750+50, 1.27+0.04,
    col="#d0d0d0", border=NA)
#arrows(log10(5750-50), log10(1.27), log10(5750+50), log10(1.27),
#    code=3, angle=90, length=0.02, lwd=2, lty=1)
#arrows(log10(5750), log10(1.27-0.04), log10(5750), log10(1.27+0.04),
#    code=3, angle=90, length=0.02, lwd=2, lty=1)

lines(10**CygB.DF$log_Teff, 10**CygB.DF$log_L, lwd=2, col="#0571b0")
points(10**CygB.DF$log_Teff[1], 10**CygB.DF$log_L[1], pch=20, cex=0.5)

lines(10**CygA.DF$log_Teff, 10**CygA.DF$log_L, lwd=2, col="#f97100")
points(10**CygA.DF$log_Teff[1], 10**CygA.DF$log_L[1], pch=20, cex=0.5)

if (ii == i.max) {
    points(10**CygA.DF$log_Teff[nrow(CygA.DF)], 
           10**CygA.DF$log_L[nrow(CygA.DF)], pch=20, cex=0.5)
    points(10**CygB.DF$log_Teff[nrow(CygB.DF)], 
           10**CygB.DF$log_L[nrow(CygB.DF)], pch=20, cex=0.5)
}
#ZAMS <- which(10**DF$log_LH / 10**DF$log_L > 0.9)[1]
#points(DF$log_Teff[ZAMS], DF$log_L[ZAMS], pch=1, cex=1)

magaxis(side=1:4, family='Times', tcl=0.25, labels=c(1,1,0,0), las=1,
    logpretty=F)#, unlog='xy')
axis(3, at=c(41000, 31000, 9500, 7240, 5920, 5300, 3850),
    labels=c("O", "B", "A", "F", "G", "K", "M"))

dev.off()



# profiles plot 
CygA.prof <- read.table(file.path('CygA', 'LOGS_MS', 
  paste0('profile', CygA.DF[nrow(CygA.DF),]$profile_number, '.data.FGONG.dat')),
    header=1)
CygB.prof <- read.table(file.path('CygB', 'LOGS_MS', 
  paste0('profile', CygB.DF[nrow(CygB.DF),]$profile_number, '.data.FGONG.dat')),
    header=1)

png(file.path('u', paste0(ifelse(ii<10, '0', ''), ii, 'u.png')), 
    width=1260.248*3/3, height=753.973*3/2, type='cairo', family='Times', 
    res=400)
par(mar=c(3, 4, 1, 1), family='Times', #cex=1.25, cex.lab=1.25, cex.axis=1.25,
    mgp=c(2, 0.25, 0))

plot(NA, axes=F, tcl=0, log='y', xaxs='i', yaxs='i',
    xlim=c(0, 1.3), 
    ylim=c(10**11, 3*10**15), #range(CygA.prof$u, CygB.prof$u), 
    ylab=expression("Sound Speed"~u), 
    xlab=expression("Radius"~r/R["☉"]))
CygA.rs <- CygA.prof$r/(6.955*10**10)
CygB.rs <- CygB.prof$r/(6.955*10**10)
lines(c(CygA.rs[1]+0.01, CygA.rs), c(CygA.prof$u[1]/1000, CygA.prof$u), 
    lwd=2, col="#f97100")
lines(c(CygB.rs[1]+0.01, CygB.rs), c(CygB.prof$u[1]/1000, CygB.prof$u), 
    lwd=2, col="#0571b0")

magaxis(side=2:4, family='Times', tcl=0.25, labels=c(1,0,0), las=1, unlog='y')
par(mgp=c(1, 0.05, 0))
magaxis(side=1, family='Times', tcl=0.25, labels=1, las=1)

dev.off()


}
