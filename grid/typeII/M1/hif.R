#### Hydrogen ionization front calculations 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source('/scratch/seismo/bellinger/asteroseismology/scripts/utils.R') 
options(scipen=10000)
#library(Bolstad)

#int_0.r <- function(x, y) -as.numeric(cumtrapz(x, y)) # int_0^r 
#int_r.R <- function(x, y) -trapz(x, y) - int_0.r(x, y) # int_r^R

logs_dir <- 'LOGS_3MS' #'.'

prof.idx <- read.table(file.path(logs_dir, 'profiles.index'), skip=1, 
    col.names=c('mdl_num', 'priority', 'prof_num'))

prof_num <- 1
hstry <- read.table(file.path(logs_dir, 'history.data'), header=1, skip=5)

for (prof_num in prof.idx$prof_num) { #c(1,100,379,695)) {

print(prof_num)

DF <- read.table(file.path(logs_dir, 
        paste0('profile', prof_num, '.data.FGONG.dat')), 
    header=1)
#DF <- read.table(paste0('profile', prof_num, '.data'), header=1, skip=5)
hstry. <- hstry[hstry$model_number == prof.idx[prof.idx$prof_num == prof_num,]$mdl_num,]

plot_hif <- function(..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, font=utils.font, mar=utils.mar, short=F) {
    
    par(mar=mar+c(0.3, -0.5, 0.1, -0.1), lwd=1.66, las=1, cex.axis=text.cex,
        mfrow=c(1,2))
    
    xlim <- log10(c(8000, 10**3.54))
    ylim <- c(0, 4.4)
    
    plot(NA, #DF$log_Teff, DF$log_L, 
        axes=F, xaxs='i', yaxs='i', 
        type='l', lwd=3, col=1, 
        #xlim=log10(c(6000, 2500)),#log10(c(13000, 2500)),#rev(range(DF$log_Teff)), 
        #ylim=log10(c(0.2, 9000)),#range(DF$log_L), 
        xlim=xlim, 
        ylim=ylim, 
        xlab="", 
        ylab="")
    
    cols <- c(
        rgb(1, 244/255, 243/255),
        rgb(1, 229/255, 207/255),
        rgb(1, 217/255, 178/255),
        rgb(1, 199/255, 142/255),
        rgb(1, 166/255, 81/255))
    rect(xlim[1],     ylim[1], log10(7500), ylim[2], col=cols[1], border=NA)#col='#ffffbf', border=NA) # F
    rect(log10(7500), ylim[1], log10(6000), ylim[2], col=cols[2], border=NA)#col='#ffffbf', border=NA) # F
    rect(log10(6000), ylim[1], log10(5200), ylim[2], col=cols[3], border=NA) #col='#fff4e8', border=NA) # G
    rect(log10(5200), ylim[1], log10(3700), ylim[2], col=cols[4], border=NA)#col='#ffddb4', border=NA) # K
    rect(log10(3700), ylim[1], xlim[2],     ylim[2], col=cols[5], border=NA) #col='#ffbd6f', border=NA) # M
    
    gradient.rect(log10(7500+10), ylim[1], log10(7500-10), ylim[2], 
        nslices=10, reds=c(1,1), greens=c(244, 229)/255, blues=c(243, 207)/255, border=NA)
    gradient.rect(log10(6000+10), ylim[1], log10(6000-10), ylim[2], 
        nslices=10, reds=c(1,1), greens=c(229, 217)/255, blues=c(207, 178)/255, border=NA)
    gradient.rect(log10(5200+10), ylim[1], log10(5200-10), ylim[2], 
        nslices=20, reds=c(1,1), greens=c(217, 199)/255, blues=c(178, 142)/255, border=NA)
    gradient.rect(log10(3700+10), ylim[1], log10(3700-10), ylim[2], 
        nslices=50, reds=c(1,1), greens=c(199, 166)/255, blues=c(142, 81)/255, border=NA)#col='#ffddb4', border=NA) # K
    
    lines(hstry$log_Teff, hstry$log_L, type='l', lwd=2, col=1)
    points(hstry.$log_Teff, hstry.$log_L, pch=21, cex=1+hstry.$log_R, col=1, bg=blue)
    
    text(log10(10**xlim[1]+60), ylim[2]-0.55, 
        bquote(tau/Gyr==.(round(signif(hstry.$star_age/10**9, 3), 4))),
        cex=0.8*text.cex, pos=4, 
        family='Helvetica LT Std Light')
    
    text(log10(10**xlim[1]+140), ylim[2]-1.25, 
        bquote(M/M["sun"]==.(signif(hstry.$star_mass, 3))),
        cex=0.8*text.cex, pos=4, 
        family='Helvetica LT Std Light')
    
    nxticks <- 6
    nyticks <- 4
    nxminor <- 4
    nyminor <- 4
    xticks <- pretty(xlim, n=nxticks)
    yticks <- pretty(ylim, n=nyticks)
    xticks.minor <- pretty(xlim, n=nxticks*nxminor)
    yticks.minor <- pretty(ylim, n=nyticks*nyminor)
    xticks.minor <- xticks.minor[!xticks.minor %in% xticks]
    yticks.minor <- yticks.minor[!yticks.minor %in% yticks]
    par(mgp=mgp+c(0, 0.3, 0))
    #magaxis()
    xpos <- seq(1000, 8000, 1000)
    xpos2 <- seq(1000, 8000, 200)
    axis(side=1, tcl=-0.346/2, at=log10(xpos2), labels=F, lwd.ticks=par()$lwd)
    axis(side=1, tcl=-0.346, at=log10(xpos), labels=xpos, cex.axis=text.cex,
        lwd.ticks=par()$lwd)
    #axis(1, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, at=xticks,
    #    labels=as.logical(make.x))
    par(mgp=mgp+c(0, 0.43, 0))
    axis(2, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, at=yticks,
        labels=as.logical(make.y))
    #axis(1, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
    #    at=xticks.minor, labels=F)
    axis(2, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
        at=yticks.minor, labels=F)
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.x) title(xlab=expression(T["eff"]/K))
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.y) title(ylab=expression(log[10](L/L["sun"])))
    
    #magaxis(side=1:4, tcl=0, labels=F)
    spectral.divs <- log10(c(30000, 10000, 7500, 6000, 5200, 3700, 2400))
    #spectral.Teffs <- log10(c(41000, 31000, 9500, 7240, 5920, 5300, 3850))
    spectral.labs <- c("A", "F", "G", "K", "M") #c("O", "B", "A", "F", "G", "K", "M")
    spectral.Teffs <- c(log10(7800), 
        (log10(7500)+log10(6000))/2, 
        (log10(6000)+log10(5200))/2, 
        (log10(5200)+log10(3700))/2, 
        log10(3550))
    axis(3, at=spectral.divs, tcl=-0.346, labels=F, cex.axis=text.cex,
        lwd.ticks=par()$lwd)
    par(mgp=mgp+c(0, -0.1, 0))
    axis(3, at=spectral.Teffs, labels=spectral.labs, 
        cex.axis=text.cex, tcl=0)
    
    
    
    
    
    
    xlim <- c(-10, 0)#5)
    ylim <- c(0, 10)
    
    plot(NA, axes=F, 
        xaxs='i',  yaxs='i', 
        xlim=xlim, ylim=ylim,
        xlab="",   ylab="")
    
    #if(F) {
    logq <- log10(1-DF$m/(hstry.$star_mass*1.988475e33))
    ##logqs <- seq(min(logq[is.finite(logq)]), -0.033, 0.001)
    ##Ts <- splinefun(logq[is.finite(logq)], DF$T[is.finite(logq)]/10**4)(logqs)
    #Teff <- 10**hstry.$log_Teff
    #surf <- splinefun(DF$T[is.finite(logq)]/10**4, logq[is.finite(logq)])(Teff/10**4)
    
    lines(logq[logq<(-0.033)], DF$T[logq<(-0.033)]/10**4, lwd=3, col=blue)
    ##abline(h=Teff/10**4, lwd=par()$lwd, lty=2)
    ##lines(logqs, Ts, lwd=3, lty=2, col=red)
    #points(surf, Teff/10**4, pch=20, cex=0.66)
    #}
    
    if (F) {
    lines(log10(1-DF$q), DF$temperature/10**4, lwd=3, col=blue)
    
    Teff <- 10**hstry.$log_Teff
    eff.idx <- min(which(DF$temperature > Teff & is.finite(log10(1-DF$q))))
    points(log10(1-DF$q)[eff.idx], Teff/10**4, pch=20, cex=0.66)
    
    ion.idx <- min(which(DF$avg_charge_H > 0.5 & is.finite(log10(1-DF$q))))
    points(log10(1-DF$q)[ion.idx], DF$temperature[ion.idx]/10**4,
        pch=20, cex=0.66)
    
    text(log10(1-DF$q)[ion.idx], DF$temperature[ion.idx]/10**4,
        expression(HIF), 
        cex=0.8*text.cex, pos=3, 
        family='Helvetica LT Std Light')
    
    text(log10(1-DF$q)[eff.idx], Teff/10**4,
        expression(tau==2/3), 
        cex=0.8*text.cex, pos=3, 
        family='Helvetica LT Std Light')
    }
    
    #R <- 10**hstry.$log_R
    #r <- DF$r/(R*6.957e10)
    #atmos <- with(DF, 
    #    sapply(x[-1], function(r) sintegral(x[x>=r], kappa[x>=r]*rho[x>=r])$value))
    #atmos <- sapply(DF$x, function(r) sintegral(DF$x[DF$x>=r], DF$kappa[DF$x>=r]*DF$rho[DF$x>=r])$value)
    #optical_depth <- sintegral
    #optical_depth <- sintegral()
    #optical_depth <- sapply(function() int_r.R, 
    #optical_depth <- with(DF, splinefun(kappa, logq)(2/3))
    #optical_depth <- 
    #temp_pt <- with(DF, splinefun(logq, T/10**4)(optical_depth))
    
    
    #points(optical_depth, temp_pt, pch=20, col=blue, cex=1)
    
    #abline(h=0, lty=2, lwd=par()$lwd)
    
    nxticks <- 6
    nyticks <- 4
    nxminor <- 4
    nyminor <- 4
    xticks <- pretty(xlim, n=nxticks)
    yticks <- pretty(ylim, n=nyticks)
    xticks.minor <- pretty(xlim, n=nxticks*nxminor)
    yticks.minor <- pretty(ylim, n=nyticks*nyminor)
    xticks.minor <- xticks.minor[!xticks.minor %in% xticks]
    yticks.minor <- yticks.minor[!yticks.minor %in% yticks]
    par(mgp=mgp+c(0, 0.3, 0))
    axis(1, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, at=xticks,
        labels=as.logical(make.x))
    par(mgp=mgp+c(0, 0.43, 0))
    axis(2, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, at=yticks,
        labels=as.logical(make.y))
    axis(1, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
        at=xticks.minor, labels=F)
    axis(2, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
        at=yticks.minor, labels=F)
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.x) title(xlab=expression(Fractional~mass~log[10](1-m/M)))
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.y) title(ylab=expression(Temp.~T/(10^4~K))) 
    
    
    #par(new=T)
    #plot(NA, axes=F, xlim=c(-10, 0), ylim=c(0, 1), xaxs='i', yaxs='i',
    #    xlab="", ylab="")
    #lines(logq[logq<(-0.033)], DF$kappa[logq<(-0.033)], lwd=3, col=blue)
}

make_plots(plot_hif, paste0('hif', sprintf("%06d", prof_num)),
        filepath=file.path('plots'),
        cex.paper=0.93, 
        wide=T, thin=F, tall=F, slides=F, make_png=T, make_pdf=F,
        make.x=T,
        make.y=T,
        font="Palatino Linotype", 
        use.cairo=T)

}

