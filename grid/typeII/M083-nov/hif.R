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

logs_dir <- 'LOGS_3MS' #'M09-nov'#'.'#

prof.idx <- read.table(file.path(logs_dir, 'profiles.index'), skip=1, 
    col.names=c('mdl_num', 'priority', 'prof_num'))

prof_num <- 1628#1100#695
hstry <- read.table(file.path(logs_dir, 'history.data'), header=1, skip=5)

for (prof_num in prof.idx$prof_num) { #c(1,100,379,695)) {#

#DF <- read.table(file.path(logs_dir, 
#        paste0('profile', prof_num, '.data.FGONG.dat')), 
#    header=1)
DF <- read.table(file.path(logs_dir, paste0('profile', prof_num, '.data')), header=1, skip=5)
hstry. <- hstry[hstry$model_number == prof.idx[prof.idx$prof_num == prof_num,]$mdl_num,]

plot_hif <- function(..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, font=utils.font, mar=utils.mar, short=F) {
    
    par(mar=mar+c(0.3, -0.5, 2, -0.1), lwd=1.66, las=1, cex.axis=text.cex,
        mfrow=c(1,2))
    
    xlim <- log10(c(15000, 3000))
    ylim <- c(0, 4.4)
    
    plot(NA, 
        axes=F, xaxs='i', yaxs='i', 
        type='l', lwd=3, col=1, 
        xlim=xlim, 
        ylim=ylim, 
        xlab="", 
        ylab="")
    
    spectral.divs <- log10(c(30000, 10000, 7500, 6000, 5200, 3700, 2400))
    rs <- c(175/255, 199/255, 1, 1, 1, 1, 1, 1)
    gs <- c(201, 216, 244, 229, 217, 199, 166)/255
    bs <- c(1, 1, 243/255, 207/255, 178/255, 142/255, 81/255)
    cols <- c(
        rgb(175/255, 201/255, 1),       # O
        rgb(199/255, 216/255, 1),       # B
        rgb(1,       244/255, 243/255), # A 
        rgb(1,       229/255, 207/255), # F 
        rgb(1,       217/255, 178/255), # G 
        rgb(1,       199/255, 142/255), # K 
        rgb(1,       166/255, 81/255))  # M
    for (ii in 1:length(spectral.divs)) {
        div <- spectral.divs[ii]
        if (div > xlim[1]) next 
        if (div < xlim[2]) div <- xlim[2]
        if (ii == 1) {
            #if (xlim[1] > div) {
                rect(xlim[1], ylim[1], div, ylim[2], col=cols[ii], border=NA)
            #}
        } else {
            prev <- spectral.divs[ii-1]
            if (prev > xlim[1]) prev <- xlim[1]
            rect(prev, ylim[1], div, ylim[2], col=cols[ii], border=NA)
        }
    }
    for (ii in 2:(length(spectral.divs)-1)) {
        div <- spectral.divs[ii]
        gradient.rect(div+0.0025, ylim[1], div-0.0025, ylim[2],
            nslices=10, border=NA, 
            reds=c(rs[ii], rs[ii+1]), 
            greens=c(gs[ii], gs[ii+1]),
            blues=c(bs[ii], bs[ii+1]))
    }
    #rect(xlim[1],     ylim[1], log10(7500), ylim[2], col=cols[1], border=NA)#col='#ffffbf', border=NA) # F
    #rect(log10(7500), ylim[1], log10(6000), ylim[2], col=cols[2], border=NA)#col='#ffffbf', border=NA) # F
    #rect(log10(6000), ylim[1], log10(5200), ylim[2], col=cols[3], border=NA) #col='#fff4e8', border=NA) # G
    #rect(log10(5200), ylim[1], log10(3700), ylim[2], col=cols[4], border=NA)#col='#ffddb4', border=NA) # K
    #rect(log10(3700), ylim[1], xlim[2],     ylim[2], col=cols[5], border=NA) #col='#ffbd6f', border=NA) # M
    
    if (F) {
    gradient.rect(log10(7500+10), ylim[1], log10(7500-10), ylim[2], 
        nslices=10, reds=c(1,1), greens=c(244, 229)/255, blues=c(243, 207)/255, border=NA)
    gradient.rect(log10(6000+10), ylim[1], log10(6000-10), ylim[2], 
        nslices=10, reds=c(1,1), greens=c(229, 217)/255, blues=c(207, 178)/255, border=NA)
    gradient.rect(log10(5200+10), ylim[1], log10(5200-10), ylim[2], 
        nslices=20, reds=c(1,1), greens=c(217, 199)/255, blues=c(178, 142)/255, border=NA)
    gradient.rect(log10(3700+10), ylim[1], log10(3700-10), ylim[2], 
        nslices=50, reds=c(1,1), greens=c(199, 166)/255, blues=c(142, 81)/255, border=NA)
    }
    
    lines(hstry$log_Teff, hstry$log_L, type='l', lwd=2, col=1)
    
    mdl.col <- cols[which.min(hstry.$log_Teff < spectral.divs)]
    mdl.cex <- 1+hstry.$log_R
    #points(hstry.$log_Teff, hstry.$log_L, pch=21, cex=mdl.cex, 
    #    col=1, bg=mdl.col, lwd=par()$lwd)
    points(hstry.$log_Teff, hstry.$log_L, pch=21, cex=mdl.cex, 
        col="#FFFFFF", lwd=par()$lwd,
        bg=with(DF[nrow(DF),], 
            rgb(red=x*.8, green=y*.8, blue=z*.8)))
    
    #in.cz <- F
    #for (ii in 1:nrow(DF)) {
    #    if (DF$mixing_type[ii]) {
    #        if (in.cz) next
    #        logr.outer <- 1+log10(DF$radius[ii]/10**hstry.$log_R)
    #        in.cz <- T
    #    } else {
    #        if (in.cz) {
    #            logr.inner <- 1+log10(DF$radius[ii]/10**hstry.$log_R)
    #            points(hstry.$log_Teff, hstry.$log_L, pch=20, col='gray',
    #                cex=logr.outer/mdl.cex)
    #            points(hstry.$log_Teff, hstry.$log_L, pch=20, col=mdl.col,
    #                cex=logr.inner/mdl.cex)
    #        }
    #        in.cz <- F
    #    }
    #}
    ##conv.env <- DF$gradr_sub_grada[1] > 0
    ##if (conv.env) {
    ##    points()
    ##}
    
    
    if(F) {
    text(log10(10**xlim[1]+60), ylim[2]-0.55, 
        bquote(tau/Gyr==.(round(signif(hstry.$star_age/10**9, 3), 4))),
        cex=0.8*text.cex, pos=4, 
        family='Helvetica')
    
    text(log10(10**xlim[1]+140), ylim[2]-1.25, 
        bquote(M/M["sun"]==.(signif(hstry.$star_mass, 3))),
        cex=0.8*text.cex, pos=4, 
        family='Helvetica')
    }
    
    age <- round(signif(hstry.$star_age/10**9, 4), 5)
    mass <- signif(hstry.$star_mass, 3)
    par(family="Helvetica")
    legend('topleft', #lty=NA, pch=NA, 
        cex=0.8*text.cex,
        bty='n', inset=c(-0.05, 0),
        #bg="#FFFFFF",
        legend=c(as.expression(bquote(tau/Gyr==.(age))),
                 as.expression(bquote(M/M["sun"]==.(mass)))))
    par(family=font)
    
    
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
    xpos <- seq(10**xlim[2], 10**xlim[1], 1000)
    xpos2 <- seq(10**xlim[2], 10**xlim[1], 200)
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
    
    
    #magaxis(side=1, lwd.ticks=par()$lwd, tcl=-0.346, mgp=mgp+c(0, 0.3, 0), unlog='x')
    #magaxis(side=2, lwd.ticks=par()$lwd, tcl=-0.346, mgp=mgp+c(0, 0.43, 0), unlog='y')
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.x) title(xlab=expression(T["eff"]/K))
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.y) title(ylab=expression(log[10](L/L["sun"])))
    if (make.x) mtext(expression("Spectral Type"), side=3, line=1.5, cex=text.cex)
    
    #magaxis(side=1:4, tcl=0, labels=F)
    #spectral.Teffs <- log10(c(41000, 31000, 9500, 7240, 5920, 5300, 3850))
    #spectral.labs <- c("A", "F", "G", "K", "M") #c("O", "B", "A", "F", "G", "K", "M")
    #spectral.Teffs <- c(log10(7800), 
    #    (log10(7500)+log10(6000))/2, 
    #    (log10(6000)+log10(5200))/2, 
    #    (log10(5200)+log10(3700))/2, 
    #    log10(3550))
    spectral.labs <- c("O", "B", "A", "F", "G", "K", "M")
    selector <- 1:(length(spectral.divs))
    spectral.Teffs <- sapply(selector, 
        function(ii) {
            #if (ii == 1)
            div <- spectral.divs[ii]
            if (div > xlim[1]) return(Inf)
            if (div < xlim[2]) div <- xlim[2]
            if (ii == 1) return((xlim[1]+div)/2)
            prev <- spectral.divs[ii-1]
            if (prev > xlim[1]) prev <- xlim[1]
            #sum(spectral.divs[(ii-1):(ii)])/2
            (div+prev)/2
        })
    axis(3, at=spectral.divs, tcl=-0.346, labels=F, cex.axis=text.cex,
        lwd.ticks=par()$lwd)
    par(mgp=mgp+c(0, -0.1, 0))
    axis(3, at=spectral.Teffs, labels=spectral.labs[selector], 
        cex.axis=text.cex, tcl=0)
    
    
    
    
    
    
    xlim <- c(-10, 0)#5)
    ylim <- c(0, 10)
    
    plot(NA, axes=F, 
        xaxs='i',  yaxs='i', 
        xlim=xlim, ylim=ylim,
        xlab="",   ylab="")
    
    if(F) {
    logq <- log10(1-DF$m/(hstry.$star_mass*1.988475e33))
    lines(logq[logq<(-0.033)], DF$T[logq<(-0.033)]/10**4, lwd=3, col=blue)
    }
    
    lines(log10(1-DF$q), DF$neutral_fraction_H, lwd=3)
    lines(log10(1-DF$q), DF$temperature/10**4, lwd=3, col=blue)
    lines(log10(1-DF$q), DF$opacity/10**1, lwd=3, col=orange)
    
    #Teff <- 10**hstry.$log_Teff
    #eff.idx <- min(which(DF$temperature > Teff & is.finite(log10(1-DF$q))))
    #points(log10(1-DF$q)[eff.idx], Teff/10**4, pch=20, cex=0.66)
    
    #ion.idx <- min(which(DF$neutral_fraction_H < 0.5 & is.finite(log10(1-DF$q))))
    #ion.idx.min <- max(which(DF$neutral_fraction_H >= 0.01 & is.finite(log10(1-DF$q))))
    #ion.idx.max <- min(which(DF$neutral_fraction_H <= 0.99 & is.finite(log10(1-DF$q))))
    
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
    axis(1, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
        at=xticks.minor, labels=F)
    
    par(mgp=mgp+c(0, 0.43, 0))
    axis(2, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, at=yticks,
        labels=as.logical(make.y))
    axis(2, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
        at=yticks.minor, labels=F)
    
    par(mgp=mgp+c(0, 0.3, 0))
    axis(3, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, 
        at=log10(1-DF$q)[sapply(xticks, function(x) 
            which.min(abs(x-log10(1-DF$radius/10**hstry.$log_R))))],
        labels=xticks)
    axis(3, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
        at=log10(1-DF$q)[sapply(xticks.minor, function(x) 
            which.min(abs(x-log10(1-DF$radius/10**hstry.$log_R))))], 
        labels=F)
    
    #par(mgp=mgp+c(0, 0.43, 0))
    #axis(4, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, at=yticks,
    #    labels=as.logical(make.y))
    #axis(4, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
    #    at=yticks.minor, labels=F)
    
    par(family="Helvetica")
    legend('topleft', col=c(orange, blue, 1), lwd=3, cex=0.8*text.cex,
        bty='n', lty=c(1,1,1),
        #bg="#FFFFFF",
        legend=c(expression(Opacity~kappa/(10~cm^2/g)),
                 expression(Temp.~T/(10^4~K)),
                 expression(Neutral~fraction~H)
                 #expression(HIF)
            ))
    par(family=font)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.x) title(xlab=expression(Fractional~mass~log[10](1-m/M)))
    
    par(mgp=mgp+c(0.55, 0, 0))
    #if (make.y) title(ylab=expression(Temp.~T/(10^4~K))) 
    if (make.y) title(ylab=expression(Envelope~structure)) 
    
    if (make.x) mtext(expression(Fractional~radius~log[10](1-r/R)), side=3, line=1.5, cex=text.cex)
    #if (make.y) mtext(expression(Opacity~kappa/(10^2~cm^2/g)), side=4, line=1.5, cex=text.cex, las=0) 
    #par(new=T)
    #plot(NA, axes=F, xlim=c(-10, 0), ylim=c(0, 1), xaxs='i', yaxs='i',
    #    xlab="", ylab="")
    #lines(logq[logq<(-0.033)], DF$kappa[logq<(-0.033)], lwd=3, col=blue)
}
plot_hif()

make_plots(plot_hif, paste0('hif', sprintf("%06d", prof_num)),
        filepath=file.path('plots'),
        cex.paper=0.93, 
        wide=T, short=F, thin=F, slides=F, make_png=T, make_pdf=F,
        make.x=T,
        make.y=T,
        font="Palatino Linotype", 
        use.cairo=T)

}

