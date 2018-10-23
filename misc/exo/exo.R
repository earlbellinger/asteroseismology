#### Plot exoplanetary uncertainties as a function of stellar uncertainties 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

require(magicaxis)
source(file.path('..', '..', 'scripts', 'utils.R'))

set.seed(0)
exo <- read.csv('exo.csv', header=1)
exo <- exo[exo$PLANETDISCMETH=='Transit',]
exo <- exo[sample(1:nrow(exo)),]
exo <- exo[-which(with(exo, UR/R>1 | URSTAR/RSTAR > 1)),]
exo <- exo[-which(with(exo, UR/R<0.004 | URSTAR/RSTAR < 0.004)),]
#method <- data.frame(method=exo$PLANETDISCMETH)

radius <- exo[,-1:-8][,1:8]
radius <- radius[complete.cases(radius) & 
    apply(radius, 1, function(x) all(x>0)),]

fu.r <- with(radius, UR/R)
fu.rstar <- with(radius, URSTAR/RSTAR)

r.lm <- lm(r~rstar, 
    data=data.frame(r=log10(fu.r*100), 
                rstar=log10(fu.rstar*100)))
#dr <- 1.6230366
#new.r <- 10**predict(r.lm, newdata=data.frame(rstar=log10(dr)))

plot_exo <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.1, -0.2, 0.1, 0.3), lwd=1.5, las=1, cex.axis=text.cex)
    
    xlim <- c(0.4, 100)
    ylim <- c(0.4, 100)
    
    plot(NA, axes=F, log='xy',
        xaxs='i',  yaxs='i', 
        xlab="",   ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- 10**pretty(log10(xlim), n=8)
    xticks <- xticks[c(-1, -length(xticks))]
    xminors <- 10**pretty(log10(xlim), n=30)
    xminors <- xminors[c(-1, -length(xminors))]
    
    if (F) {
    for (xtick in xticks)
        for (ytick in xminors) 
            points(xtick, ytick, pch=19, cex=0.1, lwd=0.5, #lwd=0, 
                bg='darkgray', col='darkgray')
    
    for (ytick in xticks) 
        for (xtick in xminors) 
            points(xtick, ytick, pch=19, cex=0.1, lwd=0.5, #lwd=0, 
                bg='darkgray', col='darkgray')
    }
    #for (ytick in xminors)
    #    points(1, ytick, pch=19, cex=0.1, lwd=0.5, 
    #        bg='black', col=adjustcolor(1, alpha.f=0.5))
    
    #rect(1.85, 0.18, 100, 0.22, col='white', border=NA)
    rect(6, 0.4, 100, 0.7, col='white', border=NA)
    
    rect(0.6977357, 0.01, 2.2955524, 1000, 
        col=adjustcolor(orange, alpha.f=0.2), border=NA)
    
    #rect(1.061222, 0.01, 1.793358, 1000, 
    #    col=adjustcolor(1, alpha.f=0.2), border=NA)
    
    #abline(a=0, b=1, lwd=1.5, lty=2)
    xs <- seq(xlim[1], xlim[2], 10)
    lines(xs, xs, lty=2, lwd=1.5)
    
    points(fu.rstar*100, fu.r*100, 
        #bg=adjustcolor('darkgreen', alpha.f=0.75), col='white',
        bg=adjustcolor(red, alpha.f=0.75), col='white',
        pch=21, cex=0.8, lwd=0.5)
    
    xs <- seq(0.15, 90, 1)
    #lines(xs, 10**predict(r.lm, newdata=data.frame(rstar=log10(xs))), 
    #    lty=3, lwd=4) #2.5)
    
    #arrows(3.15, 0.194, 5.4, 0.194, length=0.05, angle=30, code=1, lwd=0.75)
    arrows(3.1, 0.57, 5.8, 0.57, length=0.05, angle=30, code=1, lwd=0.75)
    par(family="Helvetica")
    text(5.7, 0.554,#0.192, 
        labels='asteroseismology',
        pos=4, cex=0.8*text.cex, col='black')
    par(family=font)
    #abline(r.lm, lwd=4, lty=3)
    
    
    
    par(xpd=NA)
    rect(xlim[2], ylim[1]*0.1, xlim[2]*1.1, ylim[2]*10, col='white', border=NA)
    rect(xlim[1]*0.9, ylim[1], xlim[2]*1.1, ylim[1]*0.9, col='white', border=NA)
    par(xpd=F)
    
    
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=1.3*text.cex, 
        family=font, majorn=4, labels=F, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=4, labels=F, lwd.ticks=par()$lwd)
    
    xticks <- c(0.1, 1, 10, 100)
    axis(1, xticks, paste0(xticks, '%'), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    yticks <-c(0.1, 1, 10, 100, 1000)
    axis(2, yticks, paste0(yticks, '%'), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    mtext('Exoplanet Radius Uncertainty', 2, 2.5, outer=F, las=0, cex=text.cex)
    mtext('Stellar Radius Uncertainty', 1, 1.7, outer=F, cex=text.cex)
}

make_plots(plot_exo, paste0('exo-R'), 
    #paper_pdf_width=6.97522,
    #paper_pdf_height=4.17309*1.5,
    paper_pdf_width=4.17309*1.385,
    paper_pdf_height=4.17309*1.31,
    cex.paper=0.75,
    slides=F, wide=F, tall=F, make_png=F)
