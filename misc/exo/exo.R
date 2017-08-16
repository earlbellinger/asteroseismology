#### Plot exoplanetary uncertainties as a function of stellar uncertainties 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

require(magicaxis)
source(file.path('..', '..', 'scripts', 'utils.R'))

#exo <- read.csv('uncertainties.csv', header=1)
exo <- read.csv('exo.csv', header=1)
method <- data.frame(method=exo$PLANETDISCMETH)
mass <- exo[,1:8]
mass <- cbind(mass, method)[complete.cases(mass) & 
    apply(mass, 1, function(x) all(x>0)),]
mass <- mass[-which(mass$UMASS == mass$MASS),]

radius <- exo[,-1:-8][,1:8]
radius <- cbind(radius, method)[complete.cases(radius) & 
    apply(radius, 1, function(x) all(x>0)),]

rho <- exo[,-1:-16][,1:8]
rho <- cbind(rho, method)[complete.cases(rho) & 
    apply(rho, 1, function(x) all(x>0)),]

fu.m <- with(mass, UMASS/MASS)
fu.mstar <- with(mass, UMSTAR/MSTAR)

fu.r <- with(radius, UR/R)
fu.rstar <- with(radius, URSTAR/RSTAR)

fu.rho <- with(rho, UDENSITY/DENSITY)
fu.rhostar <- with(rho, URHOSTAR/RHOSTAR)

r.lm <- lm(r~rstar, 
    data=data.frame(r=log10(fu.r*100), 
                rstar=log10(fu.rstar*100)))
dr <- 1.6230366
new.r <- 10**predict(r.lm, newdata=data.frame(rstar=log10(dr)))

m.lm <- lm(m~mstar, 
    data=data.frame(m=log10(fu.m*100), 
                mstar=log10(fu.mstar*100)))
dm <- 3.6643192
new.m <- 10**predict(m.lm, newdata=data.frame(mstar=log10(dm)))

rho.lm <- lm(rho~rhostar, 
    data=data.frame(rho=log10(fu.rho*100), 
                rhostar=log10(fu.rhostar*100)))

plot_exo <- function(add.lines=F, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
    
    par(family=font, mar=mar+c(0,0.5,0,0))
    
    plot(NA, axes=F, xaxs='i', yaxs='i', log='xy',
         xlab=expression("Stellar Uncertainty"), 
         ylab="", 
         xlim=c(0.1, 1000),#range(fu.mstar, fu.rstar),#
         ylim=c(0.1, 500))#range(fu.m, fu.r))#
    
    meth <- function(d.f) 
        as.numeric(factor(d.f$method, levels=unique(method$method))) 
    
    points(fu.mstar*100, fu.m*100, col=adjustcolor("#1E1E1E", alpha.f=0.9),
        cex=normalize(log10(mass$MASS))/3+0.33,
        pch=meth(mass))
    points(fu.rstar*100, fu.r*100, col=adjustcolor(blue, alpha.f=0.9), 
        cex=normalize(log10(radius$R))/3+0.33,
        pch=meth(radius))
    points(fu.rhostar*100, fu.rho*100, col=adjustcolor("#F97100", alpha.f=0.9),
        cex=normalize(log10(rho$DENSITY))/3+0.33,
        pch=meth(rho))
    
    #abline(m.lm, col=1, lty=2, lwd=2)
    #abline(r.lm, col="#034368", lty=2, lwd=2)
    xs <- seq(0.5, 250, 1)
    lines(xs, 10**predict(m.lm, newdata=data.frame(mstar=log10(xs))), 
        lty=2, lwd=2)
    lines(xs, 10**predict(r.lm, newdata=data.frame(rstar=log10(xs))), 
        lty=3, lwd=2, col="#034368")
    #abline(rho.lm, col="#F97100", lty=3, lwd=2)
    
    if (add.lines) {
        #segments(3.6, 0.001, 3.6, new.m, col='white', lty=2, lwd=3)
        #segments(0.01, new.m, 3.6, new.m, col='white', lty=2, lwd=3)
        
        #segments(1.7, 0.001, 1.7, new.r, col='white', lty=3, lwd=3)
        #segments(0.01, new.r, 1.7, new.r, col='white', lty=3, lwd=3)
        
        segments(dm, 0.001, dm, new.m, col=1, lty=2, lwd=2)
        segments(0.01, new.m, dm, new.m, col=1, lty=2, lwd=2)
        
        segments(dr, 0.001, dr, new.r, col="#034368", lty=3, lwd=2)
        segments(0.01, new.r, dr, new.r, col="#034368", lty=3, lwd=2)
        #
        #abline(r.lm, col=blue, lty=3, lwd=2)
    }
    
    legend("topleft", inset=c(0.02, 0.05), pch=20, #c(20, 1, 3), 
        col=c(1, blue, "#F97100"),
        legend=c('Mass', 'Radius', 'Density'), cex=0.8*text.cex) #bty='n')
    legend("bottomright", inset=c(0.02, 0.05), pch=1:6, col=1, 
        cex=0.8*text.cex,
        legend=sub("Transit Timing Variations", "TTV", unique(method$method)))
    
    magaxis(1:4, labels=F, tcl=0.25, las=1, 
            mgp=mgp, family=font, cex.axis=text.cex)
    
    xticks <- c(0.1, 1, 10, 100)
    axis(1, xticks, paste0(xticks, '%'), cex.axis=text.cex, tcl=0.25)
    axis(3, xticks, labels=F, cex.axis=text.cex, tcl=0.25)
    
    yticks <-c(0.1, 1, 10, 100, 1000)
    axis(2, yticks, paste0(yticks, '%'), cex.axis=text.cex, tcl=0.25, las=1)
    axis(4, yticks, labels=F, cex.axis=text.cex, tcl=0.25)
    
    par(mgp=mgp+c(0.8,0,0))
    title(ylab=expression("Planetary Uncertainty"))
}

make_plots(plot_exo, paste0('exo'), 
    paper_pdf_width=6.97522*1.5,
    paper_pdf_height=4.17309*1.25,
    slides=F, wide=F, tall=F)
make_plots(plot_exo, paste0('exo-lines'), 
    paper_pdf_width=6.97522*1.5,
    paper_pdf_height=4.17309*1.25,
    slides=F, wide=F, tall=F, add.lines=T)



plot_exo2 <- function(ii, add.lines=F, make.legend=F, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
    
    par(family=font, mar=mar+c(0,-2,0,0))
    
    plot(NA, axes=F, xaxs='i', yaxs='i', log='xy',
         xlab=expression("Stellar Uncertainty"), 
         ylab="", 
         xlim=c(0.1, 1000),#range(fu.mstar, fu.rstar),#
         ylim=c(0.1, 500))#range(fu.m, fu.r))#
    
    meth <- function(d.f) 
        as.numeric(factor(d.f$method, levels=unique(method$method))) 
    
    if (ii == 1) 
        points(fu.mstar*100, fu.m*100, col=adjustcolor("#323232", alpha.f=0.9),
            cex=normalize(log10(mass$MASS))+0.05,
            pch=meth(mass))
    else if (ii == 2)
        points(fu.rstar*100, fu.r*100, col=adjustcolor(blue, alpha.f=0.9), 
            cex=normalize(log10(radius$R))+0.05,#/3+0.33,
            pch=meth(radius))
    else 
        points(fu.rhostar*100, fu.rho*100, col=adjustcolor("#F97100", alpha.f=0.9),
            cex=normalize(log10(rho$DENSITY))/3+0.33,
            pch=meth(rho))
    
    #abline(m.lm, col=1, lty=2, lwd=2)
    #abline(r.lm, col="#034368", lty=2, lwd=2)
    xs <- seq(0.5, 250, 1)
    if (ii == 1)
        lines(xs, 10**predict(m.lm, newdata=data.frame(mstar=log10(xs))), 
            lty=2, lwd=2)
    else if (ii == 2)
        lines(xs, 10**predict(r.lm, newdata=data.frame(rstar=log10(xs))), 
            lty=3, lwd=2)
    #abline(rho.lm, col="#F97100", lty=3, lwd=2)
    
    if (add.lines) {
        #segments(3.6, 0.001, 3.6, new.m, col='white', lty=2, lwd=3)
        #segments(0.01, new.m, 3.6, new.m, col='white', lty=2, lwd=3)
        
        #segments(1.7, 0.001, 1.7, new.r, col='white', lty=3, lwd=3)
        #segments(0.01, new.r, 1.7, new.r, col='white', lty=3, lwd=3)
        
        if (ii == 1) {
            segments(dm, 0.001, dm, new.m, col=1, lty=2, lwd=2)
            segments(0.01, new.m, dm, new.m, col=1, lty=2, lwd=2)
        } else if (ii == 2) {
            segments(dr, 0.001, dr, new.r, col=1, lty=3, lwd=2)
            segments(0.01, new.r, dr, new.r, col=1, lty=3, lwd=2)
        }
        #
        #abline(r.lm, col=blue, lty=3, lwd=2)
    }
    
    #if (make.legend)
    #    legend("bottomright", inset=c(0.02, 0.05), pch=20, #c(20, 1, 3), 
    #        col=c(1, blue, "#F97100"),
    #        legend=c('Mass', 'Radius'), cex=0.8*text.cex) #bty='n')
    if (!make.legend)
        legend("bottomright", inset=c(0.02, 0.05), pch=1:6, col=1, 
            cex=0.8*text.cex,
            legend=sub("Transit Timing Variations", 
                "TTV", unique(method$method)))
    
    #magaxis(1:4, labels=F, tcl=0.25, las=1, 
    #        mgp=mgp, family=font, cex.axis=text.cex)
    magaxis(1:2, labels=F, tcl=0.25, las=1, 
            mgp=mgp, family=font, cex.axis=text.cex)
    
    xticks <- c(0.1, 1, 10, 100)
    axis(1, xticks, paste0(xticks, '%'), cex.axis=text.cex, tcl=0.25)
    #axis(3, xticks, labels=F, cex.axis=text.cex, tcl=0.25)
    
    yticks <-c(0.1, 1, 10, 100, 1000)
    axis(2, yticks, paste0(yticks, '%'), cex.axis=text.cex, tcl=0.25, las=1)
    #axis(4, yticks, labels=F, cex.axis=text.cex, tcl=0.25)
    
    par(mgp=mgp+c(0.8,0,0))
    title(ylab=expression("Planetary Uncertainty"))
}

make_plots(plot_exo2, paste0('exo-M'), 
    paper_pdf_width=6.97522*0.9,
    paper_pdf_height=4.17309*1.5,
    ii=1,
    slides=F, wide=F, tall=F)
make_plots(plot_exo2, paste0('exo-lines-M'), 
    paper_pdf_width=6.97522*0.9,
    paper_pdf_height=4.17309*1.5,
    ii=1,
    slides=F, wide=F, tall=F, add.lines=T)
make_plots(plot_exo2, paste0('exo-R'), 
    paper_pdf_width=6.97522*0.9,
    paper_pdf_height=4.17309*1.5,
    ii=2, make.legend=T,
    slides=F, wide=F, tall=F)
make_plots(plot_exo2, paste0('exo-lines-R'), 
    paper_pdf_width=6.97522*0.9,
    paper_pdf_height=4.17309*1.5,
    ii=2, make.legend=T,
    slides=F, wide=F, tall=F, add.lines=T)

