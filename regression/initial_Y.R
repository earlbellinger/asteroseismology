source(file.path('..', 'scripts', 'utils.R'))


cov.dir <- file.path('learn', 'covs-simulations', 'kages')
cov.files <- list.files(cov.dir)


plot_init_Y <- function(Ys, ..., text.cex=1, 
        font=utils.font, mgp=utils.mgp, mar=utils.mar) {
    par(mar=mar+c(0.5, 0, -0.1, 0))
    
    plot(NA, axes=F, xaxs='i', yaxs='i',
        ylim=c(0, 1),
        xlim=c(0.22, 0.34),
        ylab="",
        xlab="")
    dens <- density(Ys)
    lims <- dens$x >= 0.22
    dens$x <- dens$x[lims]
    dens$y <- dens$y[lims]
    par(xpd=NA)
    polygon( c(dens$x,             rev(dens$x)),  
             c(dens$y/max(dens$y), rep(0, length(dens$x))), 
         col=adjustcolor(blue, alpha.f=0.1))
    lines(dens$x, dens$y/max(dens$y), lwd=2, col=blue)
    par(xpd=F)
    #magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3)
    magaxis(2, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3)
    
    par(mgp=mgp+c(0.6, 0, 0))
    title(xlab="Initial Helium Abundance")
    par(mgp=mgp+c(1.2, 0, 0))
    title(ylab="Density")
}

for (filename in cov.files) {
    DF <- read.table(file.path(cov.dir, filename), header=1)
    make_plots(plot_init_Y, strsplit(filename, '.dat')[[1]][1],
        filepath=file.path('plots', 'init_Y'),
        Ys=DF$Y, cex.paper=1.2,
        wide=F, tall=F, slides=F)
}



