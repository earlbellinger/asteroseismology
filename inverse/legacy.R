#### Plot Kepler Legacy CDFs for KASC proceedings 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

require(magicaxis)
source(file.path('..', 'scripts', 'utils.R'))

#colors <- rev(c('#588C7E', '#F2E394', '#F2AE72', '#D96459', '#8C4646'))
colors <- c(blue, "black", red)
col.pal <- adjustcolor(colorRampPalette(colors)(7))

curr <- read.table(file.path('learn', 'tables-simulations', 
                             'legacy_curr_ascii.dat'), header=1)
init <- read.table(file.path('learn', 'tables-simulations', 
                             'legacy_init_ascii.dat'), header=1)
curr.names <- expression(tau, X["c"], log~g, L, R, Y["surf"])
init.names <- expression(M, Y[0], Z[0], alpha["MLT"], alpha["ov"], D)

mar <- utils.mar + c(0,-1,0,2.2)

exclude <- c(5774694) # fake Sun 
pred.dir <- file.path('learn', 'covs-simulations', 'legacy')
for (pred.fname in list.files(pred.dir)) {
    pred.data <- read.table(file.path(pred.dir, pred.fname), header=1)
    if (any(pred.data$X_c < 0.01)) 
        exclude <- c(exclude, strsplit(pred.fname, '.dat')[[1]])
}
print(exclude)

plot_cdfs <- function(df.i, inset=-.29, ..., 
                      text.cex=1, mgp=utils.mgp, mar=mar, font="Times") {
    par(family="Times")
    DF <- if (df.i==1) curr else init
    DF <- DF[!DF$KIC %in% exclude,]
    group <- names(DF[,-1])[-grep('^d_', names(DF[,-1]))] 
    
    unc <- DF[paste0('d_', group)] / DF[group] * 100
    sorted <- order(sapply(unc, max))
    unc <- unc[,sorted]
    
    print(sapply(unc, mean))
    print(sapply(unc, fivenum))
    
    xmin <- if (df.i==1) 0.1 else 1 
    xmax <- if (df.i==1) 100 else 200 
    xlim <- c(xmin, xmax) 
    
    xticks <- if (df.i==1) c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100) else
                           c(1, 2, 5, 10, 20, 50, 100, 200)
    
    plot(NA, axes=F, xaxs='i', yaxs='i', log='x',
         xlab=expression("Relative uncertainty"),
         ylab=expression("Cumulative count of stars"),
         xlim=xlim,
         ylim=c(0, 1.05))
    
    grid(0, NULL, lty = 6, col="cornsilk2") 
    abline(v=xticks, lty=6, col = "cornsilk2")
    
    for (ii in 1:length(group)) lines(ecdf(unc[,ii]), col=col.pal[ii], cex=0.33)
    
    magaxis(1:4, labels=F, tcl=-0.25, las=1, 
            mgp=mgp, family=font, cex.axis=text.cex)
    axis(1, xticks, paste0(xticks, '%'), cex.axis=text.cex, tcl=0)
    oh.one <- pretty(c(0, 1))
    axis(2, oh.one, round(seq(0, nrow(unc), length.out=length(oh.one))), 1,
         cex.axis=text.cex, tcl=0, las=1, mgp=mgp+c(0, 0.25, 0))
    
    par(xpd=T)
    legend("right", 
           if (df.i==1) curr.names[sorted] else init.names[sorted],
           inset=c(inset,0),
           bty='n', pch=20, lty=F, col=col.pal, cex=text.cex)
    par(xpd=F)
}

for (df.i in 1:2) 
    make_plots(plot_cdfs, 
               paste0('cdf', '-', if (df.i==1) 'curr' else 'init'),
               mar=mar,
               df.i=df.i)#,
               #inset = if (df.i==1) -.3 else -.4)

              
