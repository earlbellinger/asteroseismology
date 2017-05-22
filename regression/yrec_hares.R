
true <- read.table('data/yrec-hares/yrec_hares.dat', header=1)
unc.0 <- read.table('yrec_hares/0.dat', header=1)
unc.1 <- read.table('yrec_hares/1.dat', header=1)
unc.2 <- read.table('yrec_hares/2.dat', header=1)
unc.3 <- read.table('yrec_hares/3.dat', header=1)
uncs <- list(unc.0, unc.1, unc.2, unc.3)


for (out in names(unc.0)) {
    
    ylim <- range(unc.0[[out]], unc.1[[out]], unc.2[[out]], unc.3[[out]])
    ylim <- c(max(-2, ylim[1]), min(2, ylim[2]))
    
    plot(NA, 
        xlim=range(true[[out]]), 
        ylim=ylim,
        ylab="Relative difference")
    
    for (unc.i in 1:length(uncs)) {
        unc <- uncs[[unc.i]]
        plx <- predict(loess(unc[[out]] ~ true[[out]]), se=T)
        lines(unique(true[[out]]), unique(plx$fit), col=unc.i)
        lines(unique(true[[out]]), 
            unique(plx$fit) - unique(qt(0.975,plx$df)*plx$se), lty=2, col=unc.i)
        lines(unique(true[[out]]), 
            unique(plx$fit) + unique(qt(0.975,plx$df)*plx$se), lty=2, col=unc.i)
    }
    
    #lowess(true[[out]], unc.0[[out]])
    
    
    smooth.spline(true[[out]], unc.0[[out]])
    
    points(true[[out]], unc.3[[out]], pch=4, col=4)
    points(true[[out]], unc.2[[out]], pch=3, col=3)
    points(true[[out]], unc.1[[out]], pch=2, col=2)
    points(true[[out]], unc.0[[out]])
    dev.off()
}


