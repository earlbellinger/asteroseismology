#### Plot overlapping histograms of 16 Cyg A & B 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)

cygA <- read.table(file.path('learn_covs-simulations-amp', 
    'perturb', '16CygAamp.dat'), header=1)
cygB <- read.table(file.path('learn_covs-simulations-amp', 
    'perturb', '16CygBamp.dat'), header=1)

metcalfe <- read.table(file.path('data', 'cyg.dat'), header=1)

measA <- read.table(file.path('data', '16CygA-obs.dat'), header=1)
measB <- read.table(file.path('data', '16CygB-obs.dat'), header=1)

vermA <- data.frame(name="Y_surf", value=0.241, uncertainty=0.01)
vermB <- data.frame(name="Y_surf", value=0.242, uncertainty=0.024)

measA <- rbind(measA, vermA)
measB <- rbind(measB, vermB)

plot_cygs <- function(name, cygA, cygB, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font) {
    A <- density(cygA[[name]])
    B <- density(cygB[[name]])
    
    xlim <- range(A$x, B$x)
    ylim <- range(0, A$y, B$y)
    
    has_other <- F
    if (name %in% names(metcalfe)) {
        has_other <- T
        means <- metcalfe[[name]]
        stds <- metcalfe[[paste0('d', name)]]
        #other_A <- density(rnorm(100000, means[1], stds[1]))
        #other_B <- density(rnorm(100000, means[2], stds[2]))
        
        #xlim <- range(xlim, other_A$x, other_B$x)
        xlim <- range(xlim, means-2*stds, means+2*stds)#other_A$x, other_B$x)
        other_name <- "AMP"
        other_citation <- "Metcalfe et al. 2015"
    }
    
    if (name %in% measA$name) {
        has_other <- T
        means <- c(measA[measA$name == name,]$value, 
                   measB[measB$name == name,]$value)
        stds <- c(measA[measA$name == name,]$uncertainty, 
                  measB[measB$name == name,]$uncertainty)
        #other_A <- density(rnorm(100000, means[1], stds[1]))
        #other_B <- density(rnorm(100000, means[2], stds[2]))
        xlim <- range(xlim, means-2*stds, means+2*stds)#other_A$x, other_B$x)
        other_name <- if (name %in% vermA$name) 
            "ast" else "int"
        other_citation <- if (name %in% vermA$name) 
            "Verma et al. 2014" else "White et al. 2013"
    }
    
    show_legend <- (name == "L" || name == "radius" || name == "Y_surf" || 
           name == "age") && has_other
    
    xmean <- mean(xlim)
    
    if (name == "Y_surf") {
        xlim[2] <- xlim[2] + xmean*0.15
    }
    
    if (name == "L") {
        xlim[1] <- xlim[1] - xmean*0.2
        xlim[2] <- xlim[2] + xmean*0.1
    }
    
    if (show_legend) {
        xlim[1] <- xlim[1] - xmean * 0.2
        xlim[2] <- xlim[2] + xmean * 0.1
    }
    
    par(mar=c(2.5, 1, 1, 1), mgp=mgp-c(0.75,0,0))
    plot(A, axes=F, col=red, lwd=1.5, yaxs='i', xaxs='i', lty=2,
        xlim=xlim, ylim=c(0, ylim[2]*1.01), 
        xlab="", ylab="", main="")
    magaxis(side=1, family=font, tcl=0.25, labels=1, 
            las=1, mgp=mgp, cex.axis=text.cex)
    lines(B, col=blue, lwd=1.5, lty=2)
    
    if (has_other) {
        apos <- max(ylim) * 0.2#quantile(A$y, .55)
        bpos <- max(ylim) * 0.3#quantile(B$y, .6)
        arrows(means[1]-2*stds[1], apos, means[1]+2*stds[1], apos, 
            code=3, col=red, length=0.1)
        arrows(means[2]-2*stds[2], bpos, means[2]+2*stds[2], bpos,
            code=3, col=blue, length=0.1)
    }
    par(xpd=NA)
    eps_A <- if (mean(cygA[[name]]) > 0)
        signif(sqrt(var(cygA[[name]])) / mean(cygA[[name]]) * 100, 3)
        else 0
    eps_B <- if (mean(cygB[[name]]) > 0)
        signif(sqrt(var(cygB[[name]])) / mean(cygB[[name]]) * 100, 3)
        else 0
    uncertainties <- c(as.expression(bquote(epsilon == .( eps_A ) * "%")),
                       as.expression(bquote(epsilon == .( eps_B ) * "%")))
    legend("right", bty='n', cex=text.cex, inset=c(-0.05, 0), 
        text.col=c(red, blue), legend=uncertainties)
    if (has_other) {
        eps_Am <- signif(stds[1]/means[1]*100, 3)
        eps_Bm <- signif(stds[2]/means[2]*100, 3)
        uncertainties <- c(
            as.expression(bquote(epsilon[.(other_name)]==.( eps_Am )*"%")),
            as.expression(bquote(epsilon[.(other_name)]==.( eps_Bm )*"%"))
        )
        legend("topright", bty='n', cex=text.cex, 
            inset=c(-.05, -0.15),#inset=c(-0.075, 0),
            text.col=c(red, blue), legend=uncertainties)
    }
    
    if (show_legend) {
        legend("topleft", bty='n', cex=text.cex,
            inset=c(-0.05, -0.15),#inset=c(-0.075,0),
            pch=c(NA, NA, 20, 20),
            lty=c(2, 1, NA, NA),
            lwd=c(1.5, 1, NA, NA),
            col=c("black", "black", red, blue),
            c("Machine learning",
              other_citation,
              "16 Cyg A",
              "16 Cyg B"))
    }
    
    title(xlab=get_label(name))
}

for (name in names(cygA)) {
    make_plots(plot_cygs, paste0("cyg-", name), 
        filepath=file.path('plots', 'comparison'),
        name=name, cygA=cygA, cygB=cygB,
        mgp.paper=utils.mgp, wide=F, tall=F)
}

cygA <- read.table(file.path('learn_covs-simulations', 'hares', '16CygA.dat'), 
    header=1)
cygB <- read.table(file.path('learn_covs-simulations', 'hares', '16CygB.dat'), 
    header=1)
make_plots(plot_cygs, "cyg-radius", 
    filepath=file.path('plots', 'comparison'),
    name="radius", cygA=cygA, cygB=cygB, 
    mgp.paper=utils.mgp, wide=F, tall=F)

cygA <- read.table(file.path('learn_covs-simulations', 'kages', '16CygA.dat'), 
    header=1)
cygB <- read.table(file.path('learn_covs-simulations', 'kages', '16CygB.dat'), 
    header=1)
make_plots(plot_cygs, "cyg-L", 
    filepath=file.path('plots', 'comparison'),
    name="L", cygA=cygA, cygB=cygB, 
    mgp.paper=utils.mgp, wide=F, tall=F)

cygA <- read.table(file.path('learn_covs-simulations', 'perturb', '16CygA.dat'),
    header=1)
cygB <- read.table(file.path('learn_covs-simulations', 'perturb', '16CygB.dat'),
    header=1)
for (name in c("Y_surf", "X_c", "alpha", "overshoot", "diffusion")) {
    make_plots(plot_cygs, paste0("cyg-", name), 
        filepath=file.path('plots', 'comparison'),
        name=name, cygA=cygA, cygB=cygB,
        mgp.paper=utils.mgp, wide=F, tall=F)
}
#make_plots(plot_cygs, "cyg-Y_surf", 
#    filepath=file.path('plots', 'comparison'),
#    name="Y_surf", cygA=cygA, cygB=cygB, 
#    mgp.paper=utils.mgp, wide=F, tall=F)
#make_plots(plot_cygs, "cyg-X_c", 
#    filepath=file.path('plots', 'comparison'),
#    name="X_c", cygA=cygA, cygB=cygB, 
#    mgp.paper=utils.mgp, wide=F, tall=F)

