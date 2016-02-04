#### Plot overlapping histograms of 16 Cyg A & B 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)

cygA <- read.table(file.path('learn_covs', 'perturb', '16CygA.dat'), header=1)
cygB <- read.table(file.path('learn_covs', 'perturb', '16CygB.dat'), header=1)

metcalfe <- read.table(file.path('data', 'cyg.dat'), header=1)

measA <- read.table(file.path('data', '16CygA-obs.dat'), header=1)
measB <- read.table(file.path('data', '16CygB-obs.dat'), header=1)

#ylim <- do.call(max, Map(function(name) { 
#        print( max( density(cygA[[name]])$y, density(cygB[[name]])$y ) )
#    }, names(cygA)))

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
        other_A <- density(rnorm(100000, means[1], stds[1]))
        other_B <- density(rnorm(100000, means[2], stds[2]))
        
        xlim <- range(xlim, other_A$x, other_B$x)
        ylim <- range(ylim, other_A$y, other_B$y)
        other_name <- "AMP"
        other_citation <- "Metcalfe et al. 2015"
        
        #other_A$y <- other_A$y / max(ylim)
        #other_B$y <- other_B$y / max(ylim)
    }
    
    #has_meas <- F
    if (name %in% measA$name) {
        has_other <- T
        #has_other <- F
        #has_meas <- T
        means <- c(measA[measA$name == name,]$value, 
                   measB[measB$name == name,]$value)
        stds <- c(measA[measA$name == name,]$uncertainty, 
                  measB[measB$name == name,]$uncertainty)
        other_A <- density(rnorm(100000, means[1], stds[1]))
        other_B <- density(rnorm(100000, means[2], stds[2]))
        xlim <- range(xlim, other_A$x, other_B$x)
        ylim <- range(ylim, other_A$y, other_B$y)
        other_name <- "int"
        other_citation <- "White et al. 2013"
    }
    
    
    #A$y <- A$y / max(ylim)
    #B$y <- B$y / max(ylim)
    
    if (name == "radius") {
        xlim[1] <- xlim[1] - xlim[1]*0.05
    }
    
    if (name == "L") {
        xlim[1] <- xlim[1] - xlim[1]*0.2
    }
    
    par(mar=c(2.5, 1, 1, 1), mgp=mgp-c(0.75,0,0))
    plot(A, axes=F, col='darkred', lwd=2, 
        xlim=xlim, ylim=c(0, ylim[2]*1.01), 
        yaxs='i',
        xlab=get_label(name), ylab="", main="")
    magaxis(side=1, family=font, tcl=0.5, labels=1, 
            las=1, mgp=mgp, cex.axis=text.cex)
    lines(B, col='blue', lwd=2)
    if (has_other) {
        lines(other_A, lty=2, col='darkred')
        lines(other_B, lty=2, col='blue')
    }
    
    eps_A <- signif(sqrt(var(cygA[[name]])) / mean(cygA[[name]]) * 100, 3)
    eps_B <- signif(sqrt(var(cygB[[name]])) / mean(cygB[[name]]) * 100, 3)
    uncertainties <- c(
        as.expression(bquote(epsilon == .( eps_A ) * "%")),
        as.expression(bquote(epsilon == .( eps_B ) * "%"))
    )
    if (has_other) {
        eps_Am <- signif(stds[1]/means[1]*100, 3)
        eps_Bm <- signif(stds[2]/means[2]*100, 3)
        uncertainties <- c(uncertainties, "", 
            as.expression(bquote(epsilon[.(other_name)]==.( eps_Am )*"%")),
            as.expression(bquote(epsilon[.(other_name)]==.( eps_Bm )*"%"))
        )
    }
    legend("left", bty='n', cex=text.cex, inset=c(-0.075, 0),
        text.col=c("darkred", "blue", NA),
        legend=uncertainties,
    )
    
    if (name == "L" || (name == "age" && has_other)) {
        par(xpd=NA)
        legend("topright", bty='n', cex=text.cex,
            inset=c(-0.075,0),
            #lty=c(1,1,2,2), 
            #lwd=c(2,2,1,1), 
            #col=c("darkred", "blue", "darkred", "blue"),
            pch=c(NA, NA, 20, 20),
            lty=c(1, 2, NA, NA),
            lwd=c(2, 1, NA, NA),
            col=c("black", "black", "darkred", "blue"),
            c("Machine learning",
              other_citation,
              "16 Cyg A",
              "16 Cyg B"))
              #"16 Cyg A from Metcalfe et al. 2015",
              #"16 Cyg B from Metcalfe et al. 2015"))
    }
}

for (name in names(cygA)) {
    make_plots(plot_cygs, paste0("cyg-", name), 
        filepath=file.path('plots', 'comparison'),
        name=name, cygA=cygA, cygB=cygB)
}

cygA <- read.table(file.path('learn_covs', 'hares', '16CygA.dat'), header=1)
cygB <- read.table(file.path('learn_covs', 'hares', '16CygB.dat'), header=1)
make_plots(plot_cygs, "cyg-radius", 
    filepath=file.path('plots', 'comparison'),
    name="radius", cygA=cygA, cygB=cygB)

cygA <- read.table(file.path('learn_covs', 'kages', '16CygA.dat'), header=1)
cygB <- read.table(file.path('learn_covs', 'kages', '16CygB.dat'), header=1)
make_plots(plot_cygs, "cyg-L", 
    filepath=file.path('plots', 'comparison'),
    name="L", cygA=cygA, cygB=cygB)



