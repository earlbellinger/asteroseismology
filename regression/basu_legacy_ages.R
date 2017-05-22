#### Plot masses and ages against what we get with machine learning 
#### Also plots the diffusion factor as a function of mass for KAGES stars
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)
library(deming)

red <- adjustcolor(red, alpha.f=0.5)

col.pal <- c('black', 'black', blue, '#900090', red)

data_dir <- file.path('learn-giants', 'covs-simulations', 'legacy')
ml <- do.call(plyr:::rbind.fill, Map(function(covs) {
        name <- sub('.dat', '', basename(covs))
        DF <- read.table(covs, header=1)
        data.frame(KIC = as.integer(name),
                   Age = median(DF$age),
                   dAgeL = median(DF$age) - quantile(DF$age, .16),
                   dAgeH = quantile(DF$age, .84) - median(DF$age),
                   Age_s = sqrt(var(DF$age))
                  )
    }, file.path(data_dir, list.files(data_dir)) ))

basu <- read.table('basu_legacy_ages.txt', header=1)

DF <- merge(basu, ml, by='KIC')

plot_abs_diff <- function(..., text.cex=1, mgp=utils.mgp, font='Palatino') {
    xlims <- range(DF$tau)
    ylims <- with(DF, range(Age-dAgeL, Age+dAgeH))    
    plot(NA, axes=F, ylab="ML Age", xlab="Basu Age", 
        xlim=c(0,15), ylim=c(0,15))
    #abline(a=1, lty=2)
    segments(x0=-5,y0=-5,x1=45,y1=45, lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    arrows(DF$tau, DF$Age - DF$dAgeL, 
           DF$tau, DF$Age + DF$dAgeH, 
        length=0.01, lwd=1.5, angle=90, code=3, col=red)
    with(DF, points(Age ~ tau, pch=20, cex=0.5))
    legend('topleft', bty='n', cex=text.cex, legend=c('LEGACY Data'))
}

make_plots(plot_abs_diff, 'legacy_ml-basu')

