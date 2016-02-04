#### Plot feature importance as trained by the random forest 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)

data.cyg <- read.table(file.path('learn_covs', 'feature-importance-16CygA.dat'), 
    header=1)
data.kages <- read.table(file.path('learn_covs', 
    'feature-importance-3425851.dat'), header=1)
data.hares <- read.table(file.path('learn_covs', 
    'feature-importance-Aardvark.dat'), header=1)

make_boxplot <- function(data, ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
        label=FALSE) {
    par(mgp=c(2, 1, 0))
    sort <- order(sapply(data, median))
    DF <- data[,sort]
    boxplot(DF, horizontal=1, las=1, pch=4, cex=0.5, 
        ylim=if (label != FALSE) c(0, 0.35) else range(DF), 
        #yaxs='i',
        xlab="Feature importance", xaxt='n', tcl=0, cex.axis=text.cex,
        names=as.expression(unlist(Map(function(x) seis.labs[x], 
            names(DF))))
    )
    magaxis(side=c(1,3), family=utils.font, tcl=0.25, labels=c(1,0), 
            mgp=mgp, cex.axis=text.cex)
    if (label != FALSE) {
        legend("bottomright", bty='n', legend=label, cex=text.cex)
    }
}
make_plots(make_boxplot, 'importances-perturb', 
    filepath=file.path('plots', 'importances'),
    mar=c(3, 5, 1, 1), data=data.cyg)
make_plots(make_boxplot, 'importances-kages', 
    filepath=file.path('plots', 'importances'),
    mar=c(3, 5, 1, 1), data=data.kages, label="KAGES")
make_plots(make_boxplot, 'importances-hares', 
    filepath=file.path('plots', 'importances'),
    mar=c(3, 5, 1, 1), data=data.hares, label="Hare and Hound")

