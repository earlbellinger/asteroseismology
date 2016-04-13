#### Plot feature importance as trained by the random forest 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)
library(lattice)
library(grid)

covs_dir <- file.path('learn', 'covs-simulations')
dat.cyg <- read.table(file.path(covs_dir, 'feature-importance-perturb.dat'), 
    header=1)
dat.kages <- read.table(file.path(covs_dir, 'feature-importance-kages.dat'), 
    header=1)
dat.hares <- read.table(file.path(covs_dir, 'feature-importance-basu.dat'), 
    header=1)

make_boxplot <- function(dat, ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
        label=FALSE) {
    par(mgp=c(1.25, 1, 0))
    
    sort <- order(sapply(dat, median))
    DF <- dat[,sort]
    col.pal <- adjustcolor(colorRampPalette(c(blue, red, "black"))(1001)[1+1000*
            sapply(DF, mean)], alpha=0.75)
    
    boxplot(DF, horizontal=1, las=1, pch=20, cex=0.1, #col=col.pal, 
        range=1000,
        ylim=if (label != FALSE) 
                 range(0, round(dat.kages*1.01, 2), round(dat.hares*1.01, 2))
             else range(0, round(DF*1.01, 2)), 
        #border="black", outcol="black", medcol="white", 
        medlwd=1, yaxs='i',
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
    mar=c(3, 4, 1, 2), dat=dat.cyg)
make_plots(make_boxplot, 'importances-kages', 
    filepath=file.path('plots', 'importances'),
    mar=c(3, 5, 1, 1), dat=dat.kages, label="KOI")
make_plots(make_boxplot, 'importances-hares', 
    filepath=file.path('plots', 'importances'),
    mar=c(3, 5, 1, 1), dat=dat.hares, label="Hare and Hound")

plot_cov <- function(dat, ..., text.cex=1) {
    #par(mgp=c(1, 1, 0))
    col.pal <- colorRampPalette(c(blue, 'white', red))(1000)
    cov.mat <- cov(dat)
    sorting <- rev(order(apply(dat, 2, median)))
    cov.mat <- cov.mat[sorting, sorting]
    stand.cov.mat <- apply(cov.mat, 1, function(x) (x - mean(x)) / sqrt(var(x)))
    col.names <- as.expression(seis.labs[colnames(stand.cov.mat)])
    #colnames(stand.cov.mat) <- unicode.labs[colnames(stand.cov.mat)]#col.names
    #rownames(stand.cov.mat) <- unicode.labs[rownames(stand.cov.mat)]#col.names
        #scale(cov.mat, scale=1:ncol(cov.mat))
        #center=1:ncol(cov.mat), 
    #(cov.mat-mean(cov.mat))/sqrt(abs(var(cov.mat)))
    
    #plot(0:1, 0:1, type="n", bty="n", axes=F, xlab="",ylab="")
    #mtext("Standardized covariance", side=4, adj=0.5, line=1, 
    #    cex=text.cex, outer=T)
    a <- levelplot(stand.cov.mat, col.regions=col.pal,
        scales=list(x=list(rot=60), labels=col.names), 
        pretty=T, region=T, text.col='white',
        xlab='', ylab='',
        identifider=col.names)
    #print(a, newpage=F)
    print(a)
    grid.text('Standardized importance covariance', x=0.85, y=.56, rot=90,
    gp = gpar(cex=text.cex))
    
    
    #heatmap(cov.mat, cexRow=text.cex, cexCol=text.cex,
    #    Rowv=NA, Colv=NA, 
    #    col=col.pal,
    #    labCol=col.names, labRow=col.names)
}

make_plots(plot_cov, 'cov-perturb', 
    filepath=file.path('plots', 'importances'),
    mar=c(0, 0, 0, 0), dat=dat.cyg)

