#### Plot Kepler Legacy CDFs for KASC proceedings 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

require(magicaxis)
source(file.path('..', 'scripts', 'utils.R'))
library(RColorBrewer)

#colors <- rev(c('#588C7E', '#F2E394', '#F2AE72', '#D96459', '#8C4646'))
#colors <- c(blue, "black", red)
colors <- c(brewer.pal(11, 'Spectral')[-c(6,7)], "#512854")
#colors <- c("#003F01", "#586F7C", "#4C264F", "#F77F00", "#D62828")
col.pal <- adjustcolor(colorRampPalette(colors)(12))

curr <- rbind(
    read.table(file.path('learn', 'tables-simulations', 
                         'legacy_curr_ascii.dat'), header=1),
    read.table(file.path('learn', 'tables-simulations',
                         'kages_curr_ascii.dat'), header=1))
init <- rbind(
    read.table(file.path('learn', 'tables-simulations', 
                         'legacy_init_ascii.dat'), header=1),
    read.table(file.path('learn', 'tables-simulations',
                         'kages_init_ascii.dat'), header=1))
#curr.names <- expression(tau, X["c"], log~g, L, R, Y["surf"])
curr.names <- c("age", "core X", "log g", "luminosity", "radius", "surface Y")
#init.names <- expression(M, Y[0], Z[0], alpha["MLT"], alpha["ov"], D)
init.names <- c("mass", "initial Y", "initial Z", expression(alpha["MLT"]),
    expression(alpha["ov"]), "diffusion")

mar <- utils.mar + c(0,-1,0,3.5)

exclude <- c(5774694) # fake Sun 
pred.dir <- file.path('learn', 'covs-simulations', 'legacy')
for (pred.fname in list.files(pred.dir)) {
    pred.data <- read.table(file.path(pred.dir, pred.fname), header=1)
    if (any(pred.data$X_c < 0.01)) 
        exclude <- c(exclude, strsplit(pred.fname, '.dat')[[1]])
}
print(exclude)

DF <- merge(curr, init, by='KIC') 
DF <- DF[!DF$KIC %in% exclude,] 
group <- names(DF[,-1])[-grep('^d_', names(DF[,-1]))] 
unc <- DF[paste0('d_', group)] / DF[group] * 100 
sorted <- order(sapply(unc, median)) 
unc <- unc[,sorted] 
print(sapply(unc, mean)) 
print(sapply(unc, fivenum)) 
legend.names <- c(curr.names, init.names)[sorted]
legend.pcts <- as.expression(paste0('(', signif(sapply(unc, median), 1), '%)'))
legend.tot <- sapply(1:length(legend.names), function(ii) 
    bquote(.(legend.names[[ii]])~.(legend.pcts[[ii]])))

plot_cdfs <- function(..., text.cex=1, mgp=utils.mgp, mar=mar, font="Times") {
    
    par(family=font)
    xlim <- c(0.1, 200) 
    #xticks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
    xticks <- c(0.1, 0.5, 2, 10, 50)  
    
    plot(NA, axes=F, xaxs='i', yaxs='i', log='x',
         xlab=expression("Relative uncertainty"),
         ylab=expression("Cumulative count of stars"),
         xlim=xlim,
         ylim=c(0, 1.05))
    
    grid(0, NULL, lty = 6, col="cornsilk2") 
    abline(v=xticks, lty=6, col = "cornsilk2")
    
    #for (ii in 1:length(group)) lines(ecdf(unc[,ii]), col=1, cex=0.5)
    for (ii in 1:length(group)) {
        lines(ecdf(unc[,ii]), col=1, cex=0.5)
        lines(ecdf(unc[,ii]), col=col.pal[ii], cex=0.33)
    }
    
    par(xpd=T)
    legend("right", legend=as.expression(legend.tot),
           inset=c(-0.29, 0),#c(-0.53, 0),
           bty='n', pch=21, lty=F, 
           col=1, 
           pt.bg=col.pal,
           cex=0.8*text.cex)
    par(xpd=F)
    
    magaxis(1:4, labels=F, tcl=-0.25, las=1, 
            mgp=mgp, family=font, cex.axis=text.cex)
    axis(1, xticks, paste0(xticks, '%'), cex.axis=text.cex, tcl=-0.25)
    axis(3, xticks, labels=F, cex.axis=text.cex, tcl=-0.25)
    oh.one <- pretty(c(0, 1))
    axis(2, oh.one, round(seq(0, nrow(unc), length.out=length(oh.one))), 1,
         cex.axis=text.cex, tcl=0, las=1, mgp=mgp+c(0, 0.25, 0))
}

make_plots(plot_cdfs, paste0('cdf2'), mar=mar, 
    paper_pdf_height=4.17309*1.25,
    paper_pdf_width=6.97522*1.5,
    slides=F, wide=F, tall=F)

