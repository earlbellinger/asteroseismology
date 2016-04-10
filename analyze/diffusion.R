#### Plot diffusion 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

## Libraries
source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(GGally)
library(scales)

col.pal <- adjustcolor(colorRampPalette(brewer.pal(11, "Spectral"))(21), 
    alpha.f=0.75)

filepath <- file.path('plots', 'diffusion')

## Load data
seis.DF <- data.table(read.table(file.path('..', 
    'forward', 'simulations.dat'), header=1))
setkey(seis.DF, M, Y, Z, alpha, overshoot, diffusion)
keys <- key(seis.DF)

combos <- unique(seis.DF[,keys, with=0])
metals <- unlist(Map(function(i) min(merge(seis.DF, combos[i,])$Fe.H), 
    1:nrow(combos)))
init.FeH <- unlist(Map(function(i) merge(seis.DF, combos[i,])[1,]$Fe.H, 
    1:nrow(combos)))

trunc.metals <- metals
trunc.metals[trunc.metals < -8] <- -8

plot_diffusion <- function(Z=trunc.metals, ..., 
        mar=c(3,3,1,6), text.cex=1, mgp=utils.mgp, thin=F) {
    
    thin.mar <- mar-c(0,0,0,2)
    if (thin) par(mar=thin.mar)
    
    col.pal <- brewer.pal(9, 'Spectral')
    #colorRampPalette(brewer.pal(10, 'Spectral'))(10)
    
    plot(10**-4+combos$diffusion, combos$M,
         axes=FALSE,
         log='x', xaxs='i', yaxs='i', 
         pch=20, cex=0.5, 
         col=col.pal[floor(9*normalize(cbind(trunc.metals, init.FeH)))+1],
         xlim=c(1e-04, 1e+2), 
         ylim=c(0.7, 1.6),
         xlab="", ylab="")
    
    magaxis(side=1:4, family=utils.font, tcl=0.25, labels=c(1,1,0,0),
        mgp=mgp, cex.axis=text.cex, las=1)
    
    Z_max <- max(trunc.metals, init.FeH)
    Z_min <- min(trunc.metals, init.FeH)
    X_range <- 4000#10**diff(par()$usr)[1]
    Zlabs <- signif(quantile(seq(Z_min, Z_max, length=1000), 
                                 c(0, 0.25, 0.5, 0.75, 1)), 1)
    Zlabs[1] <- paste("<", Zlabs[1])
    Zlabs <- paste0(" ", Zlabs)
    color.legend(par()$usr[2]+0.05*X_range, par()$usr[3], 
                 par()$usr[2]+0.10*X_range, par()$usr[4], 
                 Zlabs, 
                 col.pal[1:length(col.pal)], 
                 cex=text.cex, gradient='y', align='rb')
    mtext(expression("Final Surface Metallicity [Fe/H]"["TAMS"]), 
        4, line=ifelse(thin, 2.75, 4.5), cex=text.cex)
    
    points(1, 1, pch=1, cex=1)
    points(1, 1, pch=20, cex=0.1)
    
    title(xlab=get_label('diffusion'))
    title(ylab=get_label('M'))
}
make_plots(plot_diffusion, "Fe_H_M_D", filepath=filepath,
    thin.hack=1, mar=c(3, 3, 1, 6))

plot_diffusion2 <- function(Z=init.FeH, ..., 
        mar=c(3,3,1,6), text.cex=1, mgp=utils.mgp, thin=F) {

    thin.mar <- mar-c(0,0,0,2)
    if (thin) par(mar=thin.mar)

    col.pal <- brewer.pal(9, 'Spectral')
    #colorRampPalette(brewer.pal(10, 'Spectral'))(10)
    
    Z2 <- cbind(trunc.metals, init.FeH)
    
    plot(10**-4+combos$diffusion, combos$M,
         axes=FALSE,
         log='x', xaxs='i', yaxs='i', 
         pch=20, cex=0.5, 
         col=col.pal[
             floor( 9 * ( (Z-min(Z2)) / (max(Z2)-min(Z2)) ) ) + 1
         ],
         xlim=c(1e-04, 1e+2), 
         ylim=c(0.7, 1.6),
         xlab="", ylab="")
    
    magaxis(side=1:4, family=utils.font, tcl=0.25, labels=c(1,1,0,0),
        mgp=mgp, cex.axis=text.cex, las=1)
    
    Z_max <- max(trunc.metals, init.FeH)
    Z_min <- min(trunc.metals, init.FeH)
    X_range <- 4000#10**diff(par()$usr)[1]
    Zlabs <- signif(quantile(seq(Z_min, Z_max, length=1000), 
                                 c(0, 0.25, 0.5, 0.75, 1)), 1)
    Zlabs[1] <- paste("<", Zlabs[1])
    Zlabs <- paste0(" ", Zlabs)
    color.legend(par()$usr[2]+0.05*X_range, par()$usr[3], 
                 par()$usr[2]+0.10*X_range, par()$usr[4], 
                 Zlabs, 
                 col.pal[1:length(col.pal)], 
                 cex=text.cex, gradient='y', align='rb')
    mtext(expression("Initial Surface Metallicity [Fe/H]"["ZAMS"]), 
        4, line=ifelse(thin, 2.75, 4.5), cex=text.cex)
    
    points(1, 1, pch=1, cex=1)
    points(1, 1, pch=20, cex=0.1)
    
    title(xlab=get_label('diffusion'))
    title(ylab=get_label('M'))
}
make_plots(plot_diffusion2, "FeH0_M_D", filepath=filepath,
    thin.hack=1, mar=c(3, 3, 1, 6))

plot_diffusion3 <- function(..., text.cex=1, mgp=utils.mgp) {
    col.pal <- brewer.pal(9, 'Spectral')
    Z <- log10(10**-4 + combos$diffusion)
    
    plot(10**init.FeH, 10**metals, 
         xlab=expression("Initial Metallicity"~"[Fe/H]"["ZAMS"]),
         ylab=expression("Final Metallicity"~"[Fe/H]"["TAMS"]),
         axes=FALSE,
         #xaxs='i', yaxs='i', 
         pch=20, cex=0.5, 
         col=col.pal[
             floor( 9 * ( normalize(Z) ) ) + 1 
         ])#,
         #xlim=c(1e-04, 1e+2), 
         #ylim=c(0.7, 1.6),
         #xlab="", ylab="")
    
    magaxis(side=1:4, family=utils.font, tcl=0.25, labels=c(0,0,0,0),
        mgp=mgp, cex.axis=text.cex, las=1, unlog='xy')
    at <- seq(0, 4, 0.5)
    axis(1, at=at, labels=c(signif(log10(at), 2)), tick=F, 
        cex.axis=text.cex)
    axis(2, at=at, labels=c(signif(log10(at), 2)), tick=F, 
        cex.axis=text.cex, las=1)
    
    Z_max <- max(Z)
    Z_min <- min(Z)
    X_range <- diff(par()$usr)[1]
    Zlabs <- signif(10**quantile(seq(Z_min, Z_max, length=1000), 
                                 c(0, 0.25, 0.5, 0.75, 1)), 1)
    Zlabs <- paste0(" ", Zlabs)
    color.legend(par()$usr[2]+0.05*X_range, par()$usr[3], 
                 par()$usr[2]+0.10*X_range, par()$usr[4], 
                 Zlabs, 
                 col.pal[1:length(col.pal)], 
                 cex=text.cex, gradient='y', align='rb')
    mtext("Diffusion factor D", 4, line=4.5, cex=text.cex)
}
make_plots(plot_diffusion3, "D_FeH0_FeHf", filepath=filepath,
    thin.hack=1, mar=c(3, 3, 1, 6))

