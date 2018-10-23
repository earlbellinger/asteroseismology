#### Assessing the validity of the scaling relations using the cz radius 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

#options(scipen=100000)
source('../scripts/utils.R') 
library(RColorBrewer)
library(data.table)

#set.seed(0)

DF <- fread('SG_US_step-all.dat') 
#DF2 <- DF[DF$ev_stage==1 & !is.na(DF$cz_radius),] 
DF2 <- DF[!is.na(DF$cz_radius),] 
scaling_R <- with(DF2, nu_max/3090 * (Dnu0/135.1)**(-2) * (Teff/5772)**(1/2))
scaling_M <- with(DF2, (nu_max/3090)**3 * (Dnu0/135.1)**(-4) * (Teff/5772)**(3/2))

#ev_R <- 1 - var(DF2$radius - scaling_R) / var(DF2$radius)
#ev_M <- 1 - var(DF2$M - scaling_M) / var(DF2$M)

bounds <- seq(0, 1, 0.01)

ev_M_mat <- outer(bounds, bounds, Vectorize(function(x, y) {
    idx <- DF2$cz_radius > x & DF2$cz_radius < y
    cat(paste(x, y, sum(idx)), '\n\n\n')
    if (sum(idx) < 10) return(NA) 
    max(-1, 1 - var(DF2$M[idx] - scaling_M[idx]) / var(DF2$M[idx]))
}))

ev_R_mat <- outer(bounds, bounds, Vectorize(function(x, y) {
    idx <- DF2$cz_radius > x & DF2$cz_radius < y
    cat(paste(x, y, sum(idx)), '\n\n\n')
    if (sum(idx) < 10) return(NA) 
    1 - var(DF2$radius[idx] - scaling_R[idx]) / var(DF2$radius[idx])
}))

var_M_mat <- outer(bounds, bounds, Vectorize(function(x, y) {
    idx <- DF2$cz_radius > x & DF2$cz_radius < y
    cat(paste(x, y, sum(idx)), '\n\n\n')
    if (sum(idx) < 10) return(NA) 
    var(DF2$M[idx])
}))

plot_cz_ev <- function(outermat, lab, lvls=c(-1, 0, 0.5, 0.9, 0.95, 0.99, 1),
        ...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.3, -0.3, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    xlim <- c(0, 1)
    ylim <- xlim
    
    filled.contour(bounds, bounds, outermat, 
        levels=lvls,
        xaxs='i', yaxs='i',
        key.axes={
            axis(4, cex.axis=text.cex, tcl=0, line=0, mgp=mgp)
            mtext('Explained Variance', side=4, las=3, line=2, cex=text.cex)
        },
        plot.axes={
            contour(bounds, bounds, outermat, add=TRUE, 
                levels=lvls,
                labcex=0.8*text.cex)
            abline(a=0.01, b=1, lty=2, col=1, lwd=1.5)
            
            abline(v=0.713, lwd=1.5, lty=2)
            abline(h=0.713, lwd=1.5, lty=2)
            points(0.713, 0.713, pch=20, cex=2.5, lwd=1.5, col='white')
            points(0.713, 0.713, pch=20, cex=1, lwd=1.5)
            points(0.713, 0.713, pch=1,  cex=2,  lwd=1.5)
            
            magaxis(side=1:4, family=utils.font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=mgp, cex.axis=text.cex)
            
            par(family="Helvetica", xpd=NA)
            shadowtext(0.5, 0.01, labels=lab, pos=3)
            par(family=font, xpd=F)
            par(xpd=NA)
            rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
            par(xpd=F)
        },
        plot.title={
            title(xlab='Minimum CZ Depth', cex.lab=text.cex, line=2)
            title(ylab='Maximum CZ Depth', cex.lab=text.cex, line=2)
        })
}


make_plots(plot_cz_ev, 'ev_mass', 
    outermat=ev_M_mat,
    lab="Scaling Mass", 
    filepath=file.path('plots', 'cz'), 
    slides=F, #make_png=F, #wide=F, tall=F, 
    paper_pdf_height=4.17309*1.1, 
    cex.paper=0.95) 

make_plots(plot_cz_ev, 'ev_radius', 
    outermat=ev_R_mat,
    lab="Scaling Radius", 
    filepath=file.path('plots', 'cz'), 
    slides=F, #make_png=F, #wide=F, tall=F, 
    paper_pdf_height=4.17309*1.1, 
    cex.paper=0.95) 

make_plots(plot_cz_ev, 'var_mass', 
    outermat=var_M_mat,
    lvls=seq(0, 0.15, 0.01),
    lab="Mass", 
    filepath=file.path('plots', 'cz'), 
    slides=F, #make_png=F, #wide=F, tall=F, 
    paper_pdf_height=4.17309*1.1, 
    cex.paper=0.95) 
