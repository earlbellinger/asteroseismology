#### Calculate the "relative forward problem" in asteroseismology  
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source('../../scripts/utils.R') 
library(RColorBrewer)

c2_rho <- read.table('nus/nus-c2_rho-D_noD.dat', header=1)
Gamma1_rho <- read.table('nus/nus-Gamma1_rho-D_noD.dat', header=1)
c2_Gamma1 <- read.table('nus/nus-c2_Gamma1-D_noD.dat', header=1)
u_Gamma1 <- read.table('nus/nus-u_Gamma1-D_noD.dat', header=1)
u_Y <- read.table('nus/nus-u_Y-D_noD.dat', header=1)
rho_Y <- read.table('nus/nus-rho_Y-D_noD.dat', header=1)

k.pairs <- list(c2_rho, Gamma1_rho, c2_Gamma1, u_Gamma1, u_Y, rho_Y)

col.pal <- colorRampPalette(c(blue, 'black', 'darkred'))(6)
#rev(brewer.pal(length(k.pairs), 'Spectral'))

### PLOT RESULT
plot_diff_diffs <- function(...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font) {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab="",#bquote( 'Kernel Rel. Diff. Minus Actual Rel. Diff.' ), 
        xlim=do.call(range, Map(function(k) range(k$nu.x), k=k.pairs)),
        log='y',
        ylim=c(10**-8, 10**-3))
        #ylim=do.call(range, Map(function(k) range(abs(k$m1.diff)), k=k.pairs)))
    abline(h=0, lty=2)
    for (k.pair.i in 1:length(k.pairs)) {
        k.pair <- k.pairs[k.pair.i][[1]]
        points(k.pair$nu.x, abs(k.pair$m1.diff), 
            col=col.pal[k.pair.i], cex=0.5, pch=k.pair.i)
    }
    legend('bottomright', cex=text.cex, pch=1:length(k.pairs), col=col.pal,
        inset=c(0.125, 0.1),
        legend=as.expression(c(
            bquote(c^2*','~rho),
            bquote(Gamma[1]*','~rho),
            bquote(c^2*','~Gamma[1]),
            bquote(u*','~Gamma[1]),
            bquote(u*','~Y),
            bquote(rho*','~Y)
        )))
    #magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=mgp, 
    #    las=0, cex.axis=text.cex)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=1, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=bquote( 'Kernel Rel. Diff. Minus Actual Rel. Diff.' ))
}

plot_name <- paste0('all_diffs-D_noD')
print(plot_name)
make_plots(plot_diff_diffs, plot_name)

