#### How many samples are needed to reproduce a Gaussian? 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source('../scripts/utils.R')
library(Bolstad)
library(parallel)
library(parallelMap)
library(magicaxis)
library(RColorBrewer)

kl <- function(p, q) {
    sintegral(xs., p * log(p / q))$value
}

xs <- seq(-5, 5, 0.001)
xs. <- xs[2:(length(xs)-1)]
norm <- dnorm(xs.)

samples <- density(rnorm(10**4))
sim <- splinefun(c(-5, samples$x, 5), c(0, samples$y, 0))(xs.)

plot(xs., norm, type='l', lwd=2)
lines(xs., sim, lwd=2, lty=2, col='blue')

col.pal <- brewer.pal(5, 'Spectral')

plot_norm <- function(..., 
        text.cex=1, font=utils.font, mgp=utils.mgp, mar=utils.mar) {
    par(mar=mar+c(0.3, 0, 0, 0))
    plot(xs., norm, axes=F, 
        xaxs='i', 
        xlab="", ylab="", 
        xlim=c(-5, 5), 
        ylim=c(0, 0.5), 
        type='l', lwd=2)
    nns <- c(1,2,3,4,5)
    par(xpd=NA)
    for (ii in 1:length(nns)) {
        nn <- nns[ii]
        samples <- density(rnorm(10**nn))
        sim <- splinefun(c(-5, samples$x, 5), c(0, samples$y, 0))(xs.)
        sim[sim<0]<-0
        lines(xs., sim, lwd=2, lty=3, col=col.pal[ii])
    }
    legend("topleft", pch=20, cex=0.8*text.cex, bty='n', inset=c(0.01, -0.15),
        col=c(1, col.pal), 
        legend=c("Normal", "10 Samples", "100", 
            expression(10^3), expression(10^4), expression(10^5)))
    par(xpd=F)
    abline(h=0, lty=2)
    #magaxis(1:4, tcl=0, labels=F)
    magaxis(side=1, tcl=-0.25, cex.axis=text.cex, labels=T, 
        mgp=mgp+c(0, 0.4, 0), family=font)
    magaxis(side=2, tcl=-0.25, cex.axis=text.cex, labels=T,
        mgp=mgp+c(0, 0.4, 0), las=1, family=font)
    abline(h=0, lty=2)
    #lines(fnns2, fit(fnns2), lwd=1, col=orange)
    par(mgp=mgp+c(0.4, 0, 0))
    title(xlab="x")
    par(mgp=mgp+c(0.6, 0, 0))
    title(ylab=expression(Density~psi(x)))
}

make_plots(plot_norm, "norm", tall=F, wide=F, slides=F, cex.paper=1.18, use.cairo=T,
    font='Palatino Linotype')


parallelStartMulticore(16)
#nns <- 2:6
nns <- seq(1, 6, 0.1)
fnns <- floor(10**nns)
trial <- function(nn) {
    sapply(nns, function(nn) {
        samples <- density(rnorm(floor(10**nn)))
        sim <- splinefun(c(-5, samples$x, 5), c(0, samples$y, 0))(xs.)
        sim[sim<0]<-10**-100#min(sim[sim>0])
        kl(sim, norm)
    })
}
results <- do.call(rbind, parallelMap(trial, 1:1000))
#means <- apply(results, 2, function(x) { mean(abs(x)) })
means <- apply(results, 2, mean)

plot_KL <- function(..., 
        text.cex=1, font=utils.font, mgp=utils.mgp, mar=utils.mar) {
    par(mar=mar+c(0.3, 0, 0, 0))
    fit <- splinefun(fnns, means)
    plot(NA, axes=F, log='x', 
        xaxs='i', #yaxs='i', 
        type='l', lwd=2, col=orange,
        ylim=c(0, 1),#c(0.001, 1), #
        xlim=c(10, 10**5),
        xlab="", 
        ylab="")
    fnns2 <- 10**seq(1, 6, 0.001)
    lines(fnns2, fit(fnns2), lwd=2, col=orange)
    #magaxis(1:4, tcl=0, labels=F)
    magaxis(side=1, tcl=-0.25, cex.axis=text.cex, labels=T, 
        mgp=mgp+c(0, 0.4, 0), family=font)
    magaxis(side=2, tcl=-0.25, cex.axis=text.cex, labels=T,
        mgp=mgp+c(0, 0.4, 0), las=1, family=font)
    abline(h=0, lty=2)
    #lines(fnns2, fit(fnns2), lwd=1, col=orange)
    par(mgp=mgp+c(0.4, 0, 0))
    title(xlab="Number of Samples")
    par(mgp=mgp+c(0.6, 0, 0))
    title(ylab="Relative Entropy")
}

make_plots(plot_KL, "KL", tall=F, wide=F, slides=F, cex.paper=1.18, use.cairo=T,
    font='Palatino Linotype')


quit()
\begin{equation} \label{eq:KL}
    D_{\text{KL}} (P||Q) = \int p(x) \log \frac{p(x)}{q(x)} \;\text{d}x
\end{equation}

