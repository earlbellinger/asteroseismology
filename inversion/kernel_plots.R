#### Plot kernel functions of a stellar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source('models.R')
source('frequencies.R')
source('kernels.R')
source(file.path('..', 'scripts', 'utils.R'))
#parallelStartMulticore(9)

models <- get_model_list()

k.pair <- u_Y#rho_Gamma1#rho_c2#rho_c2#rho_c2#c2_rho # 
target.name <- 'no_diffusion'
ref.mod <- 'diffusion'
mode.set <- 'BiSON'
perturb <- F
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=perturb) 
m1 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F) 

plot_kernels <- function(k, k.pair, ells=1:3, ns=c(5,5,5),
        legend.spot="topleft", ..., 
        text.cex=1, mgp=utils.mgp, font="Times", short=F) {
    modes <- paste0('l.', ells, '_n.', ns)
    plot(k$x, k[[modes[1]]], xaxs='i', 
        xlim=c(round(min(do.call(c, 
            sapply(modes, function(mode) k$x[ abs(k[[mode]]) > 0.005 ] ))), 1), 
            1),
        ylim=range(sapply(modes, function(mode) k[[mode]])),
        axes=F, type='l', lwd=1.5, 
        xlab=expression('Radius'~r/R),
        ylab=bquote( 'Kernel'~K^( .(k.pair$f1.exp)*','~.(k.pair$f2.exp)) ))
    lines(k$x, k[[modes[2]]], lty=2, col=blue, lwd=1.5)
    lines(k$x, k[[modes[3]]], lty=3, col="#F97100", lwd=2)
    abline(h=0, lty=2)
    legend(legend.spot, lty=1:3, col=c('black', blue, "#F97100"), cex=text.cex, 
        bty='n', lwd=c(1.5,1.5,2),
        legend=do.call(c, Map(function(ell, nn) 
                bquote("\u{2113}" == .(ell)*','~ n == .(nn)),
            ell=ells, nn=ns))
          #c( 
          #  expression("\u{2113}" == .(ells[1])*','~ n == .(ns[1])), 
          #  expression("\u2113" == 1*','~ n == 5), 
          #  expression("\u2113" == 2*','~ n == 5)
        )
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
}
#plot_kernels(m1$k1, k.pair);dev.off()


make_plots(plot_kernels, 
    paste0('kernel-', k.pair$f1, '_', k.pair$f2, '-', m1$short), 
    filepath=file.path('plots', 'kernels'), k=m1$k1, k.pair=u_Y, 
    legend.spot='bottomleft')

make_plots(plot_kernels, 
    paste0('kernel-', k.pair$f2, '_', k.pair$f1, '-', m1$short), 
    filepath=file.path('plots', 'kernels'), k=m1$k2, k.pair=Y_u)

make_plots(plot_kernels, 
    paste0('kernel2-', k.pair$f1, '_', k.pair$f2, '-', m1$short), 
    filepath=file.path('plots', 'kernels'), k=m1$k1, k.pair=u_Y, 
    ells=c(2,2,2), ns=c(9,6,3),
    legend.spot='bottomleft')

make_plots(plot_kernels, 
    paste0('kernel2-', k.pair$f2, '_', k.pair$f1, '-', m1$short), 
    ells=c(2,2,2), ns=c(9,6,3),
    filepath=file.path('plots', 'kernels'), k=m1$k2, k.pair=Y_u)


