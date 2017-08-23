#### Plot kernel functions of a stellar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source(file.path('..', 'scripts', 'utils.R'))
source('models.R')
source('frequencies.R')
source('kernels.R')
#parallelStartMulticore(9)

models <- get_model_list()

k.pair <- u_Y#rho_Gamma1#rho_c2#rho_c2#rho_c2#c2_rho # 
target.name <- 'no_diffusion'
ref.mod <- 'diffusion'
mode.set <- 'BiSON'
perturb <- F
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=perturb) 
m1 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, trim.ks=F) 

plot_kernels <- function(k, k.pair, ells=1:3, ns=c(5,5,5),
        legend.spot=NULL, ..., 
        text.cex=1, mgp=utils.mgp, font="CM Roman", short=F) {
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
    lines(k$x, k[[modes[3]]], lty=3, col=orange, lwd=2)
    abline(h=0, lty=2)
    #pdf.options(encoding='ISOLatin2.enc')
    if (!is.null(legend.spot)) {
        legend(legend.spot, lty=1:3, col=c('black', blue, "#F97100"), 
            cex=text.cex, bty='n', lwd=c(1.5,1.5,2),
            legend=as.expression(do.call(c, Map(function(ell, nn) 
                    bquote("\u2113" 
                        == .(ell)*','~ n == .(nn)),
                ell=ells, nn=ns)))
              #c( 
              #  expression("\u{2113}" == .(ells[1])*','~ n == .(ns[1])), 
              #  expression("\u2113" == 1*','~ n == 5), 
              #  expression("\u2113" == 2*','~ n == 5)
            )
    }
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
    wide=F, tall=F, use.cairo=T, font="Times", cex.paper=0.8)

make_plots(plot_kernels, 
    paste0('kernel-', k.pair$f2, '_', k.pair$f1, '-', m1$short), 
    filepath=file.path('plots', 'kernels'), k=m1$k2, k.pair=Y_u,
    wide=F, tall=F, use.cairo=T, font="Times", cex.paper=0.8,
    legend.spot='topleft')

make_plots(plot_kernels, 
    paste0('kernel2-', k.pair$f1, '_', k.pair$f2, '-', m1$short), 
    filepath=file.path('plots', 'kernels'), k=m1$k1, k.pair=u_Y, 
    ells=c(2,2,2), ns=c(9,6,3),
    wide=F, tall=F, use.cairo=T, font="Times", cex.paper=0.8)

make_plots(plot_kernels, 
    paste0('kernel2-', k.pair$f2, '_', k.pair$f1, '-', m1$short), 
    ells=c(2,2,2), ns=c(9,6,3),
    filepath=file.path('plots', 'kernels'), k=m1$k2, k.pair=Y_u,
    wide=F, tall=F, use.cairo=T, font="Times", cex.paper=0.8)









for (ii in 1:3) {
if (ii == 1) {
    model.list.name <- 'perturbed.model.names'
    target.name <- 'diffusion'
} else if (ii == 2) {
    model.list.name <- 'perturbed.CygA.names'
    target.name <- 'CygAwball'
} else if (ii == 3) {
    model.list.name <- 'perturbed.CygB.names'
    target.name <- 'CygBwball'
}
models <- get_model_list()
k.pair <- u_Y
#model.list.name <- 'perturbed.model.names'#'perturbed.model.names'
#target.name <- 'diffusion'
mode.set <- 'CygA'
freqs <- get_freqs(target.name=target.name, mode.set=mode.set) 
model.names <- get(model.list.name)
model.list <- list()
for (perturbed.model.name in model.names) {
    model <- get_model(freqs=freqs, 
        model.name=perturbed.model.name, 
        target.name=target.name, 
        k.pair=k.pair, square.Ks=F) 
    model.list[[perturbed.model.name]] <- model
}
proxy.star <- get_model(freqs=freqs, model.name=target.name, 
    target.name=target.name,
    k.pair=k.pair, square.Ks=F)
ell.n <- "l.1_n.14"

plot_kernels_all <- function(model.list, ell.n, ..., 
        text.cex=1, mgp=utils.mgp, font="Times", short=F) {
    
    plot(NA, axes=F, xaxs='i', 
        xlim=c(0, 0.45),
        ylim=c(-0.15, 0.6),
            #range(sapply(model.list, function(model) model$k1[[ell.n]])),
        xlab=expression('Radius'~r/R),
        ylab=bquote( 'Kernel'~K^( .(k.pair$f1.exp)*','~.(k.pair$f2.exp)) ))
    abline(h=0, lty=2)
    for (model_i in 1:length(model.list)) {
        model <- model.list[[model_i]]
        lines(xs, splinefun(model$k1$x, model$k1[[ell.n]])(xs), 
            lwd=0.5, col=red)
    }
    lines(xs, splinefun(proxy.star$k1$x, proxy.star$k1[[ell.n]])(xs),
        lwd=0.5, col=1)
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
}
#plot_kernels(m1$k1, k.pair);dev.off()

make_plots(plot_kernels_all, paste0('kernels-all-', model.list.name), 
    model.list=model.list, ell.n=ell.n,
    filepath=file.path('plots', 'kernels'), k.pair=u_Y, 
    wide=F, tall=F, use.cairo=T, font="Times", cex.paper=0.8)


plot_kernels_diffs <- function(model.list, ell.n, ..., 
        text.cex=1, mgp=utils.mgp, font="Times", short=F) {
    
    xs <- seq(0, 1, 0.001)
    mid <- splinefun(model.list[[5]]$k1$x, model.list[[5]]$k1[[ell.n]])(xs)
    plot(NA, axes=F, xaxs='i', 
        xlim=c(0, 1.01),
        ylim=c(-0.1, 0.1),
            #range(sapply(model.list, function(model) 
            #splinefun(model$k1$x, model$k1[[ell.n]])(xs) - mid)),
        xlab=expression('Radius'~r/R),
        ylab=bquote( 'Kernel'~K^( .(k.pair$f1.exp)*','~.(k.pair$f2.exp)) ))
    abline(h=0, lty=2)
    for (model_i in 1:length(model.list)) {
        model <- model.list[[model_i]]
        lines(xs, mid-splinefun(model$k1$x, model$k1[[ell.n]])(xs), 
            lwd=0.5, col=red)#col=colorRampPalette(c(red, 1, blue))(9)[model_i])
    }
    lines(xs, mid-splinefun(proxy.star$k1$x, proxy.star$k1[[ell.n]])(xs),
        lwd=1.5, col=1)
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
}
#plot_kernels_diffs(model.list, k.pair);dev.off()


make_plots(plot_kernels_diffs, paste0('kernel-diffs-', model.list.name), 
    model.list=model.list, ell.n=ell.n,
    filepath=file.path('plots', 'kernels'), k.pair=u_Y, 
    wide=F, tall=F, use.cairo=T, font="Times", cex.paper=0.8)
}




