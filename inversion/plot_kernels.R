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
target.name <- 'CygA'
ref.mod <- 'CygAwball'
mode.set <- 'CygA'
perturb <- F
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=perturb) 
m1 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, trim.ks=F) 

k.pair2 <- c2_rho 
m2 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair2, square.Ks=F, trim.ks=F) 

plot_kernels <- function(k, k.pair, ells=1:3, ns=c(5,5,5), ylim=NULL,
        legend.spot=NULL, ..., mar=utils.mar, 
        text.cex=1, mgp=utils.mgp, font="Palatino Linotype", short=F) {
    
    par(mar=mar-c(0.6, 1, 0, 0))
    
    modes <- paste0('l.', ells, '_n.', ns)
    
    if (is.null(ylim)) {
        ylim <- range(sapply(modes, function(mode) k[[mode]]))
        if (ylim[1] >= 0) ylim[1] <- -0.05
    }
    #ylim[2] <- ceil(ylim[2])
    
    xmin <- round(min(do.call(c, 
            sapply(modes, function(mode) k$x[ abs(k[[mode]]) > 0.005 ] ))), 1)
    xlim <- c(xmin, 1)
    
    plot(k$x, k[[modes[1]]], xaxs='i', yaxs='i',
        xlim=xlim,
        ylim=ylim,
        axes=F, type='l', lwd=1.5, 
        xlab="",
        ylab=bquote( 'Kernel'~K^( .(k.pair$f1.exp)*','~.(k.pair$f2.exp)) ))
    par(xpd=NA)
    lines(k$x[k$x>=xmin], k[[modes[2]]][k$x>=xmin], lty=2, col=blue, lwd=1.5)
    lines(k$x[k$x>=xmin], k[[modes[3]]][k$x>=xmin], lty=3, col=orange, lwd=2)
    par(xpd=F)
    abline(h=0, lty=2)
    #pdf.options(encoding='ISOLatin2.enc')
    if (!is.null(legend.spot)) {
        legend(legend.spot, lty=1:3, col=c('black', blue, "#F97100"), 
            cex=text.cex, bty='n', lwd=c(1.5,1.5,2), inset=c(0.02, 0.1),
            legend=#as.expression(do.call(c, Map(function(ell, nn) 
                    #bquote('\u2113' == .(ell)*','~ n == .(nn)),
                #ell=ells, nn=ns))))
              c("\u2113 = 1, n = 5", 
                "\u2113 = 2, n = 5", 
                "\u2113 = 3, n = 5"))
    }
    if (F) {
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.25,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    }
    magaxis(side=1, tcl=-0.25, labels=T,
            las=short, mgp=mgp-c(0.5,0,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=-0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(0,0.25,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    
    par(mgp=mgp-c(0.2, 0, 0))
    title(xlab=expression('Radius'~r/R))
}
#plot_kernels(m1$k1, k.pair);dev.off()


make_plots(plot_kernels, 
    paste0('kernel-', k.pair$f1, '_', k.pair$f2, '-', target.name), 
    filepath=file.path('plots', 'kernels'), k=m1$k1, k.pair=u_Y, 
    ylim=c(-4, 4),
    wide=F, tall=F, use.cairo=T, font="Palatino Linotype", cex.paper=0.75)

make_plots(plot_kernels, 
    paste0('kernel-', k.pair$f2, '_', k.pair$f1, '-', target.name), 
    filepath=file.path('plots', 'kernels'), k=m1$k2, k.pair=Y_u,
    wide=F, tall=F, use.cairo=T, font="Palatino Linotype", cex.paper=0.75,
    ylim=c(-0.05, 1.6),
    legend.spot='topleft')


make_plots(plot_kernels, 
    paste0('kernel-', k.pair2$f1, '_', k.pair2$f2, '-', target.name), 
    filepath=file.path('plots', 'kernels'), k=m2$k1, k.pair=c2_rho, 
    wide=F, tall=F, use.cairo=T, font="Palatino Linotype", cex.paper=0.75,
    ylim=c(-0.2, 6),
    legend.spot='topleft')

make_plots(plot_kernels, 
    paste0('kernel-', k.pair2$f2, '_', k.pair2$f1, '-', target.name), 
    filepath=file.path('plots', 'kernels'), k=m2$k2, k.pair=rho_c2,
    ylim=c(-4.5, 4.5),
    wide=F, tall=F, use.cairo=T, font="Palatino Linotype", cex.paper=0.75)







#make_plots(plot_kernels, 
#    paste0('kernel2-', k.pair$f1, '_', k.pair$f2, '-', m1$short), 
#    filepath=file.path('plots', 'kernels'), k=m1$k1, k.pair=u_Y, 
#    ells=c(2,2,2), ns=c(9,6,3),
#    wide=F, tall=F, use.cairo=T, font="Palatino Linotype", cex.paper=0.7)

#make_plots(plot_kernels, 
#    paste0('kernel2-', k.pair$f2, '_', k.pair$f1, '-', m1$short), 
#    ells=c(2,2,2), ns=c(9,6,3),
#    filepath=file.path('plots', 'kernels'), k=m1$k2, k.pair=Y_u,
#    wide=F, tall=F, use.cairo=T, font="Palatino Linotype", cex.paper=0.7)

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
        lines(model$k1$x, model$k1[[ell.n]], 
            lwd=0.5, col=red)
        #lines(xs, splinefun(model$k1$x, model$k1[[ell.n]])(xs), 
        #    lwd=0.5, col=red)
    }
    #lines(xs, splinefun(proxy.star$k1$x, proxy.star$k1[[ell.n]])(xs),
    #    lwd=0.5, col=1)
    
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
        xlim=c(0, 0.5),
        ylim=c(-0.005, 0.005),
            #range(sapply(model.list, function(model) 
            #splinefun(model$k1$x, model$k1[[ell.n]])(xs) - mid)),
        xlab=expression('Radius'~r/R),
        ylab=bquote( delta*K / max*"("*K*")" ))
    # ^( .(k.pair$f1.exp)*','~.(k.pair$f2.exp))
    abline(h=0, lty=2)
    for (model_i in 1:length(model.list)) {
        if (model_i == 5) next 
        model <- model.list[[model_i]]
        lines(xs, (mid-splinefun(model$k1$x, model$k1[[ell.n]])(xs)) /
                max(mid), 
            lwd=1.5, col=brewer.pal(9, 'Spectral')[[model_i]])
            #col=colorRampPalette(c(red, 1, blue))(9)[model_i])
    }
    #lines(xs, mid-splinefun(proxy.star$k1$x, proxy.star$k1[[ell.n]])(xs),
    #    lwd=2, col=1)
    
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
    wide=F, tall=F, use.cairo=T, font="Times")#, cex.paper=0.8)

}







































### LIBRARIES 
source(file.path('..', 'scripts', 'utils.R'))
source('models.R')
source('frequencies.R')
source('kernels.R')
#parallelStartMulticore(9)

models <- get_model_list()

for (k.pair in list(c2_rho, u_Y, u_Gamma1, rho_Y, Gamma1_rho, c2_Gamma1)) {

#k.pair <- u_Y#rho_Gamma1#rho_c2#rho_c2#rho_c2#c2_rho # 
target.name <- 'diffusion'
ref.mod <- 'no_diffusion'
mode.set <- 'BiSON'
perturb <- F
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=perturb) 
m1 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, trim.ks=F) 

plot_kernels <- function(k, ells=1:3, ns=c(5,5,5), ylim=NULL,
        legend.spot=NULL, make_xlab=T, k.pair=k.pair, swap=F, 
        ..., mar=utils.mar, 
        text.cex=1, mgp=utils.mgp, font="Palatino Linotype", short=F) {
    
    par(mar=mar+c(0.1, -0.2, -0.1, -0.1))
    
    modes <- paste0('l.', ells, '_n.', ns)
    
    if (is.null(ylim)) {
        ylim <- range(sapply(modes, function(mode) k[[mode]]))
        if (ylim[1] >= 0) ylim[1] <- -0.05
    }
    #ylim[2] <- ceil(ylim[2])
    
    xmin <- round(min(do.call(c, 
            sapply(modes, function(mode) k$x[ abs(k[[mode]]) > 0.005 ] ))), 1)
    xlim <- c(xmin, 1)
    xlim <- c(0, 1)
    
    plot(k$x, k[[modes[1]]], xaxs='i', yaxs='i',
        xlim=xlim,
        ylim=ylim,
        axes=F, type='l', lwd=1.5, 
        xlab="",
        ylab="")
    par(xpd=NA)
    lines(k$x[k$x>=xmin], k[[modes[2]]][k$x>=xmin], lty=2, col=blue, lwd=1.5)
    lines(k$x[k$x>=xmin], k[[modes[3]]][k$x>=xmin], lty=3, col=orange, lwd=2)
    par(xpd=F)
    abline(h=0, lty=2)
    #pdf.options(encoding='ISOLatin2.enc')
    if (!is.null(legend.spot)) {
        legend(legend.spot, lty=1:3, col=c('black', blue, "#F97100"), 
            cex=text.cex, bty='n', lwd=c(1.5,1.5,2), inset=c(0.04, 0.04),
            legend=as.expression(do.call(c, Map(function(ell, nn) 
                    bquote('\u2113' == .(ell)*','~ n == .(nn)),
                ell=ells, nn=ns))))
              #c("\u2113 = 1, n = 5", 
              #  "\u2113 = 2, n = 5", 
              #  "\u2113 = 3, n = 5"))
    }
    legend('topleft', bty='n', inset=c(-0.04, -0.04), cex=text.cex, 
        legend=if (swap) k.pair$f2.name else k.pair$f1.name)
    magaxis(side=1, tcl=-0.25, labels=make_xlab,
            las=short, mgp=mgp+c(0,0.35,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=-0.25, labels=1, majorn=4, #usepar=1,
            las=short, mgp=mgp+c(0,0.35,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    
    if (make_xlab) {
        par(mgp=mgp+c(0.6, 0, 0))
        title(xlab=expression('Radius'~r/R))
    }
    par(mgp=mgp+c(0.7, 0, 0))
    title(ylab=if (swap) {
            bquote( K^( .(k.pair$f2.exp) * ',' ~ .(k.pair$f1.exp)) )
         } else {
            bquote( K^( .(k.pair$f1.exp) * ',' ~ .(k.pair$f2.exp)) )
    })
}
#plot_kernels(m1$k1, k.pair);dev.off()

ylim <- if (k.pair$name == c2_rho$name || k.pair$name == Gamma1_rho$name) {
    c(-1, 6)
#} else if (k.pair$name == u_Y$name) {
#    c(-4, 4)
} else if (k.pair$name == rho_Y$name || k.pair$name == u_Gamma1$name || k.pair$name == c2_Gamma1$name || k.pair$name == u_Y$name) {
    c(-4, 4)
} else NULL

make_xlab <- k.pair$name == u_Y$name
legend.spot <- if (k.pair$name == u_Y$name) 'bottomleft' else NULL

make_plots(plot_kernels, 
    paste0('kernel-ell-', k.pair$f1, '_', k.pair$f2, '-', target.name), 
    filepath=file.path('plots', 'kernels2'), k=m1$k1, k.pair=k.pair, swap=F, 
    ylim=ylim, make_xlab=make_xlab, 
    wide=F, tall=F, use.cairo=T, cex.paper=1.2)

make_plots(plot_kernels, 
    paste0('kernel-n-', k.pair$f1, '_', k.pair$f2, '-', target.name), 
    filepath=file.path('plots', 'kernels2'), k=m1$k1, k.pair=k.pair, swap=F, 
    ells=c(2,2,2), ns=c(4,5,6), 
    ylim=ylim, make_xlab=make_xlab, 
    wide=F, tall=F, use.cairo=T, cex.paper=1.2)

ylim <- if (k.pair$name == c2_rho$name) {
    c(-5, 5)
} else if (k.pair$name == u_Y$name || k.pair$name == rho_Y$name) {
    c(0, 1.6)
} else if (k.pair$name == c2_Gamma1$name) {
    c(-3, 9)
} else if (k.pair$name == Gamma1_rho$name) {
    c(-4, 4)
} else if (k.pair$name == u_Gamma1$name) {
    c(-1, 6)
} else NULL

make_plots(plot_kernels, 
    paste0('kernel-ell-', k.pair$f2, '_', k.pair$f1, '-', target.name), 
    filepath=file.path('plots', 'kernels2'), k=m1$k2, k.pair=k.pair, swap=T, 
    wide=F, tall=F, use.cairo=T, cex.paper=1.2, make_xlab=make_xlab, 
    legend.spot=legend.spot, 
    ylim=ylim)#,
    #legend.spot='topleft')

make_plots(plot_kernels, 
    paste0('kernel-n-', k.pair$f2, '_', k.pair$f1, '-', target.name), 
    filepath=file.path('plots', 'kernels2'), k=m1$k2, k.pair=k.pair, swap=T, 
    wide=F, tall=F, use.cairo=T, cex.paper=1.2, make_xlab=make_xlab, 
    ells=c(2,2,2), ns=c(4,5,6), legend.spot=legend.spot, 
    ylim=ylim)#,
    #legend.spot='topleft')


}


