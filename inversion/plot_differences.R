#### Plot structural differences between a diffusion and non-diffusion model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source('../scripts/utils.R') 
source('models.R')
source('frequencies.R')
source('kernels.R')
source('OLA_invert.R')
source('OLA_plots.R')
#num_procs <- max(1, as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
#parallelStartMulticore(num_procs)

### CONSTANTS 
models  = get_model_list() 
target.name <- 'no_diffusion'
ref.mod.name <- 'diffusion'
mode.set <- 'BiSON'
error.set <- 'BiSON'
perturb <- F 

freqs <- get_freqs(target.name=target.name, mode.set=mode.set, 
    error.set=error.set, perturb=perturb) 


plot_profs <- function(r, m1.f1, 
        ..., text.cex=1, mgp=utils.mgp, font=utils.font) {
    plot(r, m1.f1, axes=F, cex=0.1, pch=20, #type='l', #log='xy',
        #xlim=c(0.9995, max(r)),
        xlim=c(0, max(r)), 
        ylim=range(m1.f1), #[r>=0.9995 & r<=max(r)]),
        xlab=expression('Radius'~r/R),
        ylab=expression(c))
    #points(r, m2.f1, cex=0.1, col='darkred')
    points(m.2$x, m.2[[k.pair$f1]], cex=0.1, pch=20, col='darkred')
    #legend("topright", lty=1:2, col=c('black', 'darkred'), cex=text.cex,
    #    legend=c('Profile file', 'FGONG'))
    abline(h=0, lty=2)
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0),
            las=0, mgp=mgp, family=font, cex.axis=text.cex)
}
#make_plots(plot_profs, 'atmosphere', filepath=file.path('plots', 'diffs'))

plot_diffs <- function(r, diffs, d.name, 
        ..., text.cex=1, mgp=utils.mgp, font=utils.font, mar=utils.mar, 
        short=F) {
    #text.cex <- text.cex*1.25
    #par(cex.lab=text.cex)
    #if (short) par(mar=mar+c(0,0.1,0,0))
    par(mar=mar+c(-0.2, 0.6, -0.2, -0.2))
    plot(r, diffs, axes=F, col=red, 
        type='l', lwd=2, #lty=2, 
        xlim=c(0, 1), #c( round(min(r[abs(diffs)>0.00005]), 1), 1 ),
        xlab="",#expression('Radius'~r/R),
        ylab="")#bquote( 'Relative difference in' ~ .(d.name) ))
        #ylab=bquote((.(d.name)[1] - .(d.name)[2]) / .(d.name)[1]))
    #points(r, diffs, pch=20, cex=0.25)
    abline(h=0, lty=2, lwd=1)
    abline(v=c(0, 1), lty=2, lwd=1)
    #points(r[c(1, length(r))],
    #       diffs[c(1, length(r))], pch=20, cex=0.3)
    magaxis(side=c(3,4), tcl=0, labels=c(0,0), 
            las=short, mgp=mgp+c(0,0.2,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=c(1), tcl=-0.25, labels=c(1), 
            las=short, mgp=mgp+c(0,0.2,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=-0.25, labels=1, #usepar=1, 
            las=short, mgp=mgp+c(0,0.3,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    #if (short) par(mgp=mgp+c(1, 0, 0))
    
    par(mgp=mgp+c(0.3, 0, 0))
    title(xlab=expression(Radius~r/R))
    par(mgp=mgp+c(1.8, 0, 0))
    #title(ylab=bquote( 'Difference in' ~ .(d.name) ))
    title(ylab=as.expression(d.name)) #bquote( delta*.(expression(d.name))/.(expression(d.name)) ))
}


for (k.pair_i in 1:length(k.pairs)) {
    k.pair <- k.pairs[[k.pair_i]]
    ref.mod <- get_model(freqs=freqs, model.name=ref.mod.name, 
        target.name=target.name, k.pair=k.pair, square.Ks=F)
    
    r <- seq(min(ref.mod$r), max(ref.mod$r), 0.001)#ref.mod$r
    diffs <- ref.mod$d.f1.spl(r)
    v.name <- k.pair$f1.exp #bquote(.(k.pair$f1.exp))
    if (v.name == bquote(Y)) d.name <- as.expression(bquote( delta*.(v.name)))
    else d.name <- as.expression(bquote( delta*.(v.name)/.(v.name) ))
    make_plots(plot_diffs, 
        paste0('d_', k.pair$f1, '-', ref.mod$short, '_', ref.mod$target.name),
        filepath=file.path('plots', 'diffs'),
        r=r, diffs=diffs, d.name=d.name,
        cex.paper=1.2, wide=F, tall=F, slides=F, make_png=F, 
        font="Palatino Linotype", 
        use.cairo=T)
    
    diffs <- ref.mod$d.f2.spl(r)
    v.name <- k.pair$f2.exp #bquote(.(k.pair$f1.exp))
    d.name <- as.expression(bquote( delta*.(v.name)/.(v.name) ))
    if (v.name == bquote(Y)) d.name <- as.expression(bquote( delta*.(v.name)))
    else d.name <- as.expression(bquote( delta*.(v.name)/.(v.name) ))
    make_plots(plot_diffs, 
        paste0('d_', k.pair$f2, '-', ref.mod$short, '_', ref.mod$target.name),
        filepath=file.path('plots', 'diffs'),
        r=r, diffs=diffs, d.name=d.name,
        cex.paper=1.2, wide=F, tall=F, slides=F, make_png=F, 
        font="Palatino Linotype", 
        use.cairo=T)
    
}
