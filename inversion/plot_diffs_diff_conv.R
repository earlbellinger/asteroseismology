#### Plot structural differences 
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
#mode.set <- 'CygB'#'BiSON'# 
#error.set <- 'BiSON' 
#perturb <- F 

plot_diffs <- function(r, f1, d.name, col, lty=1, bottext="",
        make.x=T, make.y=T, cz=NULL, 
        ..., text.cex=1, mgp=utils.mgp, font=utils.font, mar=utils.mar, 
        short=F) {
    
    par(mar=mar+c(0.3, -0.2, -0.5, -0.1), lwd=1.66, las=1, cex.axis=text.cex)
    
    xlim <- c(0, 1.00)#5)
    ylim <- c(-0.05, 0.05)
    
    plot(NA, #col=col, lty=lty, #log=log, 
        axes=F, xaxs='i', yaxs='i', 
        #type='l', lwd=3, #lty=2, 
        xlim=xlim, 
        ylim=ylim,
        xlab="",
        ylab="")
    
    in_cz <- F
    bcz <- 0
    for (ii in 1:length(cz)) {
        if (cz[ii]) {
            if (!in_cz) bcz <- r[ii]
            in_cz <- T
        } else {
            if (in_cz && !bcz %in% r[(ii-5):(ii-1)]) 
                rect(ifelse(bcz<0.01, -0.1, bcz), ylim[1], 
                     #ifelse(r[ii-1]>0.99, 1.01, r[ii-1]), 
                     r[ii-1], ylim[2], 
                    col="#DDDDDD", lwd=par()$lwd, 
                border=F)
            in_cz <- F
        }
    }
    
    lines(r[r<0.995], f1[r<0.995], col=col, lty=lty, lwd=3)
    
    abline(h=0, lty=2, lwd=par()$lwd)
    
    text(0.5, -0.049,
        bottext, 
        cex=0.965*text.cex, pos=3, 
        family='Helvetica LT Std Light')
    #abline(v=c(0, 1), lty=2, lwd=1)
    #magaxis(side=c(3,4), tcl=0, labels=c(0,0), 
    #        las=short, mgp=mgp+c(0,0.4,0), 
    #        family=font, cex.axis=text.cex)
    nxticks <- 6
    nyticks <- 4
    nxminor <- 4
    nyminor <- 4
    xticks <- pretty(xlim, n=nxticks)
    yticks <- pretty(ylim, n=nyticks)
    xticks.minor <- pretty(xlim, n=nxticks*nxminor)
    yticks.minor <- pretty(ylim, n=nyticks*nyminor)
    xticks.minor <- xticks.minor[!xticks.minor %in% xticks]
    yticks.minor <- yticks.minor[!yticks.minor %in% yticks]
    par(mgp=mgp+c(0, 0.3, 0))
    axis(1, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, at=xticks,
        labels=as.logical(make.x))
    par(mgp=mgp+c(0, 0.43, 0))
    axis(2, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, at=yticks,
        labels=as.logical(make.y))
    axis(1, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
        at=xticks.minor, labels=F)
    axis(2, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
        at=yticks.minor, labels=F)
    box(lwd=par()$lwd)
    #magaxis(side=c(1), tcl=-0.25, labels=T, 
    #        las=short, mgp=mgp+c(0,0.4,0), majorn=5, #usepar=T, 
    #        family=font, cex.axis=text.cex)
    #magaxis(side=2, tcl=-0.25, labels=T, majorn=4, #usepar=T, 
    #        las=short, mgp=mgp+c(0,0.4,0),
    #        family=font, cex.axis=text.cex)
    #par(mgp=mgp+c(0, 0.4, 0))
    #places <- pretty(f1, n=majorn)
    #labs <- paste0(places * 100, "%")
    #axis(2, at=places, lab=labs, las=TRUE, tick=F, lwd=0, cex.axis=text.cex)
    
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.x) title(xlab=expression(Radius~r/R))
    par(mgp=mgp+c(1.1, 0, 0))
    if (make.y) title(ylab=as.expression(d.name)) 
}

for (ii in c(1:2)) {
    
    if (ii == 1) {
        ref.mod.name <- 'no_diffusion_final'
        target.name <- 'diffusion_final'
    } else {
        ref.mod.name <- 'no_overshoot'
        target.name <- 'overshoot'
    }
    
    #freqs <- get_freqs(target.name=target.name) 
    
    k.pair <- k.pairs[[2]]
    
    ref.mod <- get_model(#freqs=freqs, 
        model.name=ref.mod.name, 
        target.name=target.name, 
        k.pair=k.pair, square.Ks=F)
    
    ref.mod2 <- get_model(#freqs=freqs, 
        model.name=target.name, 
        target.name=ref.mod.name, 
        k.pair=k.pair, square.Ks=F)
    
    r <- seq(min(ref.mod$r), max(ref.mod$r), 0.0001)#ref.mod$r
    #r2 <- seq(min(ref.mod2$r), max(ref.mod2$r), 0.0001)#ref.mod$r
    
    #DF <- read.table(ref.mod$fgong, header=1)
    #G <- 6.67408e-8 #6.6716823E-08
    #sigma <- 5.567063e-5 #5.67051E-5
    #a <- 7.5659122e-15
    #c <- 2.99792458e10
    #nabla_rad <- with(DF, 3 * kappa * L_r * P /
    #    (64 * pi * sigma * G * m *T**4))
    #nabla_rad <- splinefun(DF$x, nabla_rad)(r)
    #nabla_ad <- splinefun(DF$x, DF$nabla_ad)(r)
    #cz <- splinefun(DF$x, nabla_rad - DF$nabla_ad)(r) > 0
    #cz <- nabla_rad > nabla_ad
    #DF <- ref.mod$m.2
    cz2 <- splinefun(ref.mod2$fgong$x, ref.mod2$fgong$conv_stab)(r) < 0
    cz <- splinefun(ref.mod$fgong$x, ref.mod$fgong$conv_stab)(r) < 0
    cz <- cz | cz2
    #cz <- splinefun(log10(DF$T), log10(DF$P))
    
    f1 <- ref.mod$d.f1.spl(r)
    v.name <- k.pair$f1.exp #bquote(.(k.pair$f1.exp))
    if (v.name == bquote(Y)) { d.name <- as.expression(bquote( delta*.(v.name)))
    } else d.name <- as.expression(bquote( delta*.(v.name)/.(v.name) ))
    make_plots(plot_diffs,
        paste0('d_', k.pair$f1, '-', ref.mod$short, '_', ref.mod$target.name),
        filepath=file.path('plots', 'diffs-ov'),
        r=r, f1=f1, d.name=d.name,
        col=ifelse(ii==1, blue, orange),
        #lty=ifelse(ii==1, 1, 3),
        #col=ifelse(k.pair_i==1, col.pal[1], col.pal[3]),
        cex.paper=0.93, wide=F, tall=F, slides=F, make_png=F, 
        make.x=F,
        make.y=ii==1,
        font="Palatino Linotype", 
        bottext=ifelse(ii==1, 'Diffusion', 'Overshoot'),
        cz=cz,
        use.cairo=T)
    
    f1 <- ref.mod$d.f2.spl(r)
    v.name <- k.pair$f2.exp #bquote(.(k.pair$f1.exp))
    d.name <- as.expression(bquote( delta*.(v.name)/.(v.name) ))
    if (v.name == bquote(Y)) { d.name <- as.expression(bquote( delta*.(v.name)))
    } else d.name <- as.expression(bquote( delta*.(v.name)/.(v.name) ))
    make_plots(plot_diffs, 
        paste0('d_', k.pair$f2, '-', ref.mod$short, '_', ref.mod$target.name),
        filepath=file.path('plots', 'diffs-ov'),
        r=r, f1=f1, d.name=d.name,
        #col=ifelse(k.pair_i==1, col.pal[2], ifelse(k.pair_i==2, col.pal[4],
        #    col.pal[5])),
        col=ifelse(ii==1, blue, orange),
        #lty=ifelse(ii==1, 1, 3),
        cex.paper=0.93, wide=F, tall=F, slides=F, make_png=F, 
        font="Palatino Linotype", 
        make.x=1,
        make.y=ii==1,
        bottext=ifelse(ii==1, 'Diffusion', 'Overshoot'),
        cz=cz,
        use.cairo=T)

}
