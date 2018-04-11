#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

#source('../scripts/utils.R') 
library(mgcv)
library(RColorBrewer)

param.names <- c(bquote(beta), bquote(mu), bquote(Delta[f]))
param.strs <- c("beta", "mu", "deltaf")

# plot_inversion_all <- function(model, inversion, k.pair, should.cull=F,
#         legend.spot='topright', mode.set=mode.set, ylim=NULL, xlim=NULL,
#         normalize=F, sampler=T, plot_nondim=T,
#         ..., text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
#     par(mfrow=c(3,1))
#     plot_inversion(model=model, inversion=inversion, k.pair=k.pair, 
#         should.cull=should.cull, legend.spot=legend.spot, mode.set=mode.set,
#         normalize=normalize, sampler=sampler, plot_nondim=plot_nondim,
#         ylim=ylim, xlim=xlim, text.cex=text.cex, mgp=mgp, mar=mar, font=font)
#     plot_kernels(model=model, inversion=inversion, cross=F, 
#         should.cull=should.cull, legend.spot=legend.spot, normalize=normalize,
#         xlim=xlim, make.inset=F, sampler=sampler,
#         text.cex=text.cex, mgp=mgp, mar=mar, font=font)
#     plot_kernels(model=model, inversion=inversion, cross=T, 
#         xlim=xlim, make.inset=F, sampler=sampler,
#         should.cull=should.cull, legend.spot=legend.spot, normalize=normalize,
#         text.cex=text.cex, mgp=mgp, mar=mar, font=font)
# }

plot_inversion_all <- function(model, inversion, 
        wide=F, tall=F, xlim=c(0, 0.45), 
        kern.interp.xs=NULL,
        inversion_ylim=NULL, 
        avg_kern_ylim=c(-5, 18), 
        make_xlabs=c(T,T,T),
        make_ylabs=c(T,T,T),
        suppress_cross=F,
        suppress_avg=F, 
        cross_kern_ylim=NULL, ...) {
    args. <- c(as.list(environment()), list(...))
    par(mfcol=c(1,3), oma=c(4,3,0,0), mar=c(0,0,0,0))
    do.call(plot_inversion, c(args., list(
            ylim=inversion_ylim,
            make_xlab=make_xlabs[1],
            make_ylab=make_ylabs[1])))
    if (!suppress_avg) do.call(plot_kernels, c(args., list(
            cross=F, 
            ylim=avg_kern_ylim,
            make_xlab=make_xlabs[2],
            make_ylab=make_ylabs[2])))
    if (!suppress_cross) do.call(plot_kernels, c(args., list(cross=T, 
            ylim=cross_kern_ylim,
            make_xlab=make_xlabs[3],
            make_ylab=make_ylabs[3])))
}


make_plots_inversion_all <- function(model, inversion, k.str=NULL, 
        wide=F, tall=F, avg_kern_ylim=c(-5, 18), xlim=c(0, 0.45),  
        make_xlabs=c(T,T,T),
        make_ylabs=c(T,T,T),
        core.bound=NULL, 
        suppress_cross=F, 
        suppress_avg=F, 
        cross_kern_ylim=NULL, inversion_ylim=NULL, ...) {
    args. <- c(as.list(environment()), list(...))
    if (is.null(k.str)) {
        k.str <- paste0(
            '-k_', model$k.pair$short, 
            '_r-', model$name,
            '_t-', model$target.name)
    }
    do.call(make_plots, c(args., list(plot_f=plot_inversion, 
            ylim=inversion_ylim,
            make_xlab=make_xlabs[1],
            make_ylab=make_ylabs[1], 
            filename=paste0('inversion', k.str))))
    if (!suppress_avg) do.call(make_plots, 
        c(args., list(plot_f=plot_kernels, cross=F, 
            ylim=avg_kern_ylim,
            make_xlab=make_xlabs[2],
            make_ylab=make_ylabs[2], 
            filename=paste0('inversion-avg', k.str))))
    if (!suppress_avg) do.call(make_plots, 
        c(args., list(plot_f=plot_kernels, cross=T, 
            ylim=cross_kern_ylim,
            make_xlab=make_xlabs[3],
            make_ylab=make_ylabs[3], 
            filename=paste0('inversion-cross', k.str))))
}


plot_inversion <- function(model, inversion, 
        plot_nondim=F, sampler=T, should.cull=F, 
        xlim=NULL, ylim=NULL, legend.spot=NULL, log='', 
        col.pal=NULL, caption=NULL, caption2=NULL, 
        make_xlab=T, make_ylab=T,
        log.xlim.0=10e-5, core.bound=NULL, 
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    
    result <- inversion$result[sampler,]
    
    if (should.cull) result <- result[with(result, 
        fwhm.left < rs & 
        fwhm.right > rs & 
        err < 0.2 &
        fwhm.right - fwhm.left < 0.2),]
    
    #if (length(xlim)<=1) xlim <- c(
    #    max(0, 0.85*min(result$fwhm.left)),
    #    min(1, 1.15*max(result$fwhm.right)))
    if (is.null(xlim)) xlim <- range(result$fwhm.left, result$fwhm.right)
    if (grepl('x', log, ignore.case=T) && xlim[1]==0) xlim[1] <- log.xlim.0
    if (is.null(ylim)) ylim <- range(if (!is.null(d.f1.true)) {
                d.f1.true[model$r<max(result$rs) & 
                          model$r>min(result$rs) |
                          model$r<max(xlim) &
                          model$r>min(xlim)]
                } else 0, 
      result$df_dr+result$err, result$df_dr-result$err)
      #c(-0.15, 0.07),#
    
    plot(NA, axes=F, log=log, xaxs='i', yaxs='i',
         xlim=xlim, #c(0, 0.4),c(0.05, 0.35),#
         ylim=ylim,
         xlab="",
         ylab="")
    
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=2, lwd=2, col="darkgray")#"#F97100")
    }
    if ('d.f1.nondim' %in% names(model) && plot_nondim && 
            length(model$d.f1.nondim) == length(model$r)) {
        lines(model$r, model$d.f1.nondim, 
            type='l', lty=3, lwd=2, col=blue)
    }
    
    if (!is.null(core.bound)) {
        abline(v=core.bound, lty=2, col='black')
    }
    
    abline(h=0, lty=2, col='black')
    #col.pal <- c(blue)
    #col.pal <- adjustcolor(c('#f97100', blue, 'black', red, "#ffe599"), 
    #    alpha.f=0.8)
    #col.pal <- c('#f97100', blue, 'black', red, "#33a02c")
    if (is.null(col.pal)) {
        col.pal <- brewer.pal(11, "Spectral")[c(1:4,8:11)]
        if (nrow(result) == 6) col.pal <- col.pal[c(1:4, 7:8)]
        col.pal <- adjustcolor(col.pal, alpha.f=0.95)
    }
    #lty <- c(1,2,3,4)
    
    for (res_i in 1:nrow(result)) {
        res <- result[res_i,]
        with(res, arrows(fwhm.left, df_dr, fwhm.right, df_dr, 
            code=3, angle=90, length=0.01, lwd=2.66, lty=1, col=1))
        with(res, arrows(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
            code=3, angle=90, length=0.01, lwd=2.66, lty=1, col=1))
        with(res, arrows(fwhm.left, df_dr, fwhm.right, df_dr, 
            code=3, angle=90, length=0.01, lwd=2, lty=1, col=col.pal))
        with(res, arrows(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
            code=3, angle=90, length=0.01, lwd=2, lty=1, col=col.pal))
        #detach(res)
    }
    
    #with(result, arrows(fwhm.left, df_dr, fwhm.right, df_dr, 
    #    code=3, angle=90, length=0.01, lwd=2.5, lty=1,
    #    col=1))
    #with(result, arrows(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
    #    code=3, angle=90, length=0.01, lwd=2.5, lty=1,
    #    col=1))
    #with(result, arrows(fwhm.left, df_dr, fwhm.right, df_dr, 
    #    code=3, angle=90, length=0.01, lwd=2, lty=1,
    #    col=col.pal))
    #with(result, arrows(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
    #    code=3, angle=90, length=0.01, lwd=2, lty=1,
    #    col=col.pal))
    with(result, points(fwhm.mid, df_dr, col=1, pch=20, cex=0.6))
    #with(result, arrows(fwhm.left, df_dr, fwhm.right, df_dr, 
    #    code=3, angle=90, length=0.02, lwd=2, lty=1,
    #    col=col.pal))
    #with(result, arrows(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
    #    code=3, angle=90, length=0.02, lwd=2, lty=1,
    #    col=col.pal))
    #with(result, points(fwhm.mid, df_dr, col=1, pch=20, cex=1))
    
    par(mgp=mgp-c(0.2, 0, 0))
    if (make_xlab) title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.4, 0, 0))
    f1.exp <- inversion$k.pair$f1.exp
    if (make_ylab) title(ylab=bquote(delta*.(f1.exp)/.(f1.exp)))
    
    if (!is.null(caption)) legend('topleft', bty='n', cex=text.cex,
        legend=caption, inset=c(-0.04, 0.02))
    
    if (!is.null(caption2)) legend('topright', bty='n', cex=text.cex, 
        legend=caption2, inset=c(0.02, 0.02))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(make_xlab,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            #majorn=3,
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=make_ylab, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
}

plot_inversion_multi <- function(inversion.list, 
        plot_nondim=F, sampler=T, should.cull=F, 
        xlim=NULL, ylim=NULL, legend.spot=NULL, log='', 
        col.pal=NULL, caption=NULL, caption2=NULL, 
        make_xlab=T, make_ylab=T,
        log.xlim.0=10e-5, core.bound=NULL, 
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
    
    par(mar=mar+c(0.2,-0.2,-0.4,-0.5), mgp=mgp+c(0, 0.4, 0))
    plot(NA, axes=F, log=log, 
         xaxs='i', yaxs='i',
         xlim=xlim, 
         ylim=ylim,
         xlab="",
         ylab="")
    
    abline(h=0, lty=2, col='black')
    if (is.null(col.pal)) {
        col.pal <- brewer.pal(11, "Spectral")[c(1:4,8:11)]
        if (nrow(result) == 6) col.pal <- col.pal[c(1:4, 7:8)]
        col.pal <- adjustcolor(col.pal, alpha.f=0.95)
    }
    
    min.err <- min(sapply(inversion.list, function(inv.res) 
        min(inv.res$df_dr - inv.res$err)))
    
    result <- inversion.list[[1]]
    for (res_i in 1:nrow(result)) {
        df_dr <- -0.2+(res_i%%2)*0.03 #0.39 - min.err #mean(sapply(inversion.list, function(x) x[res_i,]$df_dr))
        res <- result[res_i,]
        fwhm.left <- res$fwhm.left#[res_i]
        fwhm.right <- res$fwhm.right#[res_i]
        fwhm.mid <- res$fwhm.mid#[res_i]
        err <- res$err#[res_i]
        arrows(fwhm.left, df_dr, fwhm.right, df_dr, 
            code=3, angle=90, length=0.01, lwd=1, lty=1, col=1)
        arrows(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
            code=3, angle=90, length=0.01, lwd=1, lty=1, col=1)
        #with(res, arrows(fwhm.left, df_dr, fwhm.right, df_dr, 
        #    code=3, angle=90, length=0.01, lwd=1, lty=1, col=1))
        #with(res, arrows(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
        #    code=3, angle=90, length=0.01, lwd=1, lty=1, col=1))
        #with(res, arrows(fwhm.left, df_dr, fwhm.right, df_dr, 
        #    code=3, angle=90, length=0.01, lwd=2, lty=1, col=col.pal))
        #with(res, arrows(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
        #    code=3, angle=90, length=0.01, lwd=2, lty=1, col=col.pal))
        #detach(res)
    }
    
    for (result_i in 1:length(inversion.list)) {
        result <- inversion.list[[result_i]]
        xseq <- with(result, seq(min(fwhm.mid), max(fwhm.mid), 0.001))
        if (nrow(result) > 1) {
            with(result, 
                lines(xseq, splinefun(fwhm.mid, df_dr)(xseq),
                #smooth.spline(fwhm.mid, df_dr),
                lwd=3, 
                col=adjustcolor(c(blue, blue, orange, orange)[result_i], 
                    alpha.f=0.3),
                lty=c(1,2,1,2)[result_i]))
        } #else if (nrow(result) > 1) {
          #  with(result, 
          #      lines(fwhm.mid, df_dr,
          #      lwd=1.5, 
          #      col=adjustcolor(c(blue, blue, orange, orange)[result_i], alpha.f=0.25),
          #      lty=c(1,2,1,2)[result_i]))
        
        #}
        #with(result, points(fwhm.mid, df_dr, 
        #    col=c(1,1,red,red)[result_i], #col.pal[result_i], 
        #    pch=c(1, 20, 1, 20)[result_i],
        #    cex=c(0.33, 0.5, 0.33, 0.5)[result_i]))#result_i, cex=0.5))
    }
    
    for (result_i in 1:length(inversion.list)) {
        result <- inversion.list[[result_i]]
        with(result, points(fwhm.mid, df_dr, 
            col=adjustcolor(c(blue,blue,orange,orange)[result_i], alpha.f=1),
            pch=c(1, 20, 1, 20)[result_i],
            lwd=1.5,
            cex=c(0.8, 1.2, 0.8, 1.2)[result_i]))#result_i, cex=0.5))
    }
    
    if (!is.null(caption)) legend('topleft', bty='n', cex=0.9*text.cex,
        legend=caption, inset=c(-0.1, -0.05))
    
    #if (!is.null(caption2)) legend('topleft', bty='n', inset=c(0.02, 0.14),
    if (!is.null(caption2)) legend('topright', #inset=c(0.02, 0.14),
        cex=0.5*text.cex, x.intersp=0.8,
        col=c(blue, blue, orange, orange),
        lty=c(1,2,1,2),
        pch=c(1, 20, 1, 20),
        legend=caption2)
    
    magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=make_xlab)
    magaxis(2, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F)
    axis(2, at=c(-0.2, 0.2), labels=T, tick=F, las=1, cex.axis=text.cex)
    
    par(mgp=mgp+c(0.6, 0, 0))
    if (make_xlab) title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.8, 0, 0))
    f1.exp <- inversion$k.pair$f1.exp
    if (make_ylab) title(ylab=bquote(delta*.(f1.exp)/.(f1.exp)))
    #magaxis(side=c(1,3,4), tcl=0.25, labels=c(make_xlab,0,0),
    #        las=1, mgp=mgp+c(0, 0.1, 0),#-c(0.5,0.15,0), 
    #        family=font, cex.axis=text.cex, majorn=3)
    #magaxis(side=2, tcl=0.25, labels=make_ylab, 
    #        las=1, mgp=mgp+c(0, 0.1, 0), majorn=3, 
    #        family=font, cex.axis=text.cex)
    
    #par(mgp=mgp+c(0.1, 0, 0))
    #if (make_xlab) title(xlab=expression("Radius"~r/R))
    #par(mgp=mgp+c(0.3, 0, 0))
    #f1.exp <- inversion$k.pair$f1.exp
    #if (make_ylab) title(ylab=bquote(delta*.(f1.exp)/.(f1.exp)))
}

plot_inversions <- function(model, inversions, k.pair,
        log='', xlim=NA, legend.spot='topleft', 
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
        font="Times") {
    
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    
    cull <- function(result) {
        result[with(result, 
            fwhm.left < rs & 
            fwhm.right > rs & 
            err < 0.05 & 
            fwhm.right - fwhm.left < 0.15),]
    }
    
    if (length(xlim) <= 1) xlim <- range(sapply(inversions, function(result) 
        with(cull(result), range(fwhm.left, fwhm.right))))
    
    ylim <- range(sapply(inversions, function(result) 
        with(cull(result), range(df_dr+err, df_dr-err, 
            if (!is.null(d.f1.true)) {
                d.f1.true[model$r<max(result$rs) & 
                          model$r>min(result$rs) |
                          model$r<max(xlim) &
                          model$r>min(xlim)]
            } else 0
        )))) 
    
    plot(NA, axes=F, log=log, yaxs='i', 
         xlim=xlim, 
         ylim=ylim, 
         xlab="Radius r/R",
         ylab="")
    abline(h=0, lty=2, col='black')
    
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=2, col='gray', lwd=3)
    }
    
    cols <- c('darkred', 'blue')
    for (ii in 1:length(inversions)) {
        result <- cull(inversions[[ii]])
        with(result, segments(fwhm.left, df_dr, fwhm.right, df_dr, 
            col=adjustcolor(cols[ii],alpha.f=0.5), lwd=1.5))
        with(result, segments(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
            col=adjustcolor(cols[ii],alpha.f=0.5), lwd=1.5))
        with(result, points(fwhm.mid, df_dr, col=cols[ii], pch=20, cex=0.5))
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    #legend(legend.spot, lty=c(2, 1), col=c('darkgray', 'darkred'), lwd=2, 
    #       inset=c(0.01, 0.015),
    #       pch=c(NA, 20), legend=c("Actual", "Inverted"))
    legend(legend.spot, lty=NULL, inset=c(0.01, 0.015), bty='n',
        legend=bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))))
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.exp <- model$f1.exp
    title(ylab=bquote(delta*.(f1.exp)/.(f1.exp)))
}


plot_ref_mods <- function(model.list, inversion, ref.mod=NULL, 
        sampler=T, xlim=NULL, ylim=NULL, legend.spot=NULL, log='', 
        col.pal=NULL, log.xlim.0=10e-5, c.factor=10**15, 
        make_ylab=T, caption=NULL,
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar, 
        font="Times") {
    
    par(mar=mar+c(-0.3,-0.2,0,0))
    
    result <- inversion$result[sampler,]
    
    if (is.null(xlim)) xlim <- range(result$fwhm.left, result$fwhm.right)
    if (grepl('x', log, ignore.case=T) && xlim[1]==0) xlim[1] <- log.xlim.0
    if (is.null(ylim)) ylim <- range(
      result$f+result$f.err, 
      result$f-result$f.err) / c.factor
    
    plot(NA, axes=F, log=log, xaxs='i', yaxs='i',
         xlim=xlim, 
         ylim=ylim,
         xlab="",
         ylab="")
    
    abline(h=0, lty=2, col='black')
    #col.pal <- c(blue)
    #col.pal <- adjustcolor(c('#f97100', blue, 'black', red, "#ffe599"), 
    #    alpha.f=0.8)
    #col.pal <- c('#f97100', blue, 'black', red, "#33a02c")
    if (is.null(col.pal)) {
        col.pal <- brewer.pal(11, "Spectral")[c(1:4,8:11)]
        if (nrow(result) == 6) col.pal <- col.pal[c(1:4, 7:8)]
        col.pal <- adjustcolor(col.pal, alpha.f=0.95)
    }
    lty <- c(1,2,3,4)
    
    for (model in model.list) {
        lines(model$r, model$f1 / c.factor, lwd=1, col='darkgray')
    }
    
    
    if (!is.null(ref.mod)) {
        lines(ref.mod$r, ref.mod$f1 / c.factor, lwd=1, col=1)
    }
    
    
    with(result, arrows(fwhm.left, f / c.factor, fwhm.right, f / c.factor, 
        code=3, angle=90, length=0.02, lwd=2, lty=1,
        col=col.pal))
    with(result, arrows(fwhm.mid, (f+f.err) / c.factor, 
                        fwhm.mid, (f-f.err) / c.factor, 
        code=3, angle=90, length=0.02, lwd=2, lty=1,
        col=col.pal))
    with(result, points(fwhm.mid, f / c.factor, col=col.pal, pch=20, cex=1.25))
    with(result, points(fwhm.mid, f / c.factor, col=1, pch=20, cex=1))
    
    par(mgp=mgp+c(-0.2, 0, 0))
    title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(-0.2, 0, 0))
    if (make_ylab) {
        f1.exp <- inversion$k.pair$f1.exp
        title(ylab=bquote(atop(.(inversion$k.pair$f1.name), 
            .(inversion$k.pair$f1.exp) / 
            (10^15 ~ .(inversion$k.pair$f1.units)))))
    }
    #bquote(delta*.(f1.exp)/.(f1.exp)))
    
    if (!is.null(caption)) legend('bottomleft', bty='n', cex=text.cex,
        legend=caption, inset=c(-0.04, 0.02))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=make_ylab, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
}



plot_ref_rel_diffs <- function(model.list, 
        xlim, ylim, col.pal, 
        rs, inv.lists, avg.kerns.lists, 
        cross.lists, inv.params, kern.interp.xs, 
        make_ylab=T, caption=NULL,
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar, 
        font="Times") {
    
    #text.cex <- text.cex*1.25
    #par(cex.lab=text.cex, mar=mar+c(-0.5, -0.5, 0, 0))
    
    #result <- inversion$result[sampler,]
    
    #if (is.null(xlim)) xlim <- range(result$fwhm.left, result$fwhm.right)
    #if (grepl('x', log, ignore.case=T) && xlim[1]==0) xlim[1] <- log.xlim.0
    #if (is.null(ylim)) ylim <- range(
    #  result$f+result$f.err, 
    #  result$f-result$f.err) / c.factor
    
    plot(NA, axes=F, xaxs='i', yaxs='i',
         xlim=xlim, 
         ylim=ylim,
         xlab="",
         ylab="")
    
    abline(h=0, lty=2, col='black')
    #col.pal <- c(blue)
    #col.pal <- adjustcolor(c('#f97100', blue, 'black', red, "#ffe599"), 
    #    alpha.f=0.8)
    #col.pal <- c('#f97100', blue, 'black', red, "#33a02c")
    if (is.null(col.pal)) {
        col.pal <- brewer.pal(11, "Spectral")[c(1:4,8:11)]
        if (nrow(result) == 6) col.pal <- col.pal[c(1:4, 7:8)]
        col.pal <- adjustcolor(col.pal, alpha.f=0.95)
    }
    lty <- c(1,2,3,4)
    
    
    
    for (model_i in 1:length(model.list)) {
        model <- model.list[[model_i]]
        
        inversion <- lists_to_inversion(model=model, rs=rs, 
            inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
            cross.lists=cross.lists, inv.params=inv.params, 
            kern.interp.xs=kern.interp.xs)
        result <- inversion$result
        
        with(result, arrows(fwhm.left, df_dr, fwhm.right, df_dr, 
            code=3, angle=90, length=0.02, lwd=2, lty=1,
            col=adjustcolor(col.pal[model_i], alpha.f=0.75)))
        with(result, arrows(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
            code=3, angle=90, length=0.02, lwd=2, lty=1,
            col=adjustcolor(col.pal[model_i], alpha.f=0.75)))
    }
    
    for (model_i in 1:length(model.list)) {
        model <- model.list[[model_i]]
        
        inversion <- lists_to_inversion(model=model, rs=rs, 
            inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
            cross.lists=cross.lists, inv.params=inv.params, 
            kern.interp.xs=kern.interp.xs)
        result <- inversion$result
        with(result, points(fwhm.mid, df_dr, col=1, pch=20, cex=1))
        with(result, points(fwhm.mid, df_dr, col=col.pal[model_i], 
            pch=20, cex=0.5))
    }
    
    par(mgp=mgp-c(0.2, 0, 0))
    title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.5, 0, 0))
    f1.exp <- inversion$k.pair$f1.exp
    if (make_ylab) title(ylab=bquote(delta*.(f1.exp)/.(f1.exp)))
    
    if (!is.null(caption)) legend('bottomleft', bty='n', cex=text.cex,
        legend=caption, inset=c(-0.04, 0.02))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0), 
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=make_ylab, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
}




plot_inversion_mean <- function(model, inversion.list,
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    #result <- inversion$result
    f.stds <- sapply(inversion.list, function(inversion) inversion$result$f.err)
    f.weights <- (1/f.stds) / apply(1/f.stds, 1, sum)
    f.means <- sapply(inversion.list, function(inversion) inversion$result$f)
    w.f.means <- sapply(1:nrow(f.weights), 
        function(ii) weighted.mean(f.means[ii,], f.weights[ii,]))
    w.f.stds <- apply(replicate(1000, {
            f.means2 <- sapply(1:ncol(f.means), function(ii) 
                rnorm(length(f.means[,ii]), f.means[,ii], f.stds[,ii]))
            apply(f.means2, 1, sd)
        }), 1, mean)
    fwhm.left <- apply(sapply(inversion.list, function(inversion) 
        with(inversion$result, fwhm.left)), 1, median)
    fwhm.right <- apply(sapply(inversion.list, function(inversion) 
        with(inversion$result, fwhm.right)), 1, median)
    fwhm.mid <- apply(sapply(inversion.list, function(inversion) 
        with(inversion$result, fwhm.mid)), 1, mean)
    
    #m2f <- model.list[[1]]$m2.f1.spl(fwhm.mid)
    
    m.f <- model$f1.spl(fwhm.mid)
    
    plot(NA, axes=F, 
         xlim=c(0.05, 0.35),#range(fwhm.left, fwhm.right),
         ylim=range((m.f - w.f.means-w.f.stds)/m.f, 
                    (m.f - w.f.means+w.f.stds)/m.f),
                    #c(-0.2, .1),
         xlab="Radius r/R",
         ylab="")
    
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=2, col='gray', lwd=3)
    }
    
    segments(fwhm.left,  (m.f - w.f.means) / m.f, 
             fwhm.right, (m.f - w.f.means) / m.f, lwd=2)
    segments(fwhm.mid, (m.f - w.f.means-w.f.stds)/m.f, 
             fwhm.mid, (m.f - w.f.means+w.f.stds)/m.f, lwd=2)
    
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend("bottomright", lty=c(2, NA), pch=c(NA, 3), 
           col=c('gray', 'black'), lwd=3, 
           inset=c(0.01, 0.015),
           legend=c("Actual", "Mean Inversion"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.exp <- model$f1.exp
    title(ylab=bquote('Relative difference' ~ d*.(f1.exp)/.(f1.exp)))
}

plot_inversion_lists_mean <- function(model, inversion.lists, log='',
        xlim=NULL, ylim=NULL, sampler=T, legend.spot="topleft",
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
        font="Times") {
    
    par(mgp=mgp-c(0.5,0,0))
    
    #get_values <- function(row_i, value_name) {
    #    c(sapply(inversion.lists, function(list.) {
    #        sapply(list., function(model.) model.[row_i,][[value_name]])
    #    }))
    #}
    
    #num_rs <- nrow(inv.lists[[1]][[1]])
    #rows <- 1:num_rs[sampler]
    #w.f.means <- sapply(rows, function(row_i) {
    #        fs <- get_values(row_i, 'f')
    #        f.errs <- get_values(row_i, 'f.err')
    #        weighted.mean(fs, 1/f.errs)
    #    })
    
    #fs <- sapply(rows, function(row_i) get_values(row_i, 'f'))
    #f.errs <- sapply(rows, function(row_i) get_values(row_i, 'f.err'))
    #fs <- weighted.mean(fs, 1/f.errs)
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    
    get_q <- function(q_name, lists=inversion.lists) {
        sapply(lists, function(inversion.list) {
            results <- sapply(inversion.list, function(result) result[[q_name]])
            apply(results, 1, mean)
        })
    }
    
    w.f.means <- get_q('f')
    fwhm.left <- get_q('fwhm.left')
    fwhm.mid <- get_q('fwhm.mid')
    fwhm.right <- get_q('fwhm.right')
    
    #w.f.stds <- sapply(inversion.lists, function(inversion.list) {
    #    results <- sapply(inversion.list, function(result) result$f)
    #    apply(results, 1, sd)
    #})
    
    w.f.stds <- apply(w.f.means, 1, sd)[sampler]
    w.f.stds[is.na(w.f.stds)] <- 0
    w.f.means <- apply(w.f.means, 1, mean)[sampler]
    #w.f.stds <- apply(w.f.stds, 1, mean)[sampler]
    #w.f.stds <- apply(w.f.stds, 1, mean)[sampler]
    fwhm.left <- apply(fwhm.left, 1, mean)[sampler]
    fwhm.mid <- apply(fwhm.mid, 1, mean)[sampler]
    fwhm.right <- apply(fwhm.right, 1, mean)[sampler]
    
    m.f <- model$f1.spl(fwhm.mid)
    
    if (is.null(xlim)) xlim <- range(fwhm.left, fwhm.right)
    if (is.null(ylim)) ylim <- range(if (!is.null(d.f1.true)) {
                d.f1.true[model$r<max(xlim) &
                          model$r>min(xlim)]
                } else 0, 
        (m.f - w.f.means + w.f.stds) / m.f, 
        (m.f - w.f.means - w.f.stds) / m.f) + c(-0.01, 0.01)
    
    plot(NA, axes=F, log=log, xaxs='i',
         xlim=xlim, 
         ylim=ylim,
         xlab="",
         ylab="")
    
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=1, lwd=2, col='darkgray')
    }
    if ('d.f1.nondim' %in% names(model) && length(model$d.f1.nondim) > 0) {
        lines(model$r, model$d.f1.nondim, 
            type='l', lty=2, lwd=1, col='darkgray')
    }
    
    abline(h=0, lty=2, col='black')
    col.pal <- adjustcolor(c('#f97100', blue, 'black', red, "#ffe599"), 
        alpha.f=0.8)
    lty <- c(1,2,3,4)
    
    arrows(fwhm.left, (m.f - w.f.means) / m.f, 
          fwhm.right, (m.f - w.f.means) / m.f, 
        code=3, angle=90, length=0.02,
        lwd=2, lty=1, col=col.pal)
    arrows(fwhm.mid, (m.f - w.f.means - w.f.stds) / m.f, 
           fwhm.mid, (m.f - w.f.means + w.f.stds) / m.f, 
        code=3, angle=90, length=0.02,
        lwd=2, lty=1, col=col.pal)
    points(fwhm.mid, (m.f - w.f.means) / m.f, pch=20, cex=1, col=1)
    
    if (!is.null(legend.spot)) {
        legend(legend.spot, lty=NULL, inset=c(0.01, 0.015), 
               bty='n', cex=text.cex,
              legend=as.expression(c(
                  bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))~
                          'kernel pair'),
                  bquote(.(mode.set)~'mode set')
            )))
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.5, 0, 0))
    f1.exp <- model$k.pair$f1.exp
    title(ylab=bquote(delta*.(f1.exp)/.(f1.exp)))
}

plot_inversion_lists_mean_points <- function(model, inversion.lists,
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    get_q <- function(q_name, lists=inversion.lists) {
        sapply(lists, function(inversion.list) {
            apply(sapply(inversion.list, 
                function(result) result[[q_name]]), 1, mean)
        })
    }
    
    w.f.means <- get_q('f')
    fwhm.left <- get_q('fwhm.left')
    fwhm.mid <- get_q('fwhm.mid')
    fwhm.right <- get_q('fwhm.right')
    
    w.f.stds <- sapply(inversion.lists, function(inversion.list) {
        apply(sapply(inversion.list, 
            function(result) result$f), 1, sd)
    })
    
    w.f.means <- apply(w.f.means, 1, mean)
    w.f.stds <- apply(w.f.stds, 1, mean)
    fwhm.left <- apply(fwhm.left, 1, mean)
    fwhm.mid <- apply(fwhm.mid, 1, mean)
    fwhm.right <- apply(fwhm.right, 1, mean)
    
    m.f <- model$f1.spl(fwhm.mid)
    
    
    plot(NA, axes=F, 
         xlim=c(0.05, 0.35),#range(fwhm.left, fwhm.right),
         ylim=c(-0.2, .1),
              #range((m.f - w.f.means-w.f.stds)/m.f, 
              #      (m.f - w.f.means+w.f.stds)/m.f),
              #      #
         xlab="Radius r/R",
         ylab="")
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=2, col='gray', lwd=3)
    }
    
    #segments(fwhm.left,  (m.f - w.f.means) / m.f, 
    #         fwhm.right, (m.f - w.f.means) / m.f, lwd=2)
    #segments(fwhm.mid, (m.f - w.f.means-w.f.stds)/m.f, 
    #         fwhm.mid, (m.f - w.f.means+w.f.stds)/m.f, lwd=2)
    
    for (ii in 1:length(inversion.lists)) {
        inversion.list <- inversion.lists[[ii]]
        for (jj in 1:length(inversion.list)) {
            mdl <- inversion.list[[jj]]
            points(mdl$fwhm.mid, (m.f-mdl$f)/m.f, pch=ii, col=jj, cex=0.2)
        }
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.exp <- model$f1.exp
    title(ylab=bquote(d*.(f1.exp)/.(f1.exp)))
}

plot_results_lists_mean_points <- function(model, results.lists,
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    plot(NA, axes=F, 
         xlim=c(0.05, 0.35),#range(fwhm.left, fwhm.right),
         ylim=range(-0.2, .1,
                    sapply(results.lists, function(results.list) 
                        range(results.list$df_dr))),
              #c(-0.2, .1),
              #range((m.f - w.f.means-w.f.stds)/m.f, 
              #      (m.f - w.f.means+w.f.stds)/m.f),
              #      #
         xlab="Radius r/R",
         ylab="")
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=2, col='gray', lwd=3)
    }
    
    #segments(fwhm.left,  (m.f - w.f.means) / m.f, 
    #         fwhm.right, (m.f - w.f.means) / m.f, lwd=2)
    #segments(fwhm.mid, (m.f - w.f.means-w.f.stds)/m.f, 
    #         fwhm.mid, (m.f - w.f.means+w.f.stds)/m.f, lwd=2)
    
    for (ii in 1:length(results.lists)) {
        results <- results.lists[[ii]]
        points(results$fwhm.mid, results$df_dr, pch=ii, cex=1)
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.exp <- model$f1.exp
    title(ylab=bquote('Relative difference' ~ d*.(f1.exp)/.(f1.exp)))
}

plot_results_lists <- function(model, results.lists,
                               ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    #unlisted <- do.call(rbind, Map(function(inversion.list) 
    #    do.call(rbind, inversion.list), inversion.list=inversion.lists))
    #results.lists <- Map(function(rs) {
    #    unlisted[unlisted$rs==rs,]
    #}, rs=unique(unlisted$rs))
    
    for (ii in 1:length(results.lists)) {
        results <- results.lists[[ii]]
        results <- results[ with(results, fwhm.left < rs & fwhm.right > rs &
                                          r.first_q < rs &  r.third_q > rs), ]
        results.lists[[ii]] <- results
    }
    
    #df.dr <- sapply(results.lists, function(results) 
    #    with(results, weighted.mean(df_dr, 1/err)))
    #stds <- sapply(results.lists, function(results) {
    #    perturbed <- replicate(1000, 
    #        rnorm(nrow(results), results$df_dr, results$err))
    #    mean(apply(perturbed, 2, sd))
    #})
    
    fwhm.mid <- sapply(results.lists, function(results) 
        with(results, weighted.mean(fwhm.mid, 2/(fwhm.right - fwhm.left))))
    fwhm.stds <- sapply(results.lists, function(results) {
        perturbed <- replicate(1000, with(results, 
            rnorm(nrow(results), fwhm.mid, (fwhm.right-fwhm.left)/2)))
        mean(apply(perturbed, 2, sd))
    })
    
    m.f <- model$f1.spl(fwhm.mid)
    
    df.dr <- (m.f-sapply(results.lists, function(results) 
        with(results, weighted.mean(f, 1/f.err))))/m.f
    stds <- sapply(1:length(results.lists), function(ii) {
        results <- results.lists[[ii]]
        perturbed <- (m.f[[ii]] - replicate(1000, 
            rnorm(nrow(results), results$f, results$f.err))) / m.f[[ii]]
        mean(apply(perturbed, 2, sd))
    })
    
    
    #df.dr <- sapply(results.lists, function(results) 
    #    with(results, mean(df_dr)))
    #fwhm.mid <- sapply(results.lists, function(results) 
    #    with(results, mean(fwhm.mid)))
    #stds <- sapply(results.lists, function(results) 
    #    with(results, sd(df_dr)))
    #fwhm.stds <- sapply(results.lists, function(results) 
    #    with(results, sd(fwhm.mid)))
    
    plot(NA, axes=F, 
         xlim=c(0, 0.35),#range(fwhm.left, fwhm.right),
         ylim=c(-0.2, .1),
              #range(-0.2, .1,
                    #sapply(results.lists, function(results.list) 
                    #    range(results.list$df_dr))),
              #c(-0.2, .1),
              #range((m.f - w.f.means-w.f.stds)/m.f, 
              #      (m.f - w.f.means+w.f.stds)/m.f),
              #      #
         xlab="Radius r/R",
         ylab="")
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=2, col='gray', lwd=3)
    }
    
    segments(fwhm.mid-fwhm.stds, df.dr, 
             fwhm.mid+fwhm.stds, df.dr, lwd=2)
    segments(fwhm.mid, df.dr-stds, 
             fwhm.mid, df.dr+stds, lwd=2)
    
    #for (ii in 1:length(results.lists)) {
    #    results <- results.lists[[ii]]
    #    points(results$fwhm.mid, results$df_dr, pch=ii, cex=1)
    #}
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend("bottomright", lty=c(2, NA), pch=c(NA, 3), 
           col=c('gray', 'black'), lwd=3, 
           inset=c(0.01, 0.015),
           legend=c("Actual", "Mean Inversion"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.exp <- model$f1.exp
    title(ylab=bquote('Relative difference' ~ d*.(f1.exp)/.(f1.exp)))
}

plot_results_spline <- function(model, results.lists,
                               ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    #rdf <- do.call(rbind, Map(function(inversion.list) 
    #    do.call(rbind, inversion.list), inversion.list=inversion.lists))
    
    rdf <- do.call(rbind, results.lists)
    gam.mod <- gam(df_dr~s(fwhm.mid), data=rdf)
    
    plot(NA, axes=F, 
         xlim=c(0.05, 0.35),#range(fwhm.left, fwhm.right),
         ylim=range(-0.2, .1,
                    sapply(results.lists, function(results.list) 
                        range(results.list$df_dr))),
              #c(-0.2, .1),
              #range((m.f - w.f.means-w.f.stds)/m.f, 
              #      (m.f - w.f.means+w.f.stds)/m.f),
              #      #
         xlab="Radius r/R",
         ylab="")
    
    lines(gam.mod, lwd=2)
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=2, col='gray', lwd=3)
    }
    
    #segments(fwhm.mid-fwhm.stds, df.dr, 
    #         fwhm.mid+fwhm.stds, df.dr, lwd=2)
    #segments(fwhm.mid, df.dr-stds, 
    #         fwhm.mid, df.dr+stds, lwd=2)
    
    #for (ii in 1:length(results.lists)) {
    #    results <- results.lists[[ii]]
    #    points(results$fwhm.mid, results$df_dr, pch=ii, cex=1)
    #}
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend("bottomright", lty=c(2, NA), pch=c(NA, 3), 
           col=c('gray', 'black'), lwd=3, 
           inset=c(0.01, 0.015),
           legend=c("Actual", "Mean In#version"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.exp <- model$f1.exp
    title(ylab=bquote('Relative difference' ~ d*.(f1.exp)/.(f1.exp)))
}

plot_kernels <- function(model, inversion, cross=F, 
        sampler=T, should.cull=F, normalize=F, 
        inset.plot=F, make.inset=T, cross.inset=NULL, 
        legend.spot=NULL, xlim=c(0, 0.45), ylim=NULL, log='', 
        log.xlim.0=10e-5, kern.interp.xs=NULL,
        make_xlab=T, make_ylab=T, 
        col.pal=NULL, 
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
    
    keep <- if (should.cull) which(with(inversion$result, 
        fwhm.left < rs & 
        fwhm.right > rs & 
        err < 0.05 & 
        fwhm.right - fwhm.left < 0.15)) else T
    
    kernels <- if (cross) inversion$cross_kerns else inversion$avg_kerns
    if (class(kernels) != "matrix") kernels <- matrix(kernels)
    #kernels <- matrix(kernels[,sampler])
    
    if (is.null(kern.interp.xs)) {
        xs <- if (cross) model$k2$x else model$k1$x
    } else {
        xs <- kern.interp.xs
    }
    
    if (should.cull) kernels <- kernels[,keep]
    if (normalize) {
        for (ii in 1:ncol(kernels)) 
            kernels[,ii] <- kernels[,ii] / max(abs(inversion$avg_kerns[,ii]))
    }
    
    #if (is.null(xlim)) xlim <- c(0, 1)
    if (grepl('x', log, ignore.case=T) && xlim[1]==0) xlim[1] <- log.xlim.0
    if (is.null(ylim)) {
        ylim <- if (!cross) range(Map(function(ii) 
                kernels[,ii][xs >= xlim[1] & xs <= xlim[2]], 
            ii=1:ncol(kernels))) else {
                max. <- max(abs(range(kernels)))
                c(-max., max.)
            }
    }
    if (inset.plot & make.inset) {
        if ((!cross || cross && ((abs(ylim[2]) > abs(ylim[1])) 
            || !is.null(legend.spot) && legend.spot=='topright'))
            || (!is.null(cross.inset) && cross.inset!='bottomright'))
            #par(fig=c(0.6, 0.97, 0.4, 0.96), new=T)
            par(fig=c(0.6, 0.99, 0.6, 0.98), new=T)  
        else
            par(fig=c(0.6, 0.97, 0.11, 0.67), new=T)
    } else if (make.inset) {
        par(fig=c(0, 1, 0, 1), mar=mar+c(-0.7, -0.8, -0.6, -0.3))
    }
    
    plot(NA, axes=F, log=log, xaxs='i',
         type='l', lty=2, col='gray', lwd=2, 
         ylim=ylim, ylab="", xlim=xlim, xlab="")
    abline(h=0, lty=2)
    
    if (is.null(col.pal)) {
        col.pal <- brewer.pal(11, "Spectral")[c(1:4,8:11)]
        if (ncol(kernels) == 6) col.pal <- col.pal[c(1:4, 7:8)]
        col.pal <- adjustcolor(col.pal, alpha.f=ifelse(inset.plot, 0.5, 0.9))
    }
    lty <- c(1)
    
    rs <- inversion$result$rs[keep]
    if (!inset.plot) {
        
        for (kern_i in 1:ncol(kernels)) {
            if (length(sampler) < kern_i || sampler[kern_i]) {
                lines(xs, kernels[,kern_i], 
                    lwd=2.6, lty=lty[((kern_i-1) %% length(lty))+ 1], 
                    col=1)
                lines(xs, kernels[,kern_i], 
                    lwd=2, lty=lty[((kern_i-1) %% length(lty))+ 1], 
                    col=col.pal[((kern_i-1) %% length(col.pal))+1])
            }
        }
        
        par(mgp=mgp-c(0.1, 0, 0))
        if (make_xlab) title(xlab=expression("Radius"~r/R))
        par(mgp=mgp+c(0.3 + ifelse(cross, 0.15, 0), 0, 0))
        if (make_ylab) {
            if (cross) {
                title(ylab=bquote('Cross-Term Kernel'~
                    kappa^(.(model$f2.exp)*','~.(model$f1.exp))))
            } else {
                title(ylab=bquote('Averaging Kernel'~
                    kappa^(.(model$f1.exp)*','~.(model$f2.exp))))
            }
        }
        
        #magaxis(side=c(1,3,4), tcl=0.25, labels=c(make_xlab,0,0),
        #    las=1, mgp=mgp-c(0.5,0.15,0), family=font, 
        #    cex.axis=text.cex)
        #magaxis(side=2, tcl=0.25, labels=make_ylab, 
        #    las=1, mgp=mgp+c(1,0,0), family=font, 
        #    cex.axis=text.cex)
        magaxis(side=c(1:4), tcl=0, labels=F, #c(make_xlab,0,0),
            las=1, mgp=mgp, family=font, 
            cex.axis=text.cex)
        magaxis(side=1, tcl=-0.25, labels=make_xlab, 
            las=1, mgp=mgp+c(0,0.1,0), family=font, 
            cex.axis=text.cex)
        magaxis(side=2, tcl=-0.25, labels=make_ylab, 
            las=1, mgp=mgp+c(0,0.4,0), family=font, 
            cex.axis=text.cex)
        
        
        if (make.inset) {
            args. <- c(as.list(environment()), list(...))
            args.$inset.plot <- T
            args.$xlim <- c(xlim[2], max(model$k2$x, 1.05))
            args.$ylim <- NULL
            do.call(plot_kernels, args.)
            #plot_kernels(model, inversion, cross=cross, log=log, 
            #    xlim=c(xlim[2], max(model$k2$x, 1.05)), 
            #    ylim=NULL, sampler=sampler,
            #    should.cull=should.cull, normalize=normalize, inset.plot=T,
            #    cross.inset=cross.inset,
            #    legend.spot=legend.spot, text.cex=text.cex, mgp=mgp, mar=mar,
            #    font=font)
        } 
    
    } else {
    
        rect(0, -10, 1, 10, col='white')
        
        ordering <- rev(order(sapply(1:ncol(kernels), function(kern_i) {
                kern <- kernels[,kern_i]
                kern <- kern[xs > xlim[1]]
                max(abs(kern))
            })))
        
        for (kern_i in ordering) {
            if (length(sampler) < kern_i || sampler[kern_i]) {
                lines(xs, kernels[,kern_i], 
                    lwd=2.66, lty=lty[((kern_i-1) %% length(lty))+ 1], 
                    col=1)
                lines(xs, kernels[,kern_i], 
                    lwd=2, lty=lty[((kern_i-1) %% length(lty))+ 1], 
                    col=col.pal[((kern_i-1) %% length(col.pal))+1])
            }
        }
        
        magaxis(side=1:4, tcl=0, labels=F)
        yaxislabs <- pretty(axis(2, labels=F, tick=F))
        yaxismin <- yaxislabs
        axis(2, at=yaxismin, las=1, 
            labels=F, tick=T, cex.axis=0.8*text.cex, lwd=1, tcl=-0.125/2)
        #axis(4, at=yaxismin, las=1, 
        #    labels=F, tick=T, cex.axis=0.8*text.cex, lwd=1, tcl=0.125/2)
        while (length(yaxislabs) > 3) yaxislabs <- yaxislabs[c(T, F)]
        axis(2, at=yaxislabs, las=1, 
            labels=T, tick=T, cex.axis=0.8*text.cex, lwd=1, tcl=-0.125,
            mgp=mgp+c(0, 0.1, 0))
        #axis(4, at=yaxislabs, las=1, 
        #    labels=F, tick=T, cex.axis=0.8*text.cex, lwd=1, tcl=0.125)
        par(mgp=mgp-c(1, 0.3, 0))
        axis(1, at=c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), las=1, 
            labels=F, tick=T, cex.axis=0.8*text.cex, lwd=1, tcl=-0.125/2)
        axis(1, at=c(0.6, 0.9), las=1, 
            labels=T, tick=T, cex.axis=0.8*text.cex, lwd=1, tcl=-0.125)
        #axis(3, at=c(0.5, 1), las=1, 
        #    labels=F, tick=T, cex.axis=0.8*text.cex, lwd=1, tcl=0.125)
    }
}


plot_kernel_lists <- function(model, kernel.lists, cross=F, log='', 
        xs=NULL, xlim=c(0, 0.45), ylim=NULL, normalize=F, inset.plot=F, 
        make.inset=T, legend.spot='topright', 
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
    #text.cex <- text.cex*1.25
    #par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    #par(mgp=mgp-c(0.5,0,0))
    
    if (is.null(xs)) xs <- if (cross) model$k2$x else model$k1$x
    
    #xs <- seq(0, 1, 0.001)
    kernels <- do.call(cbind, Map(function(radius_i) {
        apply(do.call(cbind, Map(function(ii) {
            trial. <- kernel.lists[[ii]]
            do.call(cbind, Map(function(jj) {
                ref.mod. <- trial.[[jj]]
                #splinefun(ref.mod.$x, ref.mod.[,radius_i])(xs)
                ref.mod.[,radius_i]
            }, jj=1:length(trial.)))
        }, ii=1:length(kernel.lists))), 1, mean)
    }, radius_i=1:ncol(kernel.lists[[1]][[1]])))[,sampler]
    
    if (normalize) {
        for (ii in 1:ncol(kernels)) {
            kernels[,ii] <- kernels[,ii] / max(abs(inversion$avg_kerns[,ii]))
        }
    }
    
    if (is.null(ylim)) {
        ylim <- if (!cross) range(Map(function(ii) 
                kernels[,ii][xs >= xlim[1] & xs <= xlim[2]], 
            ii=1:ncol(kernels))) else range(kernels)
    }
    
    if (inset.plot & make.inset) {
        if (!cross || cross && (abs(ylim[2]) > abs(ylim[1])))
            par(fig=c(0.6, 0.97, 0.3, 0.96), new=T) 
        else
            par(fig=c(0.6, 0.97, 0.11, 0.77), new=T)
    } else if (make.inset) {
        par(fig=c(0, 1, 0, 1))
    }
    
    plot(NA, log=log, xaxs=if (inset.plot) 'r' else 'i',
         type='l', lty=2, col='gray', lwd=3, 
         ylim=ylim,
         xlim=xlim,
         axes=F, 
         xlab="",
         ylab="")
    
    abline(h=0, lty=2)
    col.pal <- adjustcolor(c('#f97100', blue, 'black', red, "#ffe599"), 
        alpha.f=1)#0.8)
    lty <- c(1,2,3,4)
    for (kern_i in 1:ncol(kernels)) lines(xs, kernels[,kern_i], 
        lwd=2, lty=lty[((kern_i-1) %% length(lty))+ 1], 
        col=col.pal[((kern_i-1) %% length(col.pal))+1])
    
    if (!is.null(legend.spot) && !cross && F) {
        params <- inversion$params
        cross.term <- bquote(beta==.(params$cross.term))
        error.sup <- bquote(mu==.(params$error.sup))
        legend.txt <- c(cross.term, error.sup)
        if (!is.null(params$width)) 
            legend.txt <- c(legend.txt, bquote(Delta==.(params$width)))
        legend(legend.spot, lty=NULL, inset=c(0.01, 0.015), 
               bty='n', cex=text.cex, legend=as.expression(legend.txt))
    }
    
    if (!inset.plot) {
        
        magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
        magaxis(side=2, tcl=0.25, labels=1, 
                las=1, mgp=mgp+c(1,0,0),
                family=font, cex.axis=text.cex)
        
        title(xlab=expression("Radius"~r/R))
        par(mgp=mgp+c(0.5, 0, 0))
        if (cross) {
            title(ylab=bquote('Cross-Term Kernel'~
                                  kappa^(.(model$f2.exp)*','~.(model$f1.exp))))
        } else {
            title(ylab=bquote('Averaging Kernel'~
                                  kappa^(.(model$f1.exp)*','~.(model$f2.exp))))
        }
        
        if (make.inset) plot_kernel_lists(model, kernel.lists, cross=cross, 
                log=log, xlim=c(xlim[2], max(model$k2$x)), ylim=NULL, 
                normalize=normalize, inset.plot=T,
                legend.spot=NULL, text.cex=text.cex, mgp=mgp, mar=mar,
                font=font)
    
    } else {
        magaxis(side=c(1,3), tcl=0.125, labels=c(1,0),
            las=1, mgp=mgp-c(1,0.3,0), majorn=2, minorn=2,
            family=font, cex.axis=text.cex)
        magaxis(side=c(2,4), tcl=0.125, labels=c(1,0), 
            las=1, mgp=mgp+c(1,0,0), majorn=3, minorn=2,
            family=font, cex.axis=text.cex)
        #axis(side=1, at=xlim, cex=text.cex, tick=F)
        ##y.at <- pretty(ylim)
        ##y.at <- c(y.at[1], y.at[length(y.at)])
        #y.at <- signif(ylim, 2)
        #if (0 > ylim[1] && 0 < ylim[2]) y.at <- c(y.at[1], 0, y.at[2])
        #axis(side=2, at=y.at, cex=text.cex, tick=F)
    }
}


make_fiducial_kernels_plots <- function(model, inversion) {
    for (cross in c(F,T)) {
        for (param_i in 1:3) {
            plot.name <- paste0("SOLA_inversion-",
                                model$short, "_", 
                                target.name, "_",
                                mode.set, "-",
                                "fid_kern-",
                                if (cross) "cross-" else "avg-",
                                param.strs[param_i]
            )
            make_plots(plot_fiducial_kernel, 
                       filename=plot.name,
                       model=model, 
                       inversion=inversion,
                       param_i=param_i,
                       cross=cross,
                       make_png=F, short=T, thin=T)
        }
    }
}

plot_fiducial_kernel <- function(model, inversion, param_i, cross=F,
                                 ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                                 font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, -0.5, 0, 0))
    
    kern_type <- if (cross) "cross_kerns" else "avg_kerns"
    
    xlim <- if (cross) c(0, 0.9) else c(0, 0.5)
    plot(NA, 
         type='l', lty=2, col='gray', lwd=3, 
         ylim=1.25*range(inversion[[kern_type]][,1][model$k1$x<max(xlim)]), 
         xlim=xlim, 
         axes=F, 
         xlab="Radius r/R", 
         ylab="")
    
    abline(h=0, lty=2, col='gray')
    lines(model$k1$x, inversion[[kern_type]][,1], lty=1, lwd=2, col=1)
    
    multipliers <- if (param_i == 1) {
        c(0.1, 10)
    } else if (param_i == 2) {
        c(0.25, 1.75)
    } else c(0.75, 1.25)
    for (mult_i in 1:length(multipliers)) {
        multiplier <- multipliers[mult_i]
        params <- inversion$params 
        params[param_i] <- params[param_i] * multiplier 
        inversion2 <- invert.OLA(model=model, 
                                 rs=result$rs, 
                                 cross.term=params[[1]], 
                                 error.sup=params[[2]], 
                                 width=params[[3]]) 
        lines(model$k1$x, inversion2[[kern_type]][,1], lty=mult_i+1, 
              col=c(blue, red)[mult_i], lwd=2)
    }
    
    params <- inversion$params
    legend(if (cross) "bottomright" else "topright", 
           col=c(blue, "black", red), lty=c(2,1,3), lwd=2,
           cex=text.cex*0.8,
           inset=c(0.025, 0.05), 
           legend=as.expression(c(
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]]*multipliers[1], 2))),
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]], 2))),
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]]*multipliers[2], 2))))))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp,#+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    par(mgp=mgp+c(1.5, 0, 0))
    if (cross) {
        title(ylab=bquote('Cross-Term Kernel'~
                              kappa^(.(model$f2.exp)*','~
                                 .(model$f1.exp))),
              line=2)
    } else {
        title(ylab=bquote('Averaging Kernel'~
                              kappa^(.(model$f1.exp)*','~
                                 .(model$f2.exp))),
              line=2)
    }
}

make_sensitivities_plots <- function(model, inversion) {
    for (param_i in 1:3) {
        plot.name <- paste0("SOLA_inversion-Gauss-",
                            model$short, "_", 
                            target.name, "_",
                            mode.set, "-",
                            "sensitivity-",
                            param.strs[param_i]
        )
        make_plots(plot_sensitivity, 
                   filename=plot.name,
                   model=model, 
                   inversion=inversion,
                   param_i=param_i,
                   make_png=F, short=T, thin=T, wide=F, tall=F)
    }
}

plot_sensitivity <- function(model, inversion, param_i, 
                             ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                             font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, -0.5, 0, 0))
    
    result <- inversion$result
    
    plot(NA, axes=F, 
         xlim=c(0, 0.4),
         ylim=range(result$df_dr+result$err, result$df_dr-result$err)*1.25,
         xlab="Radius r/R",
         ylab="")
    
    abline(h=0, lty=2, col='#00000099')
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=2, col='gray', lwd=3)
    }
    
    for (row_i in 1:nrow(result)) {
        row <- result[row_i,]
        with(row, segments(fwhm.left, df_dr, fwhm.right, df_dr, 
                           col=adjustcolor('black',alpha.f=0.5)))
        with(row, segments(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err,
                           col=adjustcolor('black',alpha.f=0.5)))
        with(row, points(fwhm.mid, df_dr, col='black', pch=20, cex=0.5))
    }
    
    multipliers <- if (param_i == 1) {
        c(0.75, 1.25)
        } else if (param_i == 2) {
            c(0.75, 1.25)
        } else c(0.75, 1.25)
    for (mult_i in 1:length(multipliers)) {
        multiplier <- multipliers[mult_i]
        params <- inversion$params 
        params[param_i] <- params[param_i] * multiplier 
        inversion2 <- invert.OLA(model=model, 
                                 rs=result$rs, 
                                 cross.term=params[[1]], 
                                 error.sup=params[[2]], 
                                 width=params[[3]],
                                 targ.kern.type=inversion$targ.kern.type) 
        color <- c(blue, red)[mult_i]
        result <- inversion2$result
        for (row_i in 1:nrow(result)) {
            row <- result[row_i,]
            with(row, segments(fwhm.left, df_dr, fwhm.right, df_dr, 
                               col=adjustcolor(color,alpha.f=0.5)))
            with(row, segments(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err,
                               col=adjustcolor(color,alpha.f=0.5)))
            with(row, points(fwhm.mid, df_dr, col=color, pch=20, cex=0.5))
        }
        #lines(model$k1$x, inversion2[[kern_type]][,1], lty=mult_i+1, 
        #      col=c(blue, red)[mult_i], lwd=2)
    }
    
    params <- inversion$params
    legend("bottomright", 
           col=c(blue, "black", red), #lty=c(2,1,3), lwd=2,
           pch=3, 
           cex=text.cex*0.8,
           inset=c(0.025, 0.05), 
           legend=as.expression(c(
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]]*multipliers[1], 2))),
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]], 2))),
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]]*multipliers[2], 2))))))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    par(mgp=mgp+c(0.5, 0, 0))
    f1.exp <- model$f1.exp
    title(ylab=bquote('Relative difference' ~ d*.(f1.exp)/.(f1.exp)))
}


make_MOLA_sensitivity_plots <- function(model, inversion) {
    for (mult_i in 1:3) {
        for (fix_mu in 0:1) {
            plot.name <- paste0("MOLA_inversion-",
                                model$short, "_", 
                                target.name, "_",
                                mode.set, "-",
                                "sensitivity-",
                                fix_mu, "_", mult_i
            )
            make_plots(plot_MOLA_sensitivity,  
                       filename=plot.name, 
                       model=model, 
                       inversion=inversion, 
                       mult_i=mult_i, fix_mu=fix_mu,
                       make_png=F, short=T, thin=T, tall=F, wide=F)
        }
    }
}

plot_MOLA_sensitivity <- function(model, inversion, fix_mu, mult_i, 
                             ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                             font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, -0.5, 0, 0))
    
    result <- inversion$result
    params <- inversion$params
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    
    plot(NA, axes=F, 
         xlim=c(max(0, min(result$fwhm.left) - 0.01),
                min(1, max(result$fwhm.right) + 0.01)),
         ylim=range(if (!is.null(d.f1.true)) {
                        d.f1.true[model$r<max(result$rs) & 
                                  model$r>min(result$rs)]
                        } else 0, 
             result$df_dr+result$err, result$df_dr-result$err),
         #xlim=c(0, 0.4),
         #ylim=range(result$df_dr+result$err, result$df_dr-result$err)*1.25,
         xlab="Radius r/R",
         ylab="")
    
    abline(h=0, lty=3, col='#00000099')
    
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=2, col='gray', lwd=3)
    }
    
    mu.mults <- c(0.25, 1, 2)
    beta.mults <- c(0, 1, 2)
    
    if (fix_mu) {
        
        mu <- params[[2]] * mu.mults[mult_i]
        
        for (mult_j in 1:3) {
            
            beta <- params[[1]] * beta.mults[mult_j]
            inversion2 <- invert.OLA(model=model, rs=result$rs, 
                                     cross.term=beta, 
                                     error.sup=mu) 
            color <- c(blue, 'black', red)[mult_j]
            result <- inversion2$result
            for (row_i in 1:nrow(result)) {
                row <- result[row_i,]
                with(row, segments(fwhm.left, df_dr, fwhm.right, df_dr, 
                                   col=adjustcolor(color,alpha.f=0.5)))
                with(row, segments(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err,
                                   col=adjustcolor(color,alpha.f=0.5)))
                with(row, points(fwhm.mid, df_dr, col=color, pch=20, cex=0.5))
            }
            
        }
        
        param_i <- 1
        legend("bottomright", 
           col=c(NA, blue, "black", red), #lty=c(2,1,3), lwd=2,
           pch=3, 
           cex=text.cex*0.8,
           inset=c(0.025, 0.05), 
           legend=as.expression(c(
               bquote(.(param.names[[2]])*' = '*.(signif(mu, 2))),
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]]*beta.mults[1], 2))),
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]], 2))),
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]]*beta.mults[3], 2))))))
        
    } else {
        
        beta <- params[[1]] * beta.mults[mult_i]
        
        for (mult_j in 1:3) {
            
            mu <- params[[2]] * mu.mults[mult_j]
            inversion2 <- invert.OLA(model=model, rs=result$rs, 
                                     cross.term=beta, 
                                     error.sup=mu) 
            color <- c(blue, 'black', red)[mult_j]
            result <- inversion2$result
            for (row_i in 1:nrow(result)) {
                row <- result[row_i,]
                with(row, segments(fwhm.left, df_dr, fwhm.right, df_dr, 
                                   col=adjustcolor(color,alpha.f=0.5)))
                with(row, segments(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err,
                                   col=adjustcolor(color,alpha.f=0.5)))
                with(row, points(fwhm.mid, df_dr, col=color, pch=20, cex=0.5))
            }
            
        }
        
        param_i <- 2
        legend("bottomright", 
           col=c(NA, blue, "black", red), #lty=c(2,1,3), lwd=2,
           pch=3, 
           cex=text.cex*0.8,
           inset=c(0.025, 0.05), 
           legend=as.expression(c(
               bquote(.(param.names[[1]])*' = '*.(signif(beta, 2))),
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]]*mu.mults[1], 2))),
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]], 2))),
               bquote(.(param.names[[param_i]])*' = '*
                          .(signif(params[[param_i]]*mu.mults[3], 2))))))
        
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    par(mgp=mgp+c(0.5, 0, 0))
    f1.exp <- model$f1.exp
    title(ylab=bquote('Relative difference' ~ d*.(f1.exp)/.(f1.exp)))
}

plot_inv_coefs <- function(rs, inversion,
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    plot(NA, 
         ylim=range(inversion$inv.coefs), 
         xlim=range(rs), 
         axes=F, 
         xlab="Radius r/R", 
         ylab="")
    
    abline(h=0, lty=2, col='gray')
    
    num_rows <- nrow(inversion$inv.coefs)
    col.pal <- colorRampPalette(c(blue, 'black', red))(num_rows)
    for (ii in 1:num_rows)
        lines(rs, inversion$inv.coefs[ii,], lty=ii, lwd=1.5,
              col=col.pal[ii])
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp,#+c(1,0,0),
            family=font, cex.axis=text.cex)
    
}

plot_convolution <- function(model, inversion,
                             ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                             font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 1.5, 0, 0))
    
    d.f1.spl <- splinefun(model$k1$x[-1], model$d.f1.true)
    convs <- sapply(1:ncol(inversion$avg_kerns), function(avg_kern_i) { 
        ker.spl <- splinefun(model$k1$x, inversion$avg_kerns[,avg_kern_i])
        
        conv <- sapply(1:length(model$k1$x), 
               function(x_i) {
            x <- model$k1$x[x_i]
            
            #ker <- ker.spl(x-model$k1$x)
            #ker[x-model$k1$x < 0] <- 0
            #d.f1 <- d.f1.spl(model$k1$x)
            d.f1 <- d.f1.spl(x-model$k1$x)
            d.f1[x-model$k1$x < 0] <- 0
            ker <- ker.spl(model$k1$x)
            
            sintegral(model$k1$x, d.f1 * ker )$value
        })
        
        conv
    })
    
    plot(NA, xaxs='i',
         ylim=range(convs), 
         xlim=c(0, 1), 
         axes=F, 
         xlab="Radius r/R", 
         ylab="")
    
    abline(h=0, lty=2, col='gray')
    
    lines(model$r, model$d.f1.true, type='l', lty=2, col='gray', lwd=3)
    
    col.pal <- colorRampPalette(c(blue, 'black', red))(ncol(convs))
    for (ii in 1:ncol(convs))
        lines(model$k1$x, convs[,ii], lty=ii, lwd=1.5,
              col=col.pal[ii])
    
    legend("bottomright", lty=1:ncol(convs), lwd=1.5, 
           col=col.pal[1:ncol(convs)],
           inset=c(0.02, 0.025), 
           legend=as.expression(sapply(inversion$result$rs, 
                                       function (rr) bquote(r[0]~"="~.(rr)))))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp,#+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    title(ylab=bquote("Convolution"~bgroup("(", 
                                           kappa^(.(model$f1.exp)*','~
                                                      .(model$f2.exp))*'*'*
                                               frac(d*.(model$f1.exp),
                                                    .(model$f1.exp)),
                                           ")")),
          line=2.5)
    
}

plot_response <- function(model, inversion,
                             ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                             font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 1.5, 0, 0))
    
    d.f1.spl <- splinefun(model$k1$x[-1], model$d.f1.true)
    #convs <- apply(inversion$avg_kerns, 2, function(avg_kern) {
    #    model$d.f1.true * avg_kern[-1]
    #})
    convs <- apply(inversion$avg_kerns, 2, function(avg_kern) {
        sintegral(model$k1$x[-1], model$d.f1.true * avg_kern[-1])$value
    })
    
    plot(NA, xaxs='i',
         ylim=range(convs, 
                    model$d.f1.true[model$k1$x[-1] >= min(inversion$result$rs) &
                                    model$k1$x[-1] <= max(inversion$result$rs)]), 
         xlim=c(0, 0.4),#range(inversion$result$rs)*1.25, 
         axes=F, 
         xlab="Radius r/R", 
         ylab="")
    
    abline(h=0, lty=2, col='gray')
    
    lines(model$r, model$d.f1.true, type='l', lty=2, col='gray', lwd=3)
    
    col.pal <- colorRampPalette(c(blue, 'black', red))(length(convs))#ncol(convs))
    points(inversion$result$rs, convs, col=col.pal, pch=20, cex=1.5)
    #for (ii in 1:ncol(convs))
    #    lines(model$k1$x[-1], convs[,ii], lty=ii, lwd=1.5,
    #          col=col.pal[ii])
    
    legend("bottomright", #lty=1:ncol(convs), lwd=1.5, 
           pch=20, 
           col=col.pal[1:length(convs)],#ncol(convs)],
           inset=c(0.02, 0.025), 
           legend=as.expression(sapply(inversion$result$rs, 
                                       function (rr) bquote(r[0]~"="~.(rr)))))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp,#+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    title(ylab=bquote(integral(kappa^(.(model$f1.exp)*','~.(model$f2.exp))%.%
                               frac(d*.(model$f1.exp), .(model$f1.exp))~
                                   dr,0,R)),
          line=2.5)
    
}

plot_avg_cross <- function(model, inversion,
                          ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                          font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 1.5, 0, 0))
    
    d.f1.spl <- splinefun(model$k1$x[-1], model$d.f1.true)
    #convs <- apply(inversion$avg_kerns, 2, function(avg_kern) {
    #    model$d.f1.true * avg_kern[-1]
    #})
    convs <- sapply(1:ncol(inversion$avg_kerns), function(ii) {
        avg_kern <- inversion$avg_kerns[,ii]
        cross_kern <- inversion$cross_kerns[,ii]
        sintegral(model$k1$x[-1], model$d.f1.true * avg_kern[-1] *
                      cross_kern[-1])$value
    })
    
    plot(NA, xaxs='i',
         ylim=range(convs, 
                    model$d.f1.true[model$k1$x[-1] >= min(inversion$result$rs) &
                                    model$k1$x[-1] <= max(inversion$result$rs)]), 
         xlim=c(0, 0.4),#range(inversion$result$rs)*1.25, 
         axes=F, 
         xlab="Radius r/R", 
         ylab="")
    
    abline(h=0, lty=2, col='gray')
    
    lines(model$r, model$d.f1.true, type='l', lty=2, col='gray', lwd=3)
    
    col.pal <- colorRampPalette(c(blue, 'black', red))(length(convs))#ncol(convs))
    points(inversion$result$rs, convs, col=col.pal, pch=20, cex=1.5)
    #for (ii in 1:ncol(convs))
    #    lines(model$k1$x[-1], convs[,ii], lty=ii, lwd=1.5,
    #          col=col.pal[ii])
    
    legend("bottomright", #lty=1:ncol(convs), lwd=1.5, 
           pch=20, 
           col=col.pal[1:length(convs)],#ncol(convs)],
           inset=c(0.02, 0.025), 
           legend=as.expression(sapply(inversion$result$rs, 
                                       function (rr) bquote(r[0]~"="~.(rr)))))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp,#+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    title(ylab=bquote(integral(kappa^(.(model$f1.exp)*','~.(model$f2.exp))%.%
                                   frac(d*.(model$f1.exp), .(model$f1.exp))~
                                   dr,0,R)),
          line=2.5)
    
}

plot_error_corr <- function(model, kernels, inversion, cross=F,
                            ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                            font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, -0.5, 0, 0))
    
    cs <- inversion$inv.coefs
    rs <- inversion$result$rs
    
    r_i <- 1
    E <- sapply(1:length(rs), function(r_i) { 
        sapply(1:ncol(cs), function(r_j) {
            numerator <- sum(sapply(1:nrow(cs), function(row) {
                ln <- get_ln(model$modes[row])
                dnu <- model$nus[model$nus$l==ln$l & 
                                 model$nus$n==ln$n,]$d.r.diff
                cs[row,r_i] * cs[row,r_j] * dnu**2
            }))
            denominator <- sum(sapply(1:nrow(cs), function(row) {
                ln <- get_ln(model$modes[row])
                dnu <- model$nus[model$nus$l==ln$l & 
                                 model$nus$n==ln$n,]$d.r.diff
                sqrt(cs[row,r_i]**2 * dnu**2) * sqrt(cs[row,r_j]**2 * dnu**2)
            }))
            numerator / denominator
        })
    })
    
    plot(NA, axes=F,
         ylim=range(E), 
         xlim=range(rs),
         xlab="Radius r/R",
         ylab="Error Correlation")
    
    abline(h=0, lty=2, col='gray')
    
    col.pal <- colorRampPalette(c(blue, 'black', red))(ncol(E))
    xs <- seq(min(rs), max(rs), 0.001)
    for (ii in 1:ncol(E)) {
        lines(rs, E[,ii], 
              #xs, predict(loess(E[,ii]~rs, degree=1, span=0.3), xs),
              #smooth.spline(rs, E[,ii], spar=0.1),
              lty=ii, lwd=1.5,
              col=col.pal[ii])
        points(rs[ii], E[ii,ii], pch=20, col=col.pal[ii])
    }
    
    legend("bottomleft", lty=1:ncol(convs), lwd=1.5, 
           col=col.pal[1:ncol(convs)],
           inset=c(0.02, 0.025), 
           legend=as.expression(sapply(inversion$result$rs, 
                                       function (rr) bquote(r[0]~"="~.(rr)))))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp,#+c(1,0,0),
            family=font, cex.axis=text.cex)
    
}

plot_function <- function(model, inversion,
                          ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                          font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    result <- inversion$result
    
    xlim <- c(max(0, min(result$fwhm.left) - 0.01),
              min(1, max(result$fwhm.right) + 0.01))
    
    inside.rs <- model$r>xlim[1] & model$r<xlim[2]
    
    ylim <- range(if (!is.null(model$m2.f1)) model$m2.f1[inside.rs] else 1, 
             model$f1[inside.rs],
             result$f+result$err, result$f-result$err)
    
    plot(NA, axes=F, 
         xlim=xlim, 
         ylim=ylim,
         xlab="Radius r/R",
         ylab="")
    
    if (!is.null(model$m2.f1)) { 
        lines(model$r, model$m2.f1, type='l', lty=2, col='gray', lwd=3) 
    } 
    lines(model$r, model$f1, lty=3, col=blue, lwd=3)
    
    abline(h=0, lty=2, col='gray')
    weights <- rep(0.5, nrow(result))
    for (row_i in 1:nrow(result)) { 
        row <- result[row_i,]
        with(row, segments(fwhm.left, f, fwhm.right, f, 
                           col=adjustcolor('darkred',alpha.f=weights[row_i])))
        with(row, segments(fwhm.mid, f-f.err, fwhm.mid, f+f.err,
                           col=adjustcolor('darkred',alpha.f=weights[row_i])))
        with(row, points(fwhm.mid, f, col='darkred', pch=20, cex=0.5))
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend("topright", lty=c(2, 3, NA), 
           col=c('darkgray', blue, 'darkred'), lwd=2, 
           inset=c(0.01, 0.015),
           pch=c(NA, NA, 3), 
           legend=c("Star", "Model", "Inversion"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.name <- model$f1.name
    f1.exp <- model$f1.exp
    title(ylab=bquote(.(f1.name)~.(f1.exp)))
}

plot_freq_diff <- function(model, 
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    plot(NA, axes=F, 
         xlim=range(model$nus$nu.y), 
         ylim=with(model$nus, 
             range(nu.x-nu.y-dnu, nu.x-nu.y+dnu)),
             #range(r.diff+d.r.diff, model$r.diff-d.r.diff)), 
         xlab=expression("Frequency"~nu[T]/mu*Hz), 
         ylab="")
    
    abline(h=0, lty=2, col='#00000099')
    
    col.pal <- c('black', blue, red, '#ACA4E2')
    with(model$nus, 
        segments(nu.y, nu.x-nu.y-dnu, nu.y, nu.x-nu.y+dnu,
        #segments(nu.y, r.diff-d.r.diff, nu.y, r.diff+d.r.diff, 
        col=adjustcolor(col.pal[l], alpha.f=0.5)))
    with(model$nus, points(nu.y, 
        nu.x-nu.y, 
        #r.diff, 
        pch=l+1, col=col.pal[l], cex=0.5))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    #legend("topright", lty=c(2, 3, NA), 
    #       col=c('darkgray', blue, 'darkred'), lwd=2, 
    #       inset=c(0.01, 0.015),
    #       pch=c(NA, NA, 3), 
    #       legend=c("Star", "Model", "Inversion"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=expression("Frequency difference"~(nu[M]-nu[T])))#/nu[T]))
}

plot_freq_diffs <- function(model.list, 
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        with(model.list[[ii]]$nus[model.list[[ii]]$nus$l == 0,], 
            range((nu.x-nu.y-dnu)/nu.y, (nu.x-nu.y+dnu)/nu.y))
    }))
    
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    pch <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 20
    })
    
    ordering <- rev(order(sapply(1:length(model.list), function(ii) {
        with(model.list[[ii]]$nus[model.list[[ii]]$nus$l == 0,], nu.x-nu.y)[1]
    })))
    
    plot(NA, axes=F, 
         xlim=range(model.list[[1]]$nus[model.list[[1]]$nus$l==0,]$nu.y), 
         ylim=ylim,
         xlab=expression("Radial mode frequency"~nu/mu*Hz), 
         ylab="")
    
    abline(h=0, lty=2, col='#00000099')
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        #with(model$nus[model$nus$l == 0,],
        #    segments(nu.y, nu.x-nu.y-dnu, nu.y, nu.x-nu.y+dnu))
        with(model$nus[model$nus$l == 0,], 
            points(nu.y, (nu.x-nu.y)/nu.y, col=color[ii], pch=pch[ii]))
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend('right', pch=pch[ordering], col=color[ordering], inset=c(0.125,0),
        legend=c(
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )[ordering]
    )
    
    #legend("topright", lty=c(2, 3, NA), 
    #       col=c('darkgray', blue, 'darkred'), lwd=2, 
    #       inset=c(0.01, 0.015),
    #       pch=c(NA, NA, 3), 
    #       legend=c("Star", "Model", "Inversion"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=expression("Frequency difference"~(nu['ref']-nu)/nu))#/nu[T]))
}

plot_freq_diffs.minus_offset <- function(model.list, 
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        nus <- model$nus
        w.mean <- with(nus, weighted.mean((nu.x-nu.y)/nu.y, 1/dnu))
        with(nus[nus$l == 0,], range((nu.x-nu.y-dnu)/nu.y - w.mean, 
                                     (nu.x-nu.y+dnu)/nu.y - w.mean))
    }))
    
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    pch <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 20
    })
    
    ordering <- rev(order(sapply(1:length(model.list), function(ii) {
        with(model.list[[ii]]$nus[model.list[[ii]]$nus$l == 0,], nu.x-nu.y)[1]
    })))
    
    plot(NA, axes=F, 
         xlim=range(model.list[[1]]$nus[model.list[[1]]$nus$l==0,]$nu.y), 
         ylim=ylim,
         xlab=expression("Radial mode frequency"~nu/mu*Hz), 
         ylab="")
    
    abline(h=0, lty=2, col='#00000099')
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        nus <- model$nus
        w.mean <- with(nus, weighted.mean((nu.x-nu.y)/nu.y, 1/dnu))
        #with(model$nus[model$nus$l == 0,],
        #    segments(nu.y, nu.x-nu.y-dnu, nu.y, nu.x-nu.y+dnu))
        with(nus[model$nus$l == 0,], 
            points(nu.y, (nu.x-nu.y)/nu.y - w.mean, col=color[ii], pch=pch[ii]))
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend('topright', pch=pch[ordering], col=color[ordering], inset=c(0.125,0),
        legend=c(
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )[ordering]
    )
    
    #legend("topright", lty=c(2, 3, NA), 
    #       col=c('darkgray', blue, 'darkred'), lwd=2, 
    #       inset=c(0.01, 0.015),
    #       pch=c(NA, NA, 3), 
    #       legend=c("Star", "Model", "Inversion"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=expression("Frequency difference"~(nu['ref']-nu)/nu))#/nu[T]))
}

plot_dimless_freq_diffs <- function(model.list, 
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        with(model$nus[model$nus$l == 0,], 
            range((nu.x/sqrt(model$mass/model$radius^3) - 
                   nu.y/sqrt(solar.mass/solar.radius^3)) / 
                   (nu.y/sqrt(solar.mass/solar.radius^3))))
    }))
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    pch <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 20
    })
    
    ordering <- rev(order(sapply(1:length(model.list), function(ii) {
        with(model.list[[ii]]$nus[model.list[[ii]]$nus$l == 0,], 
            (nu.x-nu.y)/nu.y)[1]
    })))
    
    plot(NA, axes=F, 
         xlim=range(model.list[[1]]$nus[model.list[[1]]$nus$l==0,]$nu.y), 
         ylim=ylim,
         xlab=expression("Radial mode frequency"~nu/mu*Hz), 
         ylab="")
    
    abline(h=0, lty=2, col='#00000099')
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        #with(model$nus[model$nus$l == 0,],
        #    segments(nu.y, nu.x-nu.y-dnu, nu.y, nu.x-nu.y+dnu))
        with(model$nus[model$nus$l == 0,], 
            points(nu.y, (nu.x/sqrt(model$mass/model$radius^3)-
                          nu.y/sqrt(solar.mass/solar.radius^3))/
                          (nu.y/sqrt(solar.mass/solar.radius^3)), 
                   col=color[ii], pch=pch[ii]))
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend('topright', pch=pch[ordering], inset=c(0.01,0.01), 
        col=color[ordering], 
        legend=c(
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )[ordering]
    )
    
    #legend("topright", lty=c(2, 3, NA), 
    #       col=c('darkgray', blue, 'darkred'), lwd=2, 
    #       inset=c(0.01, 0.015),
    #       pch=c(NA, NA, 3), 
    #       legend=c("Star", "Model", "Inversion"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=expression("Dimensionless frequency difference"~
        (sigma['ref']-sigma)/sigma))#/nu[T]))
}

plot_dim_diffs <- function(model.list, 
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        with(model$nus[model$nus$l == 0,], 
            range((nu.x/sqrt(model$mass/model$radius^3) - 
                   nu.y/sqrt(solar.mass/solar.radius^3)) / 
                   (nu.y/sqrt(solar.mass/solar.radius^3))))
    }))
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    pch <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 20
    })
    
    ordering <- rev(order(sapply(1:length(model.list), function(ii) {
        with(model.list[[ii]]$nus[model.list[[ii]]$nus$l == 0,], 
            (nu.x-nu.y)/nu.y)[1]
    })))
    
    plot(NA, axes=F, 
         xlim=range(model.list[[1]]$nus[model.list[[1]]$nus$l==0,]$nu.y), 
         ylim=ylim,
         xlab=expression("Radial mode frequency"~nu/mu*Hz), 
         ylab="")
    
    abline(h=0, lty=2, col='#00000099')
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        #with(model$nus[model$nus$l == 0,],
        #    segments(nu.y, nu.x-nu.y-dnu, nu.y, nu.x-nu.y+dnu))
        with(model$nus[model$nus$l == 0,], 
            points(nu.y, (nu.x/sqrt(model$mass/model$radius^3)-
                          nu.y/sqrt(solar.mass/solar.radius^3))/
                          (nu.y/sqrt(solar.mass/solar.radius^3)), 
                   col=color[ii], pch=pch[ii]))
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend('topright', pch=pch[ordering], inset=c(0.01,0.01), 
        col=color[ordering], 
        legend=c(
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )[ordering]
    )
    
    #legend("topright", lty=c(2, 3, NA), 
    #       col=c('darkgray', blue, 'darkred'), lwd=2, 
    #       inset=c(0.01, 0.015),
    #       pch=c(NA, NA, 3), 
    #       legend=c("Star", "Model", "Inversion"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=expression("Dimensionless frequency difference"~
        (sigma['ref']-sigma)/sigma))#/nu[T]))
}


plot_difference <- function(model, inversion,
                            ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                            font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    result <- inversion$result
    
    xlim <- c(max(0, min(result$fwhm.left) - 0.01),
              min(1, max(result$fwhm.right) + 0.01))
    
    inside.rs <- model$r>xlim[1] & model$r<xlim[2]
    
    true.f1 <- model$m2.f1
    
    m2f <- model$m2.f1.spl(result$fwhm.mid)
    ylim <- range(model$f1[inside.rs] - true.f1[inside.rs],
                  result$f+result$f.err-m2f, 
                  result$f-result$f.err-m2f)
    
    plot(NA, axes=F, 
         xlim=xlim, 
         ylim=ylim, 
         xlab="Radius r/R", 
         ylab="")
    
    abline(h=0, lty=2, col='gray', lwd=3) 
    lines(model$r, model$f1-model$m2.f1, lty=3, col=blue, lwd=3)
    
    with(result, segments(fwhm.left, f-m2f, fwhm.right, f-m2f, 
             col=adjustcolor('darkred',alpha.f=0.5)))
    with(result, segments(fwhm.mid, f-f.err-m2f, fwhm.mid, f+f.err-m2f,
             col=adjustcolor('darkred',alpha.f=0.5)))
    with(result, points(fwhm.mid, f-m2f, col='darkred', pch=20, cex=0.5))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend("topright", lty=c(2, 3, NA), 
           col=c('darkgray', blue, 'darkred'), lwd=2, 
           inset=c(0.01, 0.015),
           pch=c(NA, NA, 3), 
           legend=c("Star", "Model", "Inversion"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.name <- model$f1.name
    f1.exp <- model$f1.exp
    title(ylab=bquote(.(f1.name)~'difference'~(.(f1.exp)['ref']-.(f1.exp))))
}

plot_function_differences <- function(model.list, inversion.list,
                             ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                             font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    xlim <- range(sapply(1:length(model.list), function(ii) {
        inversion <- inversion.list[[ii]]
        c(max(0, min(inversion$result$fwhm.left) - 0.01),
          min(1, max(inversion$result$fwhm.right) + 0.01))
    }))
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        inversion <- inversion.list[[ii]]
        result <- inversion$result
        m2f <- model$m2.f1.spl(result$fwhm.mid)
        inside.rs <- model$r>xlim[1] & model$r<xlim[2]
        range((model$f1[inside.rs] - model$m2.f1[inside.rs])/
            model$m2.f1[inside.rs]) 
    }))
    
    
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    lty <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 1
    })
    
    plot(NA, axes=F, 
         xlim=xlim, 
         ylim=ylim, 
         xlab="Radius r/R", 
         ylab="")
    
    abline(h=0, lty=2, col='gray', lwd=3) 
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        result <- inversion.list[[ii]]$result
        #pch <-   if (model$M < 1) 2 else if (model$M > 1) 3 else 1
        #color <- if (model$R < 0.99) blue else if (model$R > 1.01) red else 
        #    'black'
        lines(model$r, (model$f1-model$m2.f1)/model$m2.f1, lty=lty[ii], 
            col=adjustcolor(color[ii], alpha.f=0.75), lwd=2)
        #xs <- seq(0, 1, 0.001)
        #diffs <- splinefun(model$r, model$f1-model$m2.f1)(xs)
        #points(xs, diffs, lty=3, col=color, pch=pch, cex=0.33)
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend('bottomright', lty=lty, inset=c(0.01,0.01), col=color, 
        legend=c(
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )
    )
    
    #legend('bottomright', inset=c(0.01, 0.01), 
    #    lty=c(2, 1, 3, 1,1,1), 
    #    col=c('black', 'black', 'black', blue, 'black', red),
    #    legend=c(
    #        expression(M==0.984),
    #        expression(M==1),
    #        expression(M==1.016),
    #        expression(R==0.98),
    #        expression(R==1),
    #        expression(R==1.02)))
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.name <- model$f1.name
    f1.exp <- model$f1.exp
    title(ylab=bquote(.(f1.name)~'difference'~
        (.(f1.exp)['ref']-.(f1.exp))/.(f1.exp)))
}

plot_inv_diffs <- function(model.list, inversion.list,
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    xlim <- range(sapply(1:length(model.list), function(ii) {
        inversion <- inversion.list[[ii]]
        c(max(0, min(inversion$result$fwhm.left) - 0.01),
          min(1, max(inversion$result$fwhm.right) + 0.01))
    }))
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        result <- inversion.list[[ii]]$result
        m2f <- model$m2.f1.spl(result$fwhm.mid)
        with(result, range((f+f.err-m2f)/m2f, (f-f.err-m2f)/m2f)) 
    }))
    
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    pch <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 20
    })
    
    plot(NA, axes=F, 
         xlim=xlim+c(0, 0.125), 
         ylim=ylim, 
         xlab="Radius r/R", 
         ylab="")
    
    abline(h=0, lty=2, col='gray', lwd=3) 
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        result <- inversion.list[[ii]]$result
        #lines(model$r, model$f1-model$m2.f1, lty=3, col=blue, lwd=3)
        m2f <- model$m2.f1.spl(result$fwhm.mid)
        with(result, 
             segments(fwhm.left, (f-m2f)/m2f, fwhm.right, (f-m2f)/m2f, 
                      col=adjustcolor(color[ii], alpha.f=0.6), lend='butt'))
        with(result, 
             segments(fwhm.mid, (f-f.err-m2f)/m2f, fwhm.mid, (f+f.err-m2f)/m2f,
                      col=adjustcolor(color[ii], alpha.f=0.6), lend='butt'))
    }
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        result <- inversion.list[[ii]]$result
        m2f <- model$m2.f1.spl(result$fwhm.mid)
        with(result, 
            points(fwhm.mid, (f-m2f)/m2f, 
                   col=color[ii], pch=pch[ii], cex=1, lwd=2))
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend('right', pch=pch, inset=c(0.01,0.01), col=color, 
        legend=c(
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )
    )
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.name <- model$f1.name
    f1.exp <- model$f1.exp
    title(ylab=bquote(.(f1.name)~'difference'~
        (.(f1.exp)['inv']-.(f1.exp))/.(f1.exp)))
}

plot_inv_diffs_lines <- function(model.list, inversion.list,
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    xlim <- range(sapply(1:length(model.list), function(ii) {
        inversion <- inversion.list[[ii]]
        range(inversion$result$fwhm.mid)
        #c(max(0, min(inversion$result$fwhm.left) - 0.01),
        #  min(1, max(inversion$result$fwhm.right) + 0.01))
    }))
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        result <- inversion.list[[ii]]$result
        m2f <- model$m2.f1.spl(result$fwhm.mid)
        range((result$f-m2f)/m2f)
    }))
    
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    pch <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 20
    })
    lty <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 1
    })
    
    plot(NA, axes=F, 
         xlim=xlim+c(-0.01, 0.125), 
         ylim=ylim, 
         xlab="Radius r/R", 
         ylab="")
    
    abline(h=0, lty=2, col='gray', lwd=3) 
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        result <- inversion.list[[ii]]$result
        m2f <- model$m2.f1.spl(result$fwhm.mid)
        with(result, 
            lines(fwhm.mid, (f-m2f)/m2f, 
                  col=color[ii], lty=lty[ii], cex=1, lwd=2))
        with(result, 
            points(fwhm.mid, (f-m2f)/m2f, 
                   col=color[ii], pch=pch[ii], cex=1, lwd=2))
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend('right', pch=pch, lty=lty, inset=c(0.01,0.01), col=color, 
        legend=c(
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )
    )
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.name <- model$f1.name
    f1.exp <- model$f1.exp
    title(ylab=bquote(.(f1.name)~'difference'~
        (.(f1.exp)['inv']-.(f1.exp))/.(f1.exp)))
}

plot_inv_diffs_mean <- function(model.list, inversion.list,
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        result <- inversion.list[[ii]]$result
        m2f <- model$m2.f1.spl(result$fwhm.mid)
        range((result$f-m2f)/m2f)
    }))
    
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    pch <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 4 else 20
    })
    lty <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 1
    })
    
    f.stds <- sapply(inversion.list, function(inversion) inversion$result$f.err)
    f.weights <- (1/f.stds) / apply(1/f.stds, 1, sum)
    f.means <- sapply(inversion.list, function(inversion) inversion$result$f)
    w.f.means <- sapply(1:nrow(f.weights), 
        function(ii) weighted.mean(f.means[ii,], f.weights[ii,]))
    
    #w.f.stds <- apply(f.means, 1, sd)
    
    w.f.stds <- apply(replicate(1000, {
            f.means2 <- sapply(1:ncol(f.means), function(ii) 
                rnorm(length(f.means[,ii]), f.means[,ii], f.stds[,ii]))
            apply(f.means2, 1, sd)
        }), 1, mean)
    
    
    #w.f.stds <- apply(replicate(1000, {
    #        f.means2 <- sapply(1:ncol(f.means), function(ii) 
    #            rnorm(length(f.means[,ii]), f.means[,ii], f.stds[,ii]))
    #        sapply(1:nrow(f.weights), 
    #            function(ii) weighted.mean(f.means2[ii,], f.weights[ii,]))
    #    }), 1, sd)
    
    
    #w.f.stds <- sapply(1:nrow(f.weights), 
    #    #function(ii) sqrt(sum(f.weights[ii,]**2*f.stds[ii,]**2)))
    #    function(ii) sqrt(nrow(f.weights) * sum(1/f.stds[ii,])**-2))
    #    #(f.means[ii,]-w.f.means[ii])**2)))
    
    fwhm.left <- apply(sapply(inversion.list, function(inversion) 
        with(inversion$result, fwhm.left)), 1, median)
    fwhm.right <- apply(sapply(inversion.list, function(inversion) 
        with(inversion$result, fwhm.right)), 1, median)
    fwhm.mid <- apply(sapply(inversion.list, function(inversion) 
        with(inversion$result, fwhm.mid)), 1, mean)
    
    m2f <- model.list[[1]]$m2.f1.spl(fwhm.mid)
    
    plot(NA, axes=F, 
         xlim=c(0.05, 0.5),
              #range(fwhm.left, fwhm.right)+c(0, 0.15), 
         ylim=c(-0.08, 0.08),
              #range(ylim, (w.f.means-w.f.stds-m2f)/m2f, 
              #            (w.f.means+w.f.stds-m2f)/m2f), 
         xlab="Radius r/R", 
         ylab="")
    
    abline(h=0, lty=2, col='gray', lwd=3) 
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        result <- inversion.list[[ii]]$result
        m2f <- model$m2.f1.spl(result$fwhm.mid)
        with(result, 
            lines(fwhm.mid, (f-m2f)/m2f, lty=lty[ii], cex=1, lwd=2,
                  col=adjustcolor(color[ii], alpha.f=0.3)))
        with(result, 
            points(fwhm.mid, (f-m2f)/m2f, cex=1, lwd=2,
                   col=adjustcolor(color[ii], alpha.f=0.3), pch=pch[ii]))
    }
    
    segments(fwhm.left,  (w.f.means-m2f)/m2f, 
             fwhm.right, (w.f.means-m2f)/m2f, lwd=2)
    segments(fwhm.mid, (w.f.means-w.f.stds-m2f)/m2f, 
             fwhm.mid, (w.f.means+w.f.stds-m2f)/m2f, lwd=2)
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    legend('right', pch=c(NA, 3, pch), lty=c(2, NA, lty), 
        inset=c(0.01,0.01), col=c('gray', 1, color), lwd=c(3, 2, rep(1, 9)),
        legend=c(
            "True Profile", 
            "Mean Inversion Result", 
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )
    )
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.name <- model$f1.name
    f1.exp <- model$f1.exp
    title(ylab=bquote(.(f1.name)~'difference'~
        (.(f1.exp)['inv']-.(f1.exp))/.(f1.exp)))
}


plot_one_surfless <- function(model, k.pair, nondimensionalize=F, ...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font,
        legend.spot='left') {
    real.modes <- model$nus$nu.x > 0
    nus <- model$nus[real.modes,]
    modes <- model$modes[real.modes]
    
    if (nondimensionalize) {
        nu.x <- nus$nu.x / sqrt(cgrav * 
            solar_mass * model$M / (solar_radius * model$R)**3)
        nu.y <- nus$nu.y / sqrt(cgrav * 
            solar_mass * (model$M - model$dM) / 
            (solar_radius * (model$R - model$dR))**3)
        nus$r.diff <- (nu.x - nu.y) / nu.x
    }
    
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab="", 
        xlim=range(nus$nu.x),
        ylim=range(nus$k.diffs, nus$r.diff) + 
            max(abs(range(nus$k.diffs, nus$r.diff))) * c(-0.1, 0.1))
    abline(h=0, lty=2)
    segments(nus$nu.x, nus$k.diffs, nus$nu.x, nus$r.diff, 
        col=adjustcolor('black',alpha.f=0.25), lty=3)
    pch <- if (4 %in% nus$l) 20 else (nus$l+1)
    points(nus$nu.x, nus$k.diffs, col=blue, cex=0.33, pch=pch)
    points(nus$nu.x, nus$r.diff, col='darkred', cex=0.33, pch=pch)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=1, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    #legend(legend.spot, col=c(1,1,1,1,1, blue, 'darkred'), #bty='n',
    #    inset=c(0.01, 0.01),
    #    pch=c(NA, 1:4, 20, 20), cex=text.cex, legend=as.expression(c(
    #        bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))),
    #        expression(l==0), expression(l==1), 
    #        expression(l==2), expression(l==3),
    #        bquote('Kernels'),
    #        bquote('Exact'))))
    legend(legend.spot, col=c(1,1,1,1, blue, 'darkred'), #bty='n',
        inset=c(0.01, 0.01), 
        pch=c(1:4, 20, 20), cex=0.8*text.cex, legend=as.expression(c(
            expression(l==0), expression(l==1), 
            expression(l==2), expression(l==3),
            bquote(K^(.(k.pair$f1.exp)*','~.(k.pair$f2.exp))),
            bquote('Exact'))))
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=bquote( delta*nu/nu ))
}

plot_one_Qnl <- function(model, k.pair, ...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font,
        legend.spot='left') {
    real.modes <- model$nus$nu.x > 0
    nus <- model$nus[real.modes,]
    modes <- model$modes[real.modes]
    
    k.diffs <- nus$k.diffs / nus$Q
    r.diff <- nus$r.diff / nus$Q
    
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab="", 
        xlim=range(nus$nu.x),
        ylim=range(k.diffs, r.diff))
    abline(h=0, lty=2)
    segments(nus$nu.x, k.diffs, nus$nu.x, r.diff, 
        col=adjustcolor('black',alpha.f=0.25), lty=3)
    pch <- if (4 %in% nus$l) 20 else (nus$l+1)
    points(nus$nu.x, k.diffs, col=blue, cex=0.33, pch=pch)
    points(nus$nu.x, r.diff, col='darkred', cex=0.33, pch=pch)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=1, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    #legend(legend.spot, col=c(1,1,1,1,1, blue, 'darkred'), #bty='n',
    #    inset=c(0.01, 0.01),
    #    pch=c(NA, 1:4, 20, 20), cex=text.cex, legend=as.expression(c(
    #        bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))),
    #        expression(l==0), expression(l==1), 
    #        expression(l==2), expression(l==3),
    #        bquote('Kernels'),
    #        bquote('Exact'))))
    legend(legend.spot, col=c(1,1,1,1, blue, 'darkred'), #bty='n',
        inset=c(0.01, 0.01), 
        pch=c(1:4, 20, 20), cex=0.8*text.cex, legend=as.expression(c(
            expression(l==0), expression(l==1), 
            expression(l==2), expression(l==3),
            bquote(K^(.(k.pair$f1.exp)*','~.(k.pair$f2.exp))),
            bquote('Exact'))))
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=bquote( #'Relative frequency difference' ~ 
        delta*nu/nu / Q[n*','*l] ))
}


plot_kernel_diffs <- function(model, k.pair, nondimensionalize=F, ...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font,
        legend.spot='left') {
    real.modes <- model$nus$nu.x > 0
    nus <- model$nus[real.modes,]
    modes <- model$modes[real.modes]
    
    if (nondimensionalize) {
        nu.x <- nus$nu.x / sqrt(cgrav * 
            solar_mass * model$M / (solar_radius * model$R)**3)
        nu.y <- nus$nu.y / sqrt(cgrav * 
            solar_mass * (model$M - model$dM) / 
            (solar_radius * (model$R - model$dR))**3)
        nus$r.diff <- (nu.x - nu.y) / nu.x
    }
    
    diffs <- nus$r.diff - nus$k.diffs
    col.pal <- c("#ca0020", "#f4a582", "#0571b0", "#800080")
    
    xs <- seq(floor(min(nus$nu.x)), ceil(max(nus$nu.x)), 1)
    ys <- splinefun(nus$nu.x, diffs)(xs)
    
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 1.5, 0, 0))
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab="", 
        xlim=range(nus$nu.x),
        ylim=range(diffs))#, 0))
    abline(h=0, lty=2)
    #segments(nus$nu.x, k.diffs, nus$nu.x, nus$r.diff, 
    #    col=adjustcolor('black',alpha.f=0.25), lty=3)
    pch <- if (4 %in% nus$l) 20 else (nus$l+1)
    lines(xs, ys, col=adjustcolor('black',alpha.f=0.25), lty=3)
    points(nus$nu.x, diffs, col=col.pal[nus$l+1], pch=pch)
    #points(nus$nu.x, nus$r.diff, col='darkred', cex=0.33, pch=pch)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=1, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    legend(legend.spot, col=c(NA, col.pal), #bty='n',
        inset=c(0.02, 0.02), pch=c(NA, 1:4), cex=0.8*text.cex, 
        legend=as.expression(c(
            bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))),
            expression(l==0), expression(l==1), 
            expression(l==2), expression(l==3))))
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=bquote( #'Exact'-'Kernel' ~ 
            frac(delta*nu, nu) - 
            "<"*bold(K)%.%bold(frac(d*f, f))*">"))
    #title(xlab=expression('Frequency'~nu/mu*Hz))
}


plot_kernel_diffs_surf <- function(model, k.pair, ...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font,
        legend.spot='left') {
    real.modes <- model$nus$nu.x > 0
    nus <- model$nus[real.modes,]
    modes <- model$modes[real.modes]
    
    diffs <- nus$r.diff - nus$k.diffs
    
    #nu <- model$nus$nu.x / model$nu_ac
    #r.diff <- ( nu - m1$nus$nu.y/sqrt(m1$M/R**3) )/nu #m1$nus$r.diff #
    #inertia <- m1$nu_ac * m1$nus$Q_norm #* m1$nus$d.r.diff
    nu <- nus$nu.x / model$nu_ac
    inertia <- nus$E #Q# 
    Xpinv <- ginv( matrix(c(nu**-2, nu**2) / inertia, ncol=2) )
    a.r.1 <- Xpinv %*% ( diffs ) #/ m1$nus$d.r.diff )
    F_surf <- ( a.r.1[[1]]*nu**-2 + 
                a.r.1[[2]]*nu**2 ) / inertia
    
    #surf.lm <- lm(diffs~I(nu^2)+I(nu^-2)+I(nu^3)+I(nu^-3)+nu+I(nu^-1))
    
    diffs <- diffs - F_surf
    
    col.pal <- c("#ca0020", "#f4a582", "#0571b0", "#800080")
    
    xs <- seq(floor(min(nus$nu.x)), ceil(max(nus$nu.x)), 1)
    ys <- splinefun(nus$nu.x, diffs)(xs)
    
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 1.5, 0, 0))
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab="", 
        xlim=range(nus$nu.x),
        ylim=range(diffs))#, 0))
    abline(h=0, lty=2)
    #segments(nus$nu.x, k.diffs, nus$nu.x, nus$r.diff, 
    #    col=adjustcolor('black',alpha.f=0.25), lty=3)
    pch <- if (4 %in% nus$l) 20 else (nus$l+1)
    lines(xs, ys, col=adjustcolor('black',alpha.f=0.25), lty=3)
    points(nus$nu.x, diffs, col=col.pal[nus$l+1], pch=pch)
    #points(nus$nu.x, nus$r.diff, col='darkred', cex=0.33, pch=pch)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=1, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    legend(legend.spot, col=c(NA, col.pal), #bty='n',
        inset=c(0.02, 0.02), pch=c(NA, 1:4), cex=0.8*text.cex, 
        legend=as.expression(c(
            bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))),
            expression(l==0), expression(l==1), 
            expression(l==2), expression(l==3))))
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=bquote( #'Exact'-'Kernel' ~ 
        #bgroup("(", 
            frac(delta*nu, nu) - 
            "<"*bold(K)%.%bold(frac(d*f, f))*">"))
            #integral(bold(K)%.%bold(frac(d*f, f)),0,R)~dr, 
        #")") / mu*Hz ))
    #title(xlab=expression('Frequency'~nu/mu*Hz))
}

plot_param_hist <- function(param, xlab, use_fivenum=F, log='x', 
        majorn=5, minorn='auto', ...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font) {
    text.cex <- 1.5*text.cex
    par(cex.lab=text.cex, mar=mar+c(0.1, -0.2, 0, -0.2), mgp=mgp+c(0.4, 0, 0))
    if (log == 'x') param <- log10(param)
    if (use_fivenum) {
        fnum <- fivenum(param)
        param <- param[param>fnum[2] & param<fnum[4]]
    }
    hgram <- hist(param, plot=F)
    #xlim <- if (use_fivenum) {
    #    fnum <- fivenum(log10(param))
    #    c(ceil(fnum[2])-1, ceil(fnum[4])+1)
    #} else {
    xlim <- if (log=='x') {
        c(floor(min(hgram$mids))-1, ceil(max(hgram$mids))+1)
    } else c(0, 0.1)
    #}
    plot(hgram, main="", #hgram$mids, hgram$counts, type='S', 
        axes=F, xaxs='i', yaxs='i', col="#bfdbf7", 
        xlab=xlab,
        ylab="Count",
        xlim=xlim,
        ylim=c(0, max(hgram$counts))+1)
    #magaxis(1:4, labels=c(1,1,0,0), unlog='x')
    magaxis(side=1, tcl=-0.25, labels=T, unlog=log, mgp=mgp+c(0, 0.2, 0), 
        majorn=majorn, minorn=minorn, 
        las=1, family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=-0.25, labels=T, unlog=log, mgp=mgp+c(0, 0.25, 0), 
        las=1, family=font, cex.axis=text.cex)
    #magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0), unlog='x', 
    #    las=1, mgp=mgp-c(0.5,0.15,0), family=font, cex.axis=text.cex)
    #magaxis(side=2, tcl=0.25, labels=1, 
    #    las=1, mgp=mgp+c(1,0,0), family=font, cex.axis=text.cex)
}

plot_widths_hist <- function(widths, ...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font) {
    text.cex <- 1.5*text.cex
    par(cex.lab=text.cex, mar=mar+c(0.1, -0.2, 0, -0.2), mgp=mgp+c(0.4, 0, 0))
    hgram <- hist(widths, plot=F)
    plot(hgram, main="", #hgram$mids, hgram$counts, type='S', 
        axes=F, xaxs='i', yaxs='i', col="#bfdbf7",
        xlab=bquote("Kernel width"~Delta),
        ylab="Count",
        xlim=c(0, 0.1),
        ylim=c(0, max(hgram$counts))+2)
    magaxis(side=1, tcl=-0.25, labels=T, mgp=mgp+c(0, 0.2, 0), 
        las=1, family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=-0.25, labels=T, mgp=mgp+c(0, 0.25, 0), 
        las=1, family=font, cex.axis=text.cex)
    #magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0), 
    #    las=1, mgp=mgp-c(0.5,0.15,0), family=font, cex.axis=text.cex)
    #magaxis(side=2, tcl=0.25, labels=1, 
    #    las=1, mgp=mgp+c(1,0,0), family=font, cex.axis=text.cex)
}

plot_inversion_parameters <- function(inv.params, k.str=NULL, 
        wide=F, tall=F, ...) {
    make_plots(plot_param_hist, paste0('hist-mu', k.str),
        param=inv.params$error.sups, 
        xlab=bquote("Error suppression parameter"~mu), ...)
    make_plots(plot_param_hist, paste0('hist-beta', k.str),
        param=inv.params$cross.terms, use_fivenum=T, 
        xlab=bquote("Cross-term parameter"~beta), ...)
    make_plots(plot_widths_hist, paste0('hist-Delta', k.str),
        widths=inv.params$widths, ...)
}


plot_inversion_parameters_together <- function(inv.params, k.str=NULL, ...) {
    plot_all_params <- function(...) {
        par(mfrow=c(1,3))
        plot_param_hist(param=inv.params$error.sups, log='x', 
            xlab=bquote("Error suppression parameter"~mu), ...)
        plot_param_hist(param=inv.params$cross.terms, log='x', 
            use_fivenum=T, majorn=3, minorn=10,
            xlab=bquote("Cross-term parameter"~beta), ...)
        plot_param_hist(param=inv.params$widths, log='', 
            xlab=bquote("Kernel width"~Delta), ...)
        #plot_widths_hist(widths=inv.params$widths, ...)
    }
    make_plots(plot_all_params, paste0("hists", k.str), thin=F, tall=F)
}

