#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

#source('../scripts/utils.R') 
library(mgcv)

param.names <- c(bquote(beta), bquote(mu), bquote(Delta[f]))
param.strs <- c("beta", "mu", "deltaf")


plot_inversion <- function(model, inversion, k.pair,
        log='', xlim=NA, ylim=NA, legend.spot='topleft', 
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
        font="Times") {
    
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
    
    result <- inversion$result
    
    result <- result[with(result, 
        fwhm.left < rs & 
        fwhm.right > rs & 
        err < 0.2 &
        fwhm.right - fwhm.left < 0.2),]
    
    #if (length(xlim)<=1) xlim <- c(
    #    max(0, 0.85*min(result$fwhm.left)),
    #    min(1, 1.15*max(result$fwhm.right)))
    if (length(xlim)<=1) xlim <- range(result$fwhm.left, result$fwhm.right)
    if (length(ylim)<=1) ylim <- range(if (!is.null(d.f1.true)) {
                d.f1.true[model$r<max(result$rs) & 
                          model$r>min(result$rs) |
                          model$r<max(xlim) &
                          model$r>min(xlim)]
                } else 0, 
      result$df_dr+result$err, result$df_dr-result$err)
      #c(-0.15, 0.07),#
    
    plot(NA, axes=F, log=log,
         xlim=xlim, #c(0, 0.4),c(0.05, 0.35),#
         ylim=ylim,
         xlab="Radius r/R",
         ylab="")
    
    if (!is.null(d.f1.true)) { 
        lines(model$r, d.f1.true, type='l', lty=2, col='gray', lwd=3)
    }
    abline(h=0, lty=2, col='black')
    with(result, arrows(fwhm.left, df_dr, fwhm.right, df_dr, 
        code=3, angle=90, length=0.02,
        col=adjustcolor('darkred',alpha.f=0.5), lwd=1.5))
    with(result, arrows(fwhm.mid, df_dr-err, fwhm.mid, df_dr+err, 
        code=3, angle=90, length=0.02,
        col=adjustcolor('darkred',alpha.f=0.5), lwd=1.5))
    with(result, points(fwhm.mid, df_dr, col='darkred', pch=20, cex=0.5))
    
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
        legend=bquote(K^(.(k.pair$f1.exp)*','~.(k.pair$f2.exp))))
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.exp <- model$f1.exp
    title(ylab=bquote(#'Relative difference' ~ 
        delta*.(f1.exp)/.(f1.exp)))
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
    
    plot(NA, axes=F, log=log,
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

plot_inversion_lists_mean <- function(model, inversion.lists,
                           ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    get_q <- function(q_name, lists=inversion.lists) {
        sapply(lists, function(inversion.list) {
            results <- sapply(inversion.list, function(result) result[[q_name]])
            #mean(results)
            apply(results, 1, mean)
        })
    }
    
    w.f.means <- get_q('f')
    fwhm.left <- get_q('fwhm.left')
    fwhm.mid <- get_q('fwhm.mid')
    fwhm.right <- get_q('fwhm.right')
    
    w.f.stds <- sapply(inversion.lists, function(inversion.list) {
        results <- sapply(inversion.list, function(result) result$f)
        #sd(results)
        apply(results, 1, sd)
    })
    
    w.f.means <- apply(w.f.means, 1, mean)
    w.f.stds <- apply(w.f.stds, 1, mean)
    fwhm.left <- apply(fwhm.left, 1, mean)
    fwhm.mid <- apply(fwhm.mid, 1, mean)
    fwhm.right <- apply(fwhm.right, 1, mean)
    
    m.f <- model$f1.spl(fwhm.mid)
    
    
    for (r_i in length(inversion.lists[[1]][[1]]$rs)) {
        fs <- sapply(inversion.lists, function(inversion.list) 
            sapply(inversion.list, function(result) result$f[r_i]))
        errs <- fs <- sapply(inversion.lists, function(inversion.list) 
            sapply(inversion.list, function(result) result$f.err[r_i]))
    }
    
    plot(NA, axes=F, 
         xlim=c(0, 0.35),#range(fwhm.left, fwhm.right),
         ylim=c(-0.2, .1),
              #range((m.f - w.f.means-w.f.stds)/m.f, 
              #      (m.f - w.f.means+w.f.stds)/m.f),
              #      #c(-0.2, .1),
         xlab="Radius r/R",
         ylab="")
    
    d.f1.true <- if ('d.f1.true' %in% names(model)) model$d.f1.true else NULL
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
    title(ylab=bquote('Relative difference' ~ d*.(f1.exp)/.(f1.exp)))
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
           legend=c("Actual", "Mean Inversion"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    f1.exp <- model$f1.exp
    title(ylab=bquote('Relative difference' ~ d*.(f1.exp)/.(f1.exp)))
}

plot_kernels <- function(model, inversion, log='', xlim=NA,
         legend.spot='topright', ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
         font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    keep <- which(with(inversion$result, 
        fwhm.left < rs & 
        fwhm.right > rs & 
        err < 0.05 & 
        fwhm.right - fwhm.left < 0.15))
    
    kernels <- if (cross) inversion$cross_kerns else inversion$avg_kerns 
    kernels <- kernels[,keep]
        
    if (length(xlim)<=1) xlim <- c(0, 1)
    
    plot(NA, log=log, 
         type='l', lty=2, col='gray', lwd=3, 
         ylim=range(kernels[model$k1$x<0.9,]), 
         xlim=xlim,
         axes=F, 
         xlab="Radius r/R",
         ylab="")
    
    abline(h=0, lty=2, col='gray')
    for (kern_i in 1:ncol(kernels)) {
        if (cross) {
            lines(model$k2$x, kernels[,kern_i], lty=kern_i, lwd=2, col=kern_i)
        } else {
            lines(model$k1$x, kernels[,kern_i], lty=kern_i, lwd=2, col=kern_i)
        }
    }
    
    rs <- inversion$result$rs[keep]
    legend(legend.spot, col=1:length(rs), lty=1:length(rs), lwd=2,
           cex=0.8*text.cex,
           inset=c(if (cross) 0.15 else 0.01, 0.015), 
           legend=as.expression(sapply(rs, 
               function (rr) bquote(r[0]~"="~.(rr)))))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(1.5, 0, 0))
    if (cross) {
        title(ylab=bquote('Cross-Term Kernel'~
                              kappa^(.(model$f2.exp)*','~.(model$f1.exp))))
    } else {
        title(ylab=bquote('Averaging Kernel'~
                              kappa^(.(model$f1.exp)*','~.(model$f2.exp))))
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


plot_one_surfless <- function(model, k.pair, ...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font,
        legend.spot='left') {
    real.modes <- model$nus$nu.x > 0
    nus <- model$nus[real.modes,]
    modes <- model$modes[real.modes]
    
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab="", 
        xlim=range(nus$nu.x),
        ylim=range(nus$k.diffs, nus$r.diff))
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


plot_kernel_diffs <- function(model, k.pair, ...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font,
        legend.spot='left') {
    real.modes <- model$nus$nu.x > 0
    nus <- model$nus[real.modes,]
    modes <- model$modes[real.modes]
    
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
        bgroup("(", 
            frac(delta*nu, nu) - 
            "<"*bold(K)%.%bold(frac(d*f, f))*">", 
            #integral(bold(K)%.%bold(frac(d*f, f)),0,R)~dr, 
        ")") / mu*Hz ))
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
        bgroup("(", 
            frac(delta*nu, nu) - 
            "<"*bold(K)%.%bold(frac(d*f, f))*">", 
            #integral(bold(K)%.%bold(frac(d*f, f)),0,R)~dr, 
        ")") / mu*Hz ))
    #title(xlab=expression('Frequency'~nu/mu*Hz))
}

