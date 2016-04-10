#### Plot random forest scores from subsetter.py 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(magicaxis)
library(RColorBrewer)
source(file.path('..', 'scripts', 'utils.R'))

out.dir <- file.path('plots', 'evaluate')

legend2 <- legend
body(legend2)[[c(38,4,6)]] <- quote({
    cidx <- col(matrix(1, ncol=ncol, nrow=n.legpercol))[1:n.leg]
    wc <- tapply(rep_len(w0,n.leg), cidx, max)
    w <- sum(wc) + 0.5 * xchar
    w0 <- rep(cumsum(c(0,wc[-length(wc)])), each=n.legpercol)[1:n.leg]
})
body(legend2)[[40]] <- quote(xt <- left + xchar + xextra + w0)

plot_accuracy <- function(DF, score='ev', variable='num_tracks', ..., 
        text.cex=1, mgp=utils.mgp, font="Palatino", plotlegend=T) {
    
    if ('log_g' %in% DF$variable) {
        DF <- DF[-which(DF$variable == 'log_g'),]
    }
    
    par(xpd=FALSE)
    ylim <- if (score=='sigma') c(0.25,1.5) else c(0, 1)
    plot(NA, axes=F, 
         xlim=range(DF[[variable]]), 
         ylim=ylim, 
         #c(1, 1-0.997),#0.001),#rev(range(1-DF[[score]])),#c(0,1), 
         log='x',#'xy', 
         xaxs='i', yaxs='i', 
         xlab=if (variable == 'num_tracks') "Number of Evolutionary Tracks"
              else if (variable == 'num_points') "Number of Models Per Track"
              else "Number of Trees in Forest",
         ylab="")
    
    minor <- 2 ^ (floor(log(sqrt(max(DF[[variable]]))) / log(2)))
    tick.places <- 2^(log2(min(DF[[variable]])):log2(max(DF[[variable]])))
    axis(1, at=tick.places, labels=tick.places, 
         tcl=-0.25, cex.axis=text.cex, tick=F)
    
    par(mgp=mgp+c(1, 0, 0.15))
    
    title(ylab=if (score == 'ev') {
        expression("Explained Variance"~V[e])
    } else if (score == 'r2') {
        expression("Coefficient of Determination"~R[2])
    } else if (score == 'sigma') {
        expression("Accuracy Per Precision"~abs(widehat(y)-y)/widehat(sigma))
    })
    
    if (score == "sigma") {
        #yticks <- 0:5
        yticks <- seq(min(ylim), max(ylim), 0.25)
        axis(2, tick=F, at=yticks, cex.axis=text.cex, las=1,
            labels=yticks)
    } else {
        #yticks <- c(1, 0.316, 0.1, 0.03, 0.01, 0.003)
        yticks <- c(0, 0.25, 0.5, 0.75, 1)
        axis(2, tick=F, at=yticks, cex.axis=text.cex, las=1,
            labels=c("0%", "25%", "50%", "75%", "100%"))
            #labels=c("0%", "68%", "90%", "97%", "99%", "99.7%"))
    }
    
    #grid(0, 5, lty=6, col="cornsilk2", equilogs=F) 
    for (ytick in yticks[2:(length(yticks)-1)]) 
        abline(h=ytick, col='cornsilk2')
    for (tick in tick.places[2:(length(tick.places)-1)]) 
        abline(v=tick, col='cornsilk2')
    
    par(xpd=TRUE)
    #variables <- unique(DF[['variable']])
    variables <- c('M', 'Y', 'Z', 'alpha', 'overshoot', 'diffusion',
                   'age', 'X_c', 'mass_cc', 'Y_surf', 'radius', 'L')
    col.pal <- colorRampPalette(brewer.pal(11, 'Spectral'))(length(variables))
    for (variable_i in 1:length(variables)) {
        average <- NULL
        for (num_tracks in unique(DF[[variable]])) {
            ss <- DF[DF[[variable]] == num_tracks &
                         DF$variable == variables[variable_i],]
            new.row <- ss[1,]
            new.row[[score]] <- median(ss[[score]])
            average <- rbind(average, new.row)
        }
        vals <- average[[score]]#ifelse(average[[score]] < 0, 0, average[[score]])
        vals.outside <- if (score=="sigma") vals>max(ylim) else vals<min(ylim)
        xs <- average[[variable]]
        if (any(vals.outside)) {
            #par(xpd=FALSE)
            
            #idx <- c(vals.outside, last+1)
            #lines(xs[idx], vals[idx], 
            #    col='black', lwd=3)
            #lines(xs[idx], vals[idx], 
            #    lwd=1.5, col=adjustcolor(col.pal[variable_i], alpha=0.95))
            
            #par(xpd=TRUE)
            last <- max(which(vals.outside))
            x.trial <- seq(xs[last], xs[last+1], length.out=10000)
            y.trial <- seq(vals[last], vals[last+1], length.out=10000)
            if (score=="sigma") {
                new.val <- max(ylim)
                new.x <- x.trial[min(which(y.trial <= new.val))]
            } else {
                new.val <- min(ylim)
                new.x <- x.trial[min(which(y.trial >= new.val))]
            }
            lines(c(new.x, xs[last+1]), 
                  c(new.val, vals[last+1]), 
                col='black', lwd=3)
            lines(c(new.x, xs[last+1]), 
                  c(new.val, vals[last+1]), 
                lwd=1.5, col=adjustcolor(col.pal[variable_i], alpha=0.95))
            points(new.x, new.val, 
                col='black', cex=0.75, pch=21+(variable_i%%3), 
                bg=adjustcolor(col.pal[variable_i], 
                               red.f=1.25, green.f=1.25, blue.f=1.25))
            vals <- vals[-c(1:last)]
            xs <- xs[-c(1:last)]
        }
        lines(xs, vals, col='black', lwd=3)
        lines(xs, vals, lwd=1.5,
            col=adjustcolor(col.pal[variable_i], alpha=0.95))
        points(xs, vals, col='black', cex=0.75, 
            pch=21+(variable_i%%3), 
            bg=adjustcolor(col.pal[variable_i], 
                           red.f=1.25, green.f=1.25, blue.f=1.25))
    }
    
    if (plotlegend) {
        shapes <- 1:length(variables)%%3 + 1
        labels <- as.expression(seis.labs[variables])
        text.widths <- rep(0.15, length(labels)/2)
        text.widths[3] <- 0.25
        text.widths[4] <- 0.35
        text.widths[5] <- 0.2
        legend2("topleft", legend=labels[1:6], horiz=T, bty='n', 
               inset=c(0, -0.16), cex=text.cex, text.col='white', 
               pch=c(16, 15, 18)[shapes], col=col.pal, 
               text.width=text.widths)
        
        legend2("topleft", legend=labels[1:6], horiz=T, bty='n', 
               inset=c(0, -0.16), cex=text.cex, 
               pch=c(1, 0, 5)[shapes], col='black', 
               text.width=text.widths)
        
        legend2("topleft", legend=labels[-1:-6], horiz=T, bty='n', 
               inset=c(0, -0.11), cex=text.cex, text.col='white', 
               pch=c(16, 15, 18)[shapes], col=col.pal[-1:-6], 
               text.width=text.widths)
        
        legend2("topleft", legend=labels[-1:-6], horiz=T, bty='n', 
               inset=c(0, -0.11), cex=text.cex, 
               pch=c(1, 0, 5)[shapes], col='black', 
               text.width=text.widths)
    }
}

make_plots(plot_accuracy, "num_tracks-ev", filepath=out.dir, 
    plotlegend=T, wide=F, mar=c(3, 3.75, 2.5, 1), 
    DF=read.table(file.path('subsetter', 'num_tracks.dat'),
                  header=T, stringsAsFactors=F))

make_plots(plot_accuracy, "num_tracks-sigma", filepath=out.dir, 
    score="sigma", plotlegend=T, wide=F, mar=c(3, 3.75, 2.5, 1), 
    DF=read.table(file.path('subsetter', 'num_tracks.dat'),
                  header=T, stringsAsFactors=F))

make_plots(plot_accuracy, "num_points", filepath=out.dir, 
    variable='num_points', plotlegend=T, wide=F, mar=c(3, 3.75, 2.5, 1), 
    #mar=c(3, 3.75, 1, 1), 
    DF=read.table(file.path('subsetter', 'num_points2.dat'), 
                  header=T, stringsAsFactors=F))

make_plots(plot_accuracy, "num_trees", filepath=out.dir, 
    variable='num_trees', plotlegend=T, wide=F, mar=c(3, 3.75, 2.5, 1), 
    #mar=c(3, 3.75, 2.5, 1), 
    DF=read.table(file.path('subsetter', 'num_trees2.dat'), 
                  header=T, stringsAsFactors=F))

