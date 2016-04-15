#### Plot random forest scores from subsetter.py 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(magicaxis)
library(RColorBrewer)
source(file.path('..', 'scripts', 'utils.R'))

options(warn=2) 

out.dir <- file.path('plots', 'evaluate')

legend2 <- legend
body(legend2)[[c(38,4,6)]] <- quote({
    cidx <- col(matrix(1, ncol=ncol, nrow=n.legpercol))[1:n.leg]
    wc <- tapply(rep_len(w0,n.leg), cidx, max)
    w <- sum(wc) + 0.5 * xchar
    w0 <- rep(cumsum(c(0,wc[-length(wc)])), each=n.legpercol)[1:n.leg]
})
body(legend2)[[40]] <- quote(xt <- left + xchar + xextra + w0)

trans.black <- adjustcolor('black', alpha=0.5)

plot_accuracy <- function(DF, score='ev', variable='num_tracks', ..., 
        text.cex=1, mgp=utils.mgp, font="Palatino", plotlegend=T) {
    
    if (variable=='num_trees' && (score=='dist' || score=='sigma')) {
        DF <- DF[DF[[variable]] > 1,]
        #xlim[1] <- 2
    }
    
    if ('log_g' %in% DF$parameter) {
        DF <- DF[-which(DF$parameter == 'log_g'),]
    }
    if ('mass_cc' %in% DF$parameter) {
        DF <- DF[-which(DF$parameter == 'mass_cc'),]
    }
    
    par(xpd=FALSE)
    ylim <- if (score=='dist') {
        if (variable=='num_tracks') {
            c(0.25, 2) 
        } else if (variable=='num_points') {
            c(0.25, 1)
        } else if (variable=='num_trees') {
            c(0.25, 2)
        }
    } else c(0, 1)
    plot(NA, axes=F, 
         xlim=range(DF[[variable]]), 
         ylim=ylim, log='x',#'xy', 
         xaxs='i', yaxs='i', 
         xlab="",
         ylab="")
    
    minor <- 2 ^ (floor(log(sqrt(max(DF[[variable]]))) / log(2)))
    tick.places <- 2^(log2(min(DF[[variable]])):log2(max(DF[[variable]])))
    axis(1, at=tick.places, labels=tick.places, 
         tcl=-0.25, cex.axis=text.cex, tick=F)
    
    title(xlab=if (variable == 'num_tracks') {
        "Number of Evolutionary Tracks"
    } else if (variable == 'num_points') {
        "Number of Models Per Track"
    } else { 
        "Number of Trees in Forest"
    })
    par(mgp=mgp+c(1, 0, 0.15))
    title(ylab=if (score == 'ev') {
        expression("Explained Variance"~V[e])
    } else if (score == 'r2') {
        expression("Coefficient of Determination"~R[2])
    } else if (score == 'dist') {
        expression("Distance"~abs(widehat(y)-y)/widehat(sigma))
    } else if (score == 'sigma') {
        expression("Normalized Uncertainty"~widehat(sigma))
    } else if (score == 'diff') {
        expression("Normalized Error"~abs(widehat(y)-y))
    })
    
    if (score == "ev" || score == "r2") {
        yticks <- c(0, 0.25, 0.5, 0.75, 1)
        axis(2, tick=F, at=yticks, cex.axis=text.cex, las=1,
            labels=c("0%", "25%", "50%", "75%", "100%"))
    } else {
        yticks <- seq(min(ylim), max(ylim), ifelse(score=='dist', 0.5, 0.25))
        axis(2, tick=F, at=yticks, cex.axis=text.cex, las=1,
            labels=yticks)
    }
    
    for (ytick in yticks[2:(length(yticks)-1)]) 
        abline(h=ytick, col='cornsilk2')
    for (tick in tick.places[2:(length(tick.places)-1)]) 
        abline(v=tick, col='cornsilk2')
    
    par(xpd=TRUE)
    parameters <- c('M', 'Y', 'Z', 'alpha', 'overshoot', 'diffusion',
                   'age', 'X_c', 'Y_surf', 'radius', 'L')
    col.pal <- colorRampPalette(brewer.pal(11, 'Spectral'))(length(parameters))
    for (parameter_i in 1:length(parameters)) {
        average <- NULL
        for (value in unique(DF[[variable]])) {
            ss <- DF[DF[[variable]] == value &
                     DF$parameter == parameters[parameter_i],]
            new.row <- ss[1,]
            new.row[[score]] <- median(ss[[score]])
            average <- rbind(average, new.row)
        }
        vals <- average[[score]]
        if (score == 'diff' || score == 'sigma') {
            #param <- DF[DF$parameter==parameters[parameter_i],]
            #max.val <- max(param[[score]])
            #max.val <- max(average[[score]])
            #print(max.val)
            #vals <- normalize(vals)#vals / max.val
            vals <- normalize(vals)
        }
        vals.outside <- if (score == "ev" || score == "r2")
            vals<min(ylim) else vals>max(ylim) 
        xs <- average[[variable]]
        if (any(vals.outside)) {
            last <- max(which(vals.outside))
            if (!is.finite(vals[last])) vals[last] <- 2*max(ylim)
            x.trial <- seq(xs[last], xs[last+1], length.out=10000)
            y.trial <- seq(vals[last], vals[last+1], length.out=10000)
            if (score == "ev" || score == "r2") {
                new.val <- min(ylim)
                new.x <- x.trial[min(which(y.trial >= new.val))]
            } else {
                new.val <- max(ylim)
                new.x <- x.trial[min(which(y.trial <= new.val))]
            }
            if (parameter_i %in% 4:9) lines(c(new.x, xs[last+1]), 
                c(new.val, vals[last+1]), col=trans.black, lwd=3)
            lines(c(new.x, xs[last+1]), c(new.val, vals[last+1]), 
                lwd=1.5, col=adjustcolor(col.pal[parameter_i]))#, alpha=0.95))
            points(new.x, new.val, 
                col='black', cex=0.75, pch=21+(parameter_i%%3), 
                bg=adjustcolor(col.pal[parameter_i]))#, 
                               #red.f=1.25, green.f=1.25, blue.f=1.25))
            vals <- vals[-c(1:last)]
            xs <- xs[-c(1:last)]
        }
        if (parameter_i %in% 4:9) lines(xs, vals, col=trans.black, lwd=3)
        lines(xs, vals, lwd=1.5,
            col=adjustcolor(col.pal[parameter_i]))#, alpha=0.95))
        points(xs, vals, col='black', cex=0.75, 
            pch=21+(parameter_i%%3), 
            bg=adjustcolor(col.pal[parameter_i]))#, 
                           #red.f=1.25, green.f=1.25, blue.f=1.25))
    }
    
    if (length(plotlegend)>0 && plotlegend != F) {
        shapes <- 1:length(parameters)%%3 + 1
        labels <- as.expression(seis.labs[parameters])
        if (plotlegend == 'left') {
            text.widths <- rep(0.175, 5)
            text.widths[2] <- 0.2
            text.widths[3] <- 0.2
            text.widths[4] <- 0.4
            if (variable == "num_trees") text.widths <- text.widths / 1.5
            if (variable == "num_points") text.widths <- text.widths / 3.001
            
            legend2("topleft", legend=labels[1:5], horiz=T, bty='n', 
                inset=c(0.21, -0.11), cex=text.cex, text.col='white', 
                pch=c(16, 15, 18)[shapes], col=col.pal, 
                text.width=text.widths)
        
            legend2("topleft", legend=labels[1:5], horiz=T, bty='n', 
                inset=c(0.21, -0.11), cex=text.cex, 
                pch=c(1, 0, 5)[shapes], col='black', 
                text.width=text.widths)
        }
        
        if (plotlegend == 'right') {
            text.widths <- rep(0.175, length(labels)-4)
            text.widths[3] <- 0.2
            text.widths[4] <- 0.4
            if (variable != "num_trees") text.widths <- text.widths / 1.5
            if (variable == "num_points") text.widths <- text.widths / 2.9
            
            legend2("topleft", 
                legend=labels[-1:-5], horiz=T, bty='n', 
                inset=c(-0.32, -0.11), cex=text.cex, text.col='white', 
                pch=c(18, 16, 15)[shapes], col=col.pal[-1:-5], 
                text.width=text.widths)
        
            legend2("topleft", 
                legend=labels[-1:-5], horiz=T, bty='n', 
                inset=c(-0.32, -0.11), cex=text.cex, 
                pch=c(5, 1, 0)[shapes], col='black', 
                text.width=text.widths)
        }
    }
}

for (score in c("r2", "ev", "dist", "diff", "sigma")) {
    for (variable in c("num_tracks", "num_trees", "num_points")) {
        print(paste(score, variable))
        make_plots(plot_accuracy, paste0(variable, '-', score), 
            filepath=out.dir, 
            variable=variable, score=score, 
            plotlegend=F, 
            wide=F, mar=c(3, 3.75, 1, 1), 
            DF=read.table(file.path('subsetter', paste0(variable, '.dat')),
                  header=T, stringsAsFactors=F))
    }
}

