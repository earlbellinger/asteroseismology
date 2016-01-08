#### Mesh and scatterplot analysis of evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

## Libraries
source(file.path('..', 'scripts', 'utils.R'))
source('grid-meta.R')

library(magicaxis)
library(RColorBrewer)
library(akima)
library(parallel)
library(parallelMap)
library(data.table)
library(lattice)

## Load data
seis.DF <- data.table(read.table('grid.dat', header=1))
setkey(seis.DF, M, Y, Z, alpha)
keys <- key(seis.DF)

solar_vals <- read.table(
        file.path('..', 'inverse', 'perturb', 'Sun_perturb.dat'), 
    nrow=1, header=1)

cygA_stds <- sqrt(diag(var(read.table(
        file.path('..', 'inverse', 'perturb', '16CygA_perturb.dat'), 
    header=1))))

# Sort data
combos <- unique(seis.DF[,keys, with=0])
ages <- unlist(Map(function(i) max(merge(seis.DF, combos[i,])$age), 
    1:nrow(combos)))
combos <- combos[order(ages),]

col.pal <- colorRampPalette(brewer.pal(11, "Spectral"))(1000)

# Make inputs diagram
inputs_plot <- function(text.cex, ...) {
    H <- 1-combos$Y-combos$Z
    varmax <- max(H)
    varmin <- min(H)
    
    #print(ggpairs(data=combos, #font=font, #axisLabels="show", 
    #        columnLabels=sapply(names(seis.DF)[1:4], 
    #                            function (name) get_label(name))))
    
    print(splom(combos,
        cex=0.5, pch=20,
        col=col.pal[floor((H-varmin)/(varmax-varmin)*(length(col.pal)-1))+1],
        xlab=NULL, ylab=NULL, 
        axis.text.cex=3*text.cex/4, 
        axis.text.lineheight=0.0001,
        axis.line.tck=0.25,
        axis.half = 0,
        axis.check.overlap = 1,
        xaxs='n', yaxs='n',
        varname.cex=3*text.cex/4, 
        varnames=as.expression(seis.labs[1:4])))
        #sapply(names(seis.DF)[1:4], function (name) get_label(name))))
}
make_plots(inputs_plot, "inputs", filepath=file.path("plots", "mesh"))

## Scatter plot definitions 
scatter_plot <- function(...) { 
    Z_max <- max(seis.DF[[Z]])
    Z_min <- min(seis.DF[[Z]])
    for (simulation_i in 1:nrow(combos)) {
        track <- merge(seis.DF, combos[simulation_i,])
        color <- col.pal[floor((track[[Z]]-Z_min)/(Z_max-Z_min)*
                                   (length(col.pal)-1))+1]
        use_line <- length(unique(color))==1
        relation <- track[[Y]] ~ track[[X]]
        if (simulation_i == 1) {
            plot(relation, type=ifelse(use_line, 'l', 'p'),
                 pch=20, axes=FALSE, col=color, cex=0.01, tcl=0, 
                 ylim=ylim, xlim=xlim,
                 xlab=xlab, ylab=ylab)
            if (has_solar) {
                abline(v=solar_x, lty=3, col='black')
                abline(h=solar_y, lty=3, col='black')
            }
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.25, 0))
        } else {
            if (use_line) lines(relation, col=color[1])
            else points(relation, col=color, pch=20, cex=0.01)
        }
    }
    
    points(solar_x, solar_y, pch=1, cex=1)
    points(solar_x, solar_y, pch=20, cex=0.1)
    X_range <- diff(par()$usr)[1]
    color.legend(par()$usr[2]+0.05*X_range, par()$usr[3], 
                 par()$usr[2]+0.10*X_range, par()$usr[4], 
                 signif(quantile(seq(Z_min, Z_max, length=1000), 
                                 c(0, 0.25, 0.5, 0.75, 1)), 2), 
                 col.pal[1:length(col.pal)], gradient='y', align='rb')
    mtext(get_label(Z), 4, line=4.5, cex=1.3)
}

## Mesh plotting definitions 
mesh_plot <- function(text.cex, ...) {
    filled.contour(mesh,
                   xlim=if (X=='Teff') rev(range(mesh$x)) else range(mesh$x),
                   levels=color_levels[[Z]], 
                   color=colorRampPalette(brewer.pal(11, "Spectral")),
                   key.axes={
                       axis(4, cex.axis=text.cex, tcl=0, line=0)
                       mtext(get_label(Z), side=4, las=3, line=3, cex=text.cex)
                   },
                   plot.axes={
                       contour(mesh, add=TRUE, labcex=0.5, levels=Z_levels[[Z]])
                       if (has_solar) {
                           points(solar_x, solar_y, pch=1, cex=1)
                           points(solar_x, solar_y, pch=20, cex=0.1)
                       }
                       magaxis(side=1:4, family=font, tcl=0.25, 
                               labels=c(1,1,0,0), mgp=utils.mgp, 
                               cex.axis=text.cex)
                   },
                   plot.title={
                       title(xlab=xlab, cex.lab=text.cex, line=3)
                       title(ylab=ylab, cex.lab=text.cex, line=3)
                   })
}

for (Z in names(seis.DF)[1:9]) {
    xy_names <- names(seis.DF)[-1:-9]
    for (X_i in 1:length(xy_names)) {
        X <- xy_names[X_i]
        for (Y_i in (1+X_i):length(xy_names)) {
            Y <- xy_names[Y_i]
            if (Y == 'Teff' && X == 'L') {
                Y = 'L'
                X = 'Teff'
            }
            xlab <- get_label(X)
            ylab <- get_label(Y)
            ylim <- quantile(seis.DF[[Y]], c(0.0001, 0.9999))
            xlim <- quantile(seis.DF[[X]], c(0.0001, 0.9999))
            if (X == 'Teff') xlim <- rev(xlim)
            has_solar <- any(grepl(X, names(solar_vals))) && 
                any(grepl(Y, names(solar_vals)))
            if (has_solar) {
                solar_x <- solar_vals[[X]]
                solar_y <- solar_vals[[Y]]
            }
            
            make_plots(scatter_plot, 
                       paste0(paste(Z, Y, X, sep="_"), "-scatter"),
                       filepath=file.path("plots", "mesh"),
                       mar=c(3,4,1,6), thin=FALSE)
            
            xx <- normalize(seis.DF[[X]])
            yy <- normalize(seis.DF[[Y]])
            
            cygx <- cygA_stds[[X]]/(max(seis.DF[[X]])-min(seis.DF[[X]]))
            cygy <- cygA_stds[[Y]]/(max(seis.DF[[Y]])-min(seis.DF[[Y]]))
            
            mesh <- interp(xx, yy, seis.DF[[Z]],
               xo=seq(0, 1, cygx), yo=seq(0, 1, cygy))
            
            mesh$x <- seq(min(seis.DF[[X]]), max(seis.DF[[X]]), cygA_stds[[X]])
            mesh$y <- seq(min(seis.DF[[Y]]), max(seis.DF[[Y]]), cygA_stds[[Y]])
            
            make_plots(mesh_plot, 
                       paste0(paste(Z, Y, X, sep="_"), "-mesh"),
                       filepath=file.path("plots", "mesh"),
                       mar=c(5, 6, 1, 0), thin=FALSE)
        }
    }
}

