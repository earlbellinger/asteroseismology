#### Collect nearly-evenly-spaced points from all evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(lpSolve)
library(parallel)
library(parallelMap)
library(magicaxis)
library(mblm)

args <- commandArgs(TRUE)
sim_dir <- if (length(args)>0) args[1] else 'simulations'
output_fname <- paste0(sim_dir, '.dat')

simulations <- file.path(sim_dir, list.files(sim_dir))
simulations <- simulations[grep('.dat', simulations)]

# Load data
load_data <- function(filename, num_points=32, space_var='X_c') {
    DF <- read.table(filename, header=1, check.names=0)
    #DF[DF['Fe/H'] < -10,]['Fe/H'] <- -10
    #if (any(DF['Fe/H'] < -10)) DF[DF['Fe/H'] < -10,]['Fe/H'] <- -10
    
    # clip PMS
    decreasing_L <- which(diff(DF$L) < 0 & DF$age[-1] < 0.25)
    if (any(decreasing_L)) {
        goes_back_up <- diff(decreasing_L) > 1
        pms <- ifelse(any(goes_back_up), 
                   which(goes_back_up)[1] + 1, 
                   max(decreasing_L))
        print(paste(filename, "Clipping", pms, "points"))
        DF <- DF[-1:-pms,]
    }
    
    # detect outliers and remove them
    while (T) {
       dd1 <- diff(DF$Dnu0_median)
       dd2 <- diff(DF$dnu02_median)
       outliers <- dd1 > 30 | dd2 > 10
       if (any(outliers)) {
           print(paste(filename, "Rejecting", sum(outliers), "outliers"))
           DF <- DF[-(1+which(outliers)),]
       } else break
    }
    #for (name in names(DF)[grep("median", names(DF))]) {
    #    while (T) {
    #       #resids <- resid(lm(DF[[name]] ~ DF$age + I(DF$age^2) + I(DF$age^3)))
    #        #ages <- DF$age
    #        #other <- DF[[name]]
    #        dd <- diff(DF[[name]])
    #        outliers <- dd > median(dd) + 100*mad(dd)
    #        #outliers <- abs(resid(mblm(other~ages))) > 10
    #        #outliers <- abs(resids) > 10
    #        if (any(outliers)) {
    #            print(paste("Rejecting", sum(outliers), "outliers"))
    #            DF <- DF[-(1+which(outliers)),]
    #        } else break
    #    }
    #}
    
    # set ZAMS age 
    DF$age <- DF$age - min(DF$age)
    DF <- DF[DF$age <= 15,]
    
    # solve linear transport problem to get equally-spaced points 
    x <- DF[[space_var]]
    nrow.DF <- length(x)
    if (nrow.DF < num_points) {
        print(paste(filename, "has too few points"))
        return(NULL)
    }
    ideal <- seq(max(x), min(x), length=num_points)
    cost.mat  <- outer(ideal, x, function(x, y) abs(x-y))
    row.signs <- rep("==", num_points)
    row.rhs   <- rep(1, num_points)
    col.signs <- rep("<=", nrow.DF)
    col.rhs   <- rep(1, nrow.DF)
    sol <- lp.transport(cost.mat, "min", row.signs, row.rhs,
        col.signs, col.rhs)$solution
    new.DF <- DF[apply(sol, 1, which.max),]
    
    y <- new.DF[[space_var]]
    
    make_plots(plot_spacing, paste0(basename(filename), "-spacing"), 
        filepath=file.path('plots', sim_dir, "spacing"), 
        x=x, y=y, ideal=ideal, num_points=num_points,
        paper_pdf_height=3.25, slides_pdf_height=3.25,
        mar=c(3, 1, 1, 1), thin=F, tall=F, slides=F)
    
    return(new.DF)
}

plot_spacing <- function(x, y, ideal, num_points, ...,
        text.cex=1, mgp=utils.mgp, font="Palatino") {
    plot(NA, ylim=c(0, 1), axes=F, ann=F, 
        xlim=c(0, round(max(ideal) + 0.05, 1)), 
        xaxs='i', yaxs='i')
    #axis(1, tck=-0.05, text.cex)
    magaxis(1, tck=-0.05, ratio=0.1, family=font, cex.axis=text.cex, usepar=1)
    title(xlab=expression("Core-hydrogen abundance"~X[c]))
    points(x, rep(0.35, length(x)), pch=4, cex=0.25, col='darkblue', xpd=NA)
    points(y, rep(0.225, length(y)), pch=3, cex=0.25, col='darkred', xpd=NA)
    points(ideal, rep(0.1, num_points), pch=20, cex=0.25, xpd=NA)
    legend("top", text.width=c(0.15, 0.15, 0.2), 
        legend=c("All models", 
                 "Models selected", 
                 "Equidistant spacing"), 
       pch=c(4, 3, 20), col=c('darkblue', 'darkred', 'black'), horiz=1,
       cex=text.cex)
}

# Load and combine data over all simulations 
parallelStartMulticore(max(1, min(detectCores(), 62)))
seis.DF <- do.call(rbind, parallelMap(load_data, simulations))
seis.DF <- seis.DF[complete.cases(seis.DF),]
print(paste(nrow(seis.DF), "rows"))
print(sapply(seis.DF, fivenum))

# Save data
write.table(seis.DF, output_fname, quote=FALSE, sep='\t', row.names=FALSE)

