#### Histogram of sub-giant masses 
#### Also plots the diffusion factor as a function of mass for KAGES stars 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(parallel)
library(parallelMap)
parallelStartMulticore(max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))

source(file.path('..', 'scripts', 'utils.R'))

#hare.path <- file.path('learn-giants', 'covs-simulations(1|2|3)', 'sg-hares-mesa')
hare.path <- file.path('learn-giants', 'covs-simulations0', 'sg-hares-mesa')
hare.fnames <- list.files(hare.path)
#hare.fnames <- hare.fnames[grepl('Model', hare.fnames)]

#hares <- read.table(file.path('data', 'sg-hares', 'ana_16_2', 'track.dat'), 
#    header=1)
hares <- read.table(file.path('data', 'sg-hares-mesa.dat'), header=1)
names(hares)[1] <- 'model'

get_masses <- function(directory) {
    as.numeric(parallelMap(function(fname) 
        mean(read.table(file.path(directory, fname), header=1)$M),
    list.files(directory)))
}

#masses.x <- get_masses(file.path('learn-giants', 'covs-simulationsx', 'sg-hares-mesa'))
#masses.1 <- get_masses(file.path('learn-giants', 'covs-simulations1', 'sg-hares-mesa'))
#masses.2 <- get_masses(file.path('learn-giants', 'covs-simulations(1|2)', 'sg-hares-mesa'))
#masses.3 <- get_masses(file.path('learn-giants', 'covs-simulations(1|2|3)', 'sg-hares-mesa'))

masses.x <- get_masses(file.path('learn-giants', 'covs-simulations1', 'sg-hares-mesa'))
masses.1 <- get_masses(file.path('learn-giants', 'covs-simulations0', 'sg-hares-mesa'))

plot_hist <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font) {
    par(cex.lab=text.cex,cex.axis=text.cex,cex.main=text.cex,cex.sub=text.cex,
        mgp=mgp+c(0, 0.4, 0))
    hist(masses, xlab="Predicted Hare Mass", ylab="Counts", main='', las=1)
}

make_plots(plot_hist, 'hare-sg-mass-hist', 
    filepath=file.path('plots', 'hare-sg'))


plot_density <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font) {
    #par(cex.lab=text.cex,cex.axis=text.cex,cex.main=text.cex,cex.sub=text.cex,
    #    mgp=mgp+c(0, 0.4, 0))
    #hist(masses, xlab="Predicted Hare Mass", ylab="Counts", main='', las=1)
    dens <- density(masses.x)
    plot(dens$x, dens$y/max(dens$y), axes=F, tcl=0, type='l', 
        xlab=expression("Predicted Hare Mass"~M/M["Sun"]),
        ylab="Normalized Density",
        main="",
        xlim=c(0.7, 3),
        ylim=c(0, 1))
    abline(v=1.6, lty=2)
    abline(h=0, lty=1, col="#00000066")
    with(density(masses.1), lines(x, y/max(y), col=2))
    #with(density(masses.2), lines(x, y/max(y), col=3))
    #with(density(masses.3), lines(x, y/max(y), col=4))
	legend("topright", col=1:4, lty=c(1,1,2),
	    legend=c('Without mixed modes', 
		         'With closest mixed mode', 
				 #'Closest two',
				 #'Closest three',
				 'True mass'))
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp+c(0,-0.1,0), family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1,
            las=1, mgp=mgp+c(0,0.1,0), family=font, cex.axis=text.cex)
}

make_plots(plot_density, 'hare-sg-mass-pdf-new', 
    filepath=file.path('plots', 'hare-sg'))

get_radii <- function(directory) {
    as.numeric(parallelMap(function(fname) 
        mean(read.table(file.path(directory, fname), header=1)$radius),
    list.files(directory)))
}

radii.x <- get_radii(file.path('learn-giants', 'covs-simulations1', 'sg-hares-mesa'))
radii.1 <- get_radii(file.path('learn-giants', 'covs-simulations0', 'sg-hares-mesa'))

#radii.x <- get_radii(file.path('learn-giants', 'covs-simulationsx', 'sg-hares-mesa'))
#radii.1 <- get_radii(file.path('learn-giants', 'covs-simulations1', 'sg-hares-mesa'))
#radii.2 <- get_radii(file.path('learn-giants', 'covs-simulations(1|2)', 'sg-hares-mesa'))
#radii.3 <- get_radii(file.path('learn-giants', 'covs-simulations(1|2|3)', 'sg-hares-mesa'))

#radii.x <- radii.x[1:104]
#radii.1 <- radii.1[1:104]
#radii.2 <- radii.2[1:104]
#radii.3 <- radii.3[1:104]
#radii <- hares$radius[1:104]



get_radii <- function(directory) {
    do.call(rbind, parallelMap(function(fname) {
        DF <- read.table(file.path(directory, fname), header=1)
        data.frame(model=as.numeric(sub('.dat', '', sub('Model', '', fname))),
                   radius=median(DF$radius),
                   dRL=quantile(DF$radius, 0.16),
                   dRH=quantile(DF$radius, 0.84))
    }, list.files(directory)))
}
radii.1 <- get_radii(file.path('learn-giants', 'covs-simulations1', 'sg-hares-mesa'))
radii.x <- get_radii(file.path('learn-giants', 'covs-simulations0', 'sg-hares-mesa'))

#radii.3 <- radii.3[-nrow(radii.3),]
#radii.x <- radii.x[-nrow(radii.x),]

plot_radii <- function(..., text.cex=1, mgp=utils.mgp, font='Palatino') {
    plot(NA, axes=F, 
        ylab=expression("ML Radius"~R/R["Sun"]), 
        xlab=expression("Hare Radius"~R/R["Sun"]), 
        xlim=range(hares$radius), 
        ylim=range(radii.1[,-1], radii.x[,-1]))
    segments(x0=-5,y0=-5,x1=45,y1=45, lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    
    comb <- merge(hares, radii.1, by='model')
    with(comb, arrows(radius.x, dRL, 
                      radius.x, dRH, 
        length=0.01, lwd=1.5, angle=90, code=3, col="#ca002066"))
    with(comb, points(radius.y ~ radius.x, pch=20, cex=0.5))
    
    comb <- merge(hares, radii.x, by='model')
    with(comb, arrows(radius.x, dRL, 
                      radius.x, dRH, 
        length=0.01, lwd=1.5, angle=90, code=3, col="#0571b066"))
    with(comb, points(radius.y ~ radius.x, pch=4, cex=0.5))
    
    legend('bottomright', col=c(red, blue), lty=1, pch=c(20, 4), cex=text.cex,
        legend=c("No mixed modes", "One mixed mode"))
}
make_plots(plot_radii, 'hare-radii-new', filepath=file.path('plots', 'hare-sg'))



comb <- merge(hares, radii.x, by='model')
make_plots(plot_radii, 'hare-radii-x', filepath=file.path('plots', 'hare-sg'))

plot_density <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font) {
    dens <- density(radii.x)
    plot(dens$x, dens$y/max(dens$y), axes=F, tcl=0, type='l', 
        xlab="Predicted Hare Radius",
        ylab="Normalized Density",
        main="",
        xlim=c(2, 4),
        ylim=c(0, 1))
    #abline(v=1.6, lty=2)
    abline(h=0, lty=1, col="#00000066")
    with(density(radii.1), lines(x, y/max(y), col=2))
    with(density(radii.2), lines(x, y/max(y), col=3))
    with(density(radii.3), lines(x, y/max(y), col=4))
    with(density(hares$radius[-1]), lines(x, y/max(y), lty=2))
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp+c(0,-0.1,0), family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1,
            las=1, mgp=mgp+c(0,0.1,0), family=font, cex.axis=text.cex)
}

make_plots(plot_density, 'hare-sg-radii-pdf', 
    filepath=file.path('plots', 'hare-sg'))

#logL <- do.call(rbind, parallelMap(function(fname) {
#    DF <- read.table(file.path(hare.path, fname), header=1)
#    data.frame(logL=log10(mean(DF$L)), 
#              dlogL=sqrt(var(DF$L))/(mean(DF$L)*log(10)))
#    }, hare.fnames))



ages <- do.call(rbind, parallelMap(function(fname) {
    DF <- read.table(file.path(hare.path, fname), header=1)
    data.frame(model=as.numeric(sub('.dat', '', fname)),
               age=median(DF$age), 
              dageL=quantile(DF$age, .16),
              dageH=quantile(DF$age, .84))
    }, hare.fnames))
comb <- merge(hares, ages, by='model')



get_ages <- function(directory) do.call(rbind, parallelMap(function(fname) {
    DF <- read.table(file.path(directory, fname), header=1)
    data.frame(model=as.numeric(sub('.dat', '', fname)),
               age=median(DF$age), 
              dageL=quantile(DF$age, .16),
              dageH=quantile(DF$age, .84))
    }, list.files(directory)))
ages.1 <- get_ages(file.path('learn-giants', 'covs-simulations0', 'sg-hares-mesa'))
ages.x <- get_ages(file.path('learn-giants', 'covs-simulations1', 'sg-hares-mesa'))

#ages.3 <- ages.3[-nrow(ages.3),]
#ages.x <- ages.x[-nrow(ages.x),]



get_ages <- function(directory) do.call(rbind, parallelMap(function(fname) {
    DF <- read.table(file.path(directory, fname), header=1)
    data.frame(model=as.numeric(sub('.dat', '', fname)),
               age=median(DF$age), 
              dageL=quantile(DF$age, .16),
              dageH=quantile(DF$age, .84))
    }, list.files(directory)))
ages.3 <- get_ages(file.path('learn-giants', 'covs-simulations(1|2|3)', 'sg-hares-mesa'))
ages.x <- get_ages(file.path('learn-giants', 'covs-simulationsx', 'sg-hares-mesa'))


plot_ages <- function(..., text.cex=1, mgp=utils.mgp, font='Palatino') {    
    plot(NA, axes=F, 
        ylab=expression("ML Age"~tau/"Gyr"), 
        xlab=expression("Hare Age"~tau/"Gyr"), 
        xlim=with(merge(hares, ages.x, by='model'), range(age.x)), 
        ylim=range(ages.1[,-1], ages.x[,-1]))
    segments(x0=-5,y0=-5,x1=45,y1=45, lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    with(merge(hares, ages.1, by='model'), 
	    arrows(age.x, dageL, age.x, dageH, 
        length=0.01, lwd=1.5, angle=90, code=3, col="#ca002066"))
    with(merge(hares, ages.x, by='model'), 
	    arrows(age.x, dageL, age.x, dageH, 
        length=0.01, lwd=1.5, angle=90, code=3, col="#0571b066"))
    with(merge(hares, ages.1, by='model'), points(age.y ~ age.x, pch=20, cex=0.5))
	with(merge(hares, ages.x, by='model'), points(age.y ~ age.x, pch=20, cex=0.5))
	legend('topleft', lty=1, cex=text.cex, col=c(red, blue),
	    legend=c("One Mixed Mode", "No Mixed Modes"))
}

make_plots(plot_ages, 'hare-ages-comb-mesa-new', filepath=file.path('plots', 'hare-sg'))


#comb <- merge(hares, ages, by='model')

plot_ages <- function(..., text.cex=1, mgp=utils.mgp, font='Palatino') {    
    plot(NA, axes=F, 
        ylab=expression("ML Age"~tau/"Gyr"), 
        xlab=expression("Hare Age"~tau/"Gyr"), 
        xlim=range(comb$age.x[-1]), 
        ylim=range(comb$dageL[-1], comb$dageH[-1]))
    segments(x0=-5,y0=-5,x1=45,y1=45, lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    with(comb, arrows(age.x, dageL, 
                      age.x, dageH, 
        length=0.01, lwd=1.5, angle=90, code=3, col=red))
    with(comb, points(age.y ~ age.x, pch=20, cex=0.5))
}

make_plots(plot_ages, 'hare-ages', filepath=file.path('plots', 'hare-sg'))

plot_age_density <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font) {
    dens <- density(comb$age.x[-1])
    plot(dens$x, dens$y/max(dens$y), axes=F, tcl=0, type='l', 
        xlab="Hare Age",
        ylab="Normalized Density",
        main="",
        xlim=c(1, 3),
        ylim=c(0, 1))
    abline(h=0, lty=1, col="#00000066")
    dens <- density(comb$age.y[-1])
    lines(dens$x, dens$y/max(dens$y), lty=2, col='darkred')
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp+c(0,-0.1,0), family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1,
            las=1, mgp=mgp+c(0,0.1,0), family=font, cex.axis=text.cex)
}

make_plots(plot_age_density, 'hare-sg-age-pdf', 
    filepath=file.path('plots', 'hare-sg'))



