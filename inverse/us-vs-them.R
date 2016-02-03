#### Plot KAGES masses and ages against what we get with machine learning 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)

kages <- read.table(file.path('data', 'kages.dat'), header=1)

data_dir <- file.path('learn_covs', 'kages')
ml <- do.call(plyr:::rbind.fill, Map(function(cov) {
        name <- sub('.dat', '', basename(cov))
        if (!name %in% kages$KIC) 
            return(NULL)
        DF <- read.table(cov, header=1)
        data.frame(Name = as.integer(name),
                   Age = median(DF$age),
                   dAgeL = median(DF$age) - quantile(DF$age, .16),
                   dAgeH = quantile(DF$age, .84) - median(DF$age),
                   Mass = median(DF$M),
                   dMassL = median(DF$M) - quantile(DF$M, .16),
                   dMassH = quantile(DF$M, .84) - median(DF$M)
                   #dAge = sqrt(var(DF$age))
                  )
    }, file.path(data_dir, list.files(data_dir)) ))

kages <- kages[kages$KIC %in% ml$Name,]
kages <- kages[order(kages$KIC),]
ml <- ml[order(ml$Name),]

plot_ages <- function(..., text.cex=1, mgp=utils.mgp, mar=utils.mar) {
    lims <- range(1, kages$Age-kages$dAgeL, kages$Age+kages$dAgeH,
                  ml$Age-ml$dAge, ml$Age+ml$dAge)
    plot(NA, axes=F,
        ylab=expression("Age from KAGES"~tau/"Gyr"),
        xlab=expression("Age from Machine Learning"~tau/"Gyr"),
        xlim=lims, ylim=lims)
    abline(coef=c(0,1), lty=2)
    magaxis(side=1:4, family="Palatino", tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    #with(kages, arrows(Age-dAgeL, ml$Age, Age+dAgeH, ml$Age, 
    #    length=0, angle=90, code=3, col="darkgray"))
    #with(ml, arrows(kages$Age, Age-dAgeL, kages$Age, Age+dAgeH, 
    #    length=0, angle=90, code=3, col="darkgray"))
    with(ml, arrows(Age-dAgeL, kages$Age, Age+dAgeH, kages$Age, 
        length=0, angle=90, code=3, col="darkgray"))
    with(kages, arrows(ml$Age, Age-dAgeL, ml$Age, Age+dAgeH, 
        length=0, angle=90, code=3, col="darkgray"))
    points(kages$Age ~ ml$Age, pch=1)
}

plot_masses <- function(..., text.cex=1, mgp=utils.mgp, mar=utils.mar) {
    lims <- range(0.7, kages$Mass-kages$dMassL, kages$Mass+kages$dMassH,
                  ml$Mass-ml$dMass, ml$Mass+ml$dMass)
    plot(NA, axes=F,
        ylab=expression("Mass from KAGES"~M/M["\u0298"]),
        xlab=expression("Mass from Machine Learning"~
            M/M["\u0298"]),
        xlim=lims, ylim=lims)
    abline(coef=c(0,1), lty=2)
    magaxis(side=1:4, family="Palatino", tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    #with(kages, arrows(Mass-dMassL, ml$Mass, Mass+dMassH, ml$Mass, 
    #    length=0, angle=90, code=3, col="darkgray"))
    #with(ml, arrows(kages$Mass, Mass-dMassL, kages$Mass, Mass+dMassH, 
    #    length=0, angle=90, code=3, col="darkgray"))
    with(ml, arrows(Mass-dMassL, kages$Mass, Mass+dMassH, kages$Mass, 
        length=0, angle=90, code=3, col="darkgray"))
    with(kages, arrows(ml$Mass, Mass-dMassL, ml$Mass, Mass+dMassH, 
        length=0, angle=90, code=3, col="darkgray"))
    points(kages$Mass ~ ml$Mass, pch=1)
}

make_plots(plot_ages, 'kages', filepath=file.path('plots', 'comparison'))
make_plots(plot_masses, 'kmasses', filepath=file.path('plots', 'comparison'))

