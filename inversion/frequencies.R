#### Helio- and astero-seismic inversions 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

#source('models.R')
library(plyr)

Sun.MDI <- list(freq.path=file.path('data', 'SolarFreq_MDI.txt'),
                freq.col.names=c('l', 'n', 'nu', 'dnu'))

MDI <- list(freq.path=file.path('data', 'Rhodes-MDI.dat'),
                freq.col.names=F)

BiSON <- list(freq.path=file.path('data', 'Sun-freqs.dat'),
              freq.col.names=F)
              
BiSON.MDI <- list(freq.path=file.path('data', 'BiSON-MDI.dat'),
              freq.col.names=F)

Tagesstern <- list(freq.path=file.path('data', 'Tagesstern-freqs.dat'),
              freq.col.names=F)

CygA <- list(freq.path=file.path('data', '16CygA-freqs.dat'),
                   freq.col.names=F)

CygB <- list(freq.path=file.path('data', '16CygB-freqs.dat'),
                   freq.col.names=F)

subgiant <- list(freq.path=file.path('data', '16CygA-freqs.dat'),
                   freq.col.names=F)

freqs.dir <- file.path('..', 'regression', 'data', 'inversions')
for (filename in list.files(freqs.dir)) {
    star.name <- paste0('KIC_', strsplit(filename, '-freqs.dat')[[1]][1])
    assign(star.name, 
        list(freq.path=file.path(freqs.dir, filename), freq.col.names=F))
}
if (F) {
KIC_8006161 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '8006161-freqs.dat'), freq.col.names=F)
KIC_6106415 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '6106415-freqs.dat'), freq.col.names=F)
KIC_12258514 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '12258514-freqs.dat'), freq.col.names=F)
KIC_6225718 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '6225718-freqs.dat'), freq.col.names=F)
KIC_10068307 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '10068307-freqs.dat'), freq.col.names=F)
KIC_6116048 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '6116048-freqs.dat'), freq.col.names=F)
KIC_3632418 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '3632418-freqs.dat'), freq.col.names=F)
KIC_7510397 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '7510397-freqs.dat'), freq.col.names=F)
KIC_8938364 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '8938364-freqs.dat'), freq.col.names=F)
KIC_5774694 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '5774694-freqs.dat'), freq.col.names=F)
KIC_8760414 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '8760414-freqs.dat'), freq.col.names=F)
KIC_7970740 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '7970740-freqs.dat'), freq.col.names=F)
KIC_5184732 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '5184732-freqs.dat'), freq.col.names=F)
KIC_7940546 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '7940546-freqs.dat'), freq.col.names=F)
KIC_10162436 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '10162436-freqs.dat'), freq.col.names=F)
KIC_8379927 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '8379927-freqs.dat'), freq.col.names=F)
KIC_8228742 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '8228742-freqs.dat'), freq.col.names=F)
KIC_8694723 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '8694723-freqs.dat'), freq.col.names=F)
#KIC_6933899 <- list(freq.path=file.path('..', 'regression', 'data', 
#    'legacy', '6933899-freqs.dat'), freq.col.names=F)
#KIC_3656476 <- list(freq.path=file.path('..', 'regression', 'data', 
#    'legacy', '3656476-freqs.dat'), freq.col.names=F)
KIC_9414417 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '9414417-freqs.dat'), freq.col.names=F)
KIC_4914923 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '4914923-freqs.dat'), freq.col.names=F)
KIC_12317678 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '12317678-freqs.dat'), freq.col.names=F)
KIC_10454113 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '10454113-freqs.dat'), freq.col.names=F)
KIC_9139163 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '9139163-freqs.dat'), freq.col.names=F)
KIC_7680114 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '7680114-freqs.dat'), freq.col.names=F)
KIC_10516096 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '10516096-freqs.dat'), freq.col.names=F)
KIC_10963065 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '10963065-freqs.dat'), freq.col.names=F)
KIC_12009504 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '12009504-freqs.dat'), freq.col.names=F)
KIC_9098294 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '9098294-freqs.dat'), freq.col.names=F)
KIC_5773345 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '5773345-freqs.dat'), freq.col.names=F)
KIC_8394589 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '8394589-freqs.dat'), freq.col.names=F)
KIC_6679371 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '6679371-freqs.dat'), freq.col.names=F)
KIC_7103006 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '7103006-freqs.dat'), freq.col.names=F)
KIC_12069424 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '12069424-freqs.dat'), freq.col.names=F)
KIC_12069449 <- list(freq.path=file.path('..', 'regression', 'data', 
    'legacy', '12069424-freqs.dat'), freq.col.names=F)
}

get_freqs <- function(target.name="BiSON", mode.set=NA, error.set=NA, 
        perturb=F, half=F) {
    if (exists(target.name)) {
        target <- get(target.name)
    } else if (exists('models') && target.name %in% names(models)) {
        target <- models[[target.name]]
    } else {
        cat(paste("Error: can't find target", target.name, '\n'))
        return(NULL)
    }
    freqs <- parse_freqs(path=target$freq.path, col.names=target$freq.col.names)
    
    if (!is.na(mode.set)) {
        mode.set <- get(mode.set)
        m.set <- with(mode.set, parse_freqs(path=freq.path, 
                                                 col.names=freq.col.names))
        m.set <- rename(m.set, c("dnu"="m.set.dnu"))
        freqs <- with(merge(freqs, m.set, by=c('l', 'n')),
             data.frame(l=l, n=n, nu=nu.x, dnu=m.set.dnu))
    }
    if (!is.na(error.set)) {
        error.set <- get(error.set)
        e.set <- with(error.set, parse_freqs(path=freq.path, 
                                             col.names=freq.col.names))
        e.set <- rename(e.set, c("dnu"="e.set.dnu"))
        freqs <- with(merge(freqs, e.set, by=c('l', 'n')),
             data.frame(l=l, n=n, nu=nu.x, dnu=e.set.dnu))
    }
    
    if (half) freqs$dnu <- freqs$dnu / 2
    if (perturb) freqs$nu <- rnorm(nrow(freqs), freqs$nu, freqs$dnu)
    
    freqs[order(freqs$l, freqs$n),]
}
