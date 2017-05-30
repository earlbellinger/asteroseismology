#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

#source('models.R')
library(plyr)

Sun.MDI <- list(freq.path=file.path('data', 'SolarFreq_MDI.txt'),
                freq.col.names=c('l', 'n', 'nu', 'dnu'))

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

get_freqs <- function(target.name="BiSON", mode.set=NA, error.set=NA, 
        perturb=F, half=F) {
    target <- get(target.name)
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
