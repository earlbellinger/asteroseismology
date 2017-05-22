#### Model all frequencies of the Sun 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

fgongs.dir <- 'high_ell'

freqs.mdi <- read.table(
    file.path('..', '..', 'inverse', 'data', 'SolarFreq_MDI.txt'),
    col.names=c('l', 'n', 'nu', 'dnu'))

freqs.bison <- read.table(
    file.path('..', '..', 'inverse', 'data', 'Sun-freqs.dat'),
    header=1)

freqs <- rbind(freqs.bison, freqs.mdi)
freqs <- freqs[freqs$l <= 200,]

ells <- sort(unique(freqs$l))
ells. <- seq(min(ells), max(ells), 10)

#ker.type <- 1

for (FGONG in c("diffusion")) { #, "no_diffusion")) { 
    for (ii in 1:(length(ells.)-1)) { 
        for (ker.type in 2:6) { #c(1)) { #
            ell.min <- ells.[ii] 
            ell.max <- ells.[ii+1] 
            freqs. <- freqs[freqs$l >= ell.min & freqs$l <= ell.max,] 
            nu.min <- max(min(freqs.$nu - 100), 1000) 
            nu.max <- min(max(freqs.$nu + 100), 4000) 
            n.min <- min(freqs.$n) 
            n.max <- max(freqs.$n) 
            cmd <- paste("cd", file.path('high_ell', FGONG), ';', 
                "maybe_sub.sh -p 1 kerexact.sh", 
                paste0(FGONG, '.FGONG'), 
                paste0(ell.min,'_',ell.max), 
                ker.type, 
                paste0(10,',',ell.min), 
                paste0(nu.min,',',nu.max), 
                paste0(n.min,',',n.max)) 
            print(cmd) 
            system(cmd) 
        }
    }
}

