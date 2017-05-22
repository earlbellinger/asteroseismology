#### Collect evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(parallel)
library(parallelMap)

args <- commandArgs(TRUE)
sim_dir <- if (length(args)>0) args[1] else 'simulations'
output_fname <- paste0(sim_dir, '.dat')

simulations <- file.path(sim_dir, list.files(sim_dir))
simulations <- simulations[grep('.dat', simulations)]

load_data <- function(filename) {
    DF <- read.table(filename, header=1, check.names=0)
    if (any(DF['Fe/H'] < -7.5)) DF[DF['Fe/H'] < -7.5,]['Fe/H'] <- -7.5
    DF[DF < 10e-20 & DF > -10e-20 & !is.na(DF)] <- 0
    return(DF)
}

# Load and combine data over all simulations 
parallelStartMulticore(max(1, min(detectCores(), 62)))
seis.DF <- do.call(plyr:::rbind.fill, parallelMap(load_data, simulations))
#seis.DF <- seis.DF[complete.cases(seis.DF),]
print(paste(nrow(seis.DF), "rows"))
print(sapply(seis.DF, fivenum))

# Save data
write.table(seis.DF, output_fname, quote=FALSE, sep='\t', row.names=FALSE)

