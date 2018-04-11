#### Runs OLA_optimize_MR2 for a given star 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

set.seed(0)
options(scipen=99999)

### LIBRARIES 
args <- commandArgs(TRUE)
star <- if (length(args)>0) args[1] else '8006161'
covs.directory <- if (length(args)>1) args[2] else file.path('learn-inversions', 
    'covs-simulations', 'inversions')
obs.directory <- if (length(args)>2) args[3] else 'inversions'

DF <- read.table(file.path('..', 'regression', covs.directory, 
    paste0(star, '.dat')), header=1)

inv.file <- file.path('..', 'regression', 'data', 'inversions', 
    paste0(star, '-obs.dat'))

params <- if (file.exists(inv.file)) {
    read.table(inv.file, header=1)
} else {
    read.table(file.path('..', 'regression', 'data', obs.directory, 
        paste0(star, '-obs.dat')), header=1)
}

arr <- star
if (star == '12069424') arr <- 'CygA'
if (star == '12069449') arr <- 'CygB'

command <- paste0('maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R KIC_', 
    star, ' KIC_', star, ' KIC_', star, ' ', star, '_meanRmeanM', ' ', arr,
    ' mod_Gauss ')

M.mean <- mean(DF$M)
M.sd <- sd(DF$M)

if ('radius' %in% params$name) {
    R.mean <- params$value[params$name == 'radius']
    R.sd <- params$uncertainty[params$name == 'radius']
    command <- paste0(command, 
        M.mean, ' ', R.mean, ' ', M.sd, ' ', R.sd, ' 128 F')
} else {
    R.mean <- mean(DF$radius)
    R.sd <- sd(DF$radius)
    
    command <- paste0(command, 
        M.mean, ' ', R.mean, ' ', M.sd, ' ', R.sd, ' 128 F ',
        cov(data.frame(DF$M, DF$radius))[[2]], ' ', 
        cov(data.frame(DF$M, DF$radius))[[3]])
}
print(command) 
system(command) 

