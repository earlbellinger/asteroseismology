#### Collect nearly-evenly-spaced points from all evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(lpSolve)
library(parallel)
library(parallelMap)

args <- commandArgs(TRUE)
sim_dir <- if (length(args)>0) args[1] else 'simulations'
output_fname <- paste0(sim_dir, '.dat')

simulations <- file.path(sim_dir, list.files(sim_dir))
simulations <- simulations[grep('.dat', simulations)]

# Load data
load_data <- function(filename, num_points=150, space_var='H') {
    DF <- read.table(filename, header=1, check.names=0)
    DF$age <- DF$age - min(DF$age)
    
    nrow.DF <- length(DF[[space_var]])
    if (nrow.DF < num_points) return(NULL)
    ideal <- seq(max(DF[[space_var]]), min(DF[[space_var]]), length=num_points)
    cost.mat  <- outer(ideal, DF[[space_var]], function(x, y) abs(x-y))
    row.signs <- rep("==", num_points)
    row.rhs   <- rep(1, num_points)
    col.signs <- rep("<=", nrow.DF)
    col.rhs   <- rep(1, nrow.DF)
    sol <- lp.transport(cost.mat, "min", row.signs, row.rhs,
        col.signs, col.rhs)$solution
    DF[apply(sol, 1, which.max),]
}
parallelStartMulticore(max(1, detectCores()))
seis.DF <- do.call(rbind, parallelMap(load_data, simulations))
print(paste(nrow(seis.DF), "rows"))
print(sapply(seis.DF, fivenum))

# Save data
write.table(seis.DF, output_fname, quote=FALSE, sep='\t', row.names=FALSE)

