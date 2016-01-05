#### Collect nearly-evenly-spaced points from all evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(lpSolve)
library(parallel)
library(parallelMap)

output_fname <- 'grid.dat'
sim_dir <- 'simulations'
simulations <- file.path(sim_dir, list.files(sim_dir))
simulations <- simulations[grep('.dat', simulations)]

# Load data
load_data <- function(filename, num_points=100) {
    DF <- read.table(filename, header=1, check.names=0)
    DF <- DF[, -which(grepl("mass|Dnu_", names(DF)))]
    DF$age <- DF$age - min(DF$age)
    DF <- DF[DF$age <= 13.8,]
    
    nrow.DF <- length(DF$Hc)
    if (nrow.DF < num_points) return(NULL)
    ideal <- seq(max(DF$Hc), min(DF$Hc), length=num_points)
    cost.mat  <- outer(ideal, DF$Hc, function(x, y) abs(x-y))
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
