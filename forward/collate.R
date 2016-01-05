#### Collect nearly-evenly-spaced points from all evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(lpSolve)
library(parallel)
library(parallelMap)

sim_dir <- 'simulations'
simulations <- file.path(sim_dir, list.files(sim_dir))
simulations <- simulations[grep('.dat', simulations)]

# Load data
load_data <- function(f) {
    DF <- read.table(f, header=1, check.names=0)
    DF <- DF[,-which(grepl("mass|Dnu_", names(DF)))]
    DF$age <- DF$age - min(DF$age)
    DF <- DF[DF$age <= 13.8,]
    
    n <- 100
    if (nrow(DF) < n) return(NULL)
    N <- length(DF$Hc)
    ideal <- seq(max(DF$Hc), min(DF$Hc), length=n)
    cost.mat  <- outer(ideal, DF$Hc, function(x, y) abs(x-y))
    row.signs <- rep("==", n)
    row.rhs   <- rep(1, n)
    col.signs <- rep("<=", N)
    col.rhs   <- rep(1, N)
    sol <- lp.transport(cost.mat, "min", row.signs, row.rhs,
        col.signs, col.rhs)$solution
    DF[apply(sol, 1, which.max),]
}
parallelStartMulticore(max(1, detectCores()/2))
seis.DF <- do.call(rbind, parallelMap(load_data, simulations))
print(nrow(seis.DF))
print(sapply(seis.DF, fivenum))

# Save data
write.table(seis.DF, file.path('grids', 'deleter.dat'), quote=FALSE, 
    sep='\t', row.names=FALSE)
