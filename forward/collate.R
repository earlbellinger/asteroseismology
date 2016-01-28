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
load_data <- function(filename, num_points=70, space_var='Hc') {
    DF <- read.table(filename, header=1, check.names=0)
    
    hs <- 1-DF$Y[1]-DF$Z[1]
    dteff <- diff(DF$Teff)
    dl <- diff(DF$L)
    pms <- with(DF, 
        #which(100*(hs-H[-1])/(hs) < 1.5
        which(DF$age[-1] < 0.1
            & dteff %in% boxplot.stats(dteff, coef=3)$out
            & dl %in% boxplot.stats(dl, coef=3)$out
        )
    )
    if (any(pms)) {
        print(paste(filename, "Clipping", 1+max(pms), "points"))
        DF <- DF[-1:-(1+max(pms)),]
    } else {
        print(paste(filename, "No PMS"))
    }
    DF$age <- DF$age - min(DF$age) # set ZAMS age
    DF <- DF[DF$age <= 14,]
    
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
parallelStartMulticore(max(1, detectCores()/2))
seis.DF <- do.call(rbind, parallelMap(load_data, simulations))
seis.DF <- seis.DF[complete.cases(seis.DF),]
print(paste(nrow(seis.DF), "rows"))
print(sapply(seis.DF, fivenum))

# Save data
write.table(seis.DF, output_fname, quote=FALSE, sep='\t', row.names=FALSE)

