#### Determine which profile files should have pulsation frequencies calculated 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(lpSolve)

args <- commandArgs(TRUE)
num_points <- if (length(args)>0) as.numeric(args[1]) else 64
log_dir <- if (length(args)>1) args[2] else 'LOGS'

DF <- read.table(file.path(log_dir, 'history.data'), header=1, skip=5)[-1,]

## if main sequence, crop PMS
if (grepl('LOGS_MS', log_dir)) {
    decreasing_L <- which(diff(DF$L) < 0 & DF$center_h1[-1] > 0.55)
    if (any(decreasing_L)) {
        goes_back_up <- diff(decreasing_L) > 1
        pms <- max(decreasing_L)
        DF <- DF[-1:-pms,]
    }
}

DF <- DF[with(DF, 10**log_R>0.94, 
                  10**log_R<1.06, 
                  10**log_L>0.85, 
                  10**log_L<1.15),]#$model_number

model_nums <- if (nrow(DF) < num_points) {
    DF$model_number
} else {
    ## solve linear transport problem to get equally-spaced points 
    space_var <- ifelse(grepl('LOGS_MS', log_dir), 'center_h1', 'star_age')
    x <- DF[[space_var]]
    nrow.DF <- length(x)
    ideal <- seq(max(x), min(x), length=num_points)
    cost.mat  <- outer(ideal, x, function(x, y) abs(x-y))
    row.signs <- rep("==", num_points)
    row.rhs   <- rep(1, num_points)
    col.signs <- rep("<=", nrow.DF)
    col.rhs   <- rep(1, nrow.DF)
    sol <- lp.transport(cost.mat, "min", row.signs, row.rhs,
        col.signs, col.rhs)$solution
    new.DF <- DF[apply(sol, 1, which.max),]

    new.DF$model_number
}

profiles.index <- read.table(file.path(log_dir, 'profiles.index'), 
    header=F, skip=1, col.names=c('model_num', 'priority', 'profile_num'))

cat(file.path(log_dir, paste0('profile',
    profiles.index[profiles.index$model_num %in% model_nums,]$profile_num,
    '.data.FGONG\n')))

