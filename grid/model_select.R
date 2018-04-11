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
#if (grepl('LOGS_MS', log_dir)) {
#    decreasing_L <- which(diff(DF$L) < 0 & DF$center_h1[-1] > 0.55)
#    if (any(decreasing_L)) {
#        goes_back_up <- diff(decreasing_L) > 1
#        pms <- max(decreasing_L)
#        DF <- DF[-1:-pms,]
#    }
#}

## solve linear transport problem to get equally-spaced points 
DF$lg_X_c <- log10(DF$center_h1)
DF$lg_Y_c <- log10(DF$center_he3 + DF$center_he4)
## solve linear transport problem to get equally-spaced points 
space_var <- ifelse(grepl('LOGS_MS', log_dir), 'lg_X_c', 
             ifelse(grepl('LOGS_HEB', log_dir), 'lg_Y_c', 
             #ifelse(grepl('LOGS_BUMP', log_dir), 'center_degeneracy', 
                'star_age'))
#space_var <- ifelse(grepl('LOGS_MS', log_dir), 'center_h1', 'star_age')
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

if (grepl('LOGS_BUMP', log_dir)) {
    new.DF <- rbind(new.DF, DF[which.max(DF$log_L),])#[-1,])
}

model_nums <- new.DF$model_number

profiles.index <- read.table(file.path(log_dir, 'profiles.index'), 
    header=F, skip=1, col.names=c('model_num', 'priority', 'profile_num'))

cat(file.path(log_dir, paste0('profile',
    profiles.index[profiles.index$model_num %in% model_nums,]$profile_num,
    '.data.GYRE\n')))

