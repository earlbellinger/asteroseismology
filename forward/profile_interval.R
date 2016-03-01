#### Determine how often a profile should be written. Also returns max year. 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

DF <- read.table(file.path('LOGS', 'history.data'), header=1, skip=5)

# remove radiative regions 
cutoff <- 10**-3
mx1_diff <- abs(DF$mx1_top - DF$mx1_bot)
mx2_diff <- abs(DF$mx2_top - DF$mx2_bot)
radiative <- with(DF, 
    mass_conv_core <= 0 & mx1_diff <= cutoff & mx2_diff <= 0 |
    mass_conv_core  > 0 & mx2_diff <= 0)
if (any(radiative)) DF <- DF[1:min(which(radiative)),]

cat(paste(nrow(DF), max(DF$star_age), '\n'))

