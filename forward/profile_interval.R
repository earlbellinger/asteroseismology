#### Determine how often a profile should be written. Also returns max year. 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

DF <- read.table('LOGS/history.data', skip=5, header=1)
decreasing_L <- which(diff(DF$log_L) < 0 & DF$star_age[-1] < 0.25*10**9)
if (any(decreasing_L)) {
    goes_back_up <- diff(decreasing_L) > 1
    pms <- ifelse(any(goes_back_up), 
                  which(goes_back_up)[1] + 1, 
                  max(decreasing_L))
    DF <- DF[-1:-pms,]
}
DF <- DF[-1:-floor(nrow(DF)*0.1),] # in case more PMS leaks in

# remove radiative regions 
cutoff <- 0
mx1_diff <- abs(DF$conv_mx1_top - DF$conv_mx1_bot)
mx2_diff <- abs(DF$conv_mx2_top - DF$conv_mx2_bot)
radiative <- with(DF, 
    mass_conv_core <= 0 & mx1_diff <= cutoff & mx2_diff <= cutoff |
    mass_conv_core  > 0 & mx2_diff <= cutoff)
if (any(radiative)) DF <- DF[1:min(which(radiative)),]

cat(paste(nrow(DF), max(DF$star_age), '\n'))

