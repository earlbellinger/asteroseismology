#### Get M, R, rho, P, and T from FGONG 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### PARSE COMMAND LINE ARGS
args <- commandArgs(TRUE)
directory <- if (length(args)>0) args[1] else "G2V"
print(directory)

hdr <- read.table(
        file.path('models', directory, 'LOGS_MS', 'profile1.data'), 
    header=1, nrow=1, skip=1)

prof <- read.table(
        file.path('models', directory, 'LOGS_MS', 'profile1.data.FGONG.dat'), 
    header=1)

output <- with(prof, data.frame(r, x, rho, P, T))
colnames(output)[1] <- 'r/cm'
colnames(output)[2] <- 'r/R_star'

write.table(output,
    file=file.path('models', paste0(directory, '.dat')), 
    row.names=F, quote=F, sep='\t')

