#### Makes 9 models spanning the mass and alpha_MLT of each LEGACY star 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

set.seed(0)
options(scipen=99999)
#library(randomForest)

modes <- read.table(file.path('..', 'regression', 
        paste0('learn-inv'), paste0('tables-inv'), 
        'inversions_modes.dat'), 
    header=1)

star.dir <- file.path('..', 'regression', 'learn-inv', 'covs-inv', 'inversions')
for (star in list.files(star.dir)) {
    star.name <- strtoi(strsplit(star, '.dat')[[1]][1])
    if (is.na(star.name)) next 
    
    star.mode <- modes[modes$Name == star.name,]
    DF <- read.table(file.path(star.dir, star), header=1)
    
    M.mean <- mean(DF$M)
    M.sd <- sd(DF$M)
    alpha.mean <- mean(DF$alpha)
    alpha.sd <- sd(DF$alpha)
    
    lowhigh <- c('low', 'mean', 'high') 
    for (M in 1:3) { 
        for (R in 1:3) { 
            mass  <- c(    M.mean - M.sd,     star.mode$M,         M.mean + M.sd)[M] 
            alpha <- c(alpha.mean + alpha.sd, star.mode$alpha, alpha.mean - alpha.sd)[R] 
            command <- paste0("maybe_sub.sh -e -p 1 ./dispatch.sh -f",
                " -d models/", star.name, 
                " -n ", paste0(lowhigh[R], 'R', lowhigh[M], 'M'), 
                " -M ", mass, 
                " -Y ", star.mode$Y, 
                " -Z ", star.mode$Z, 
                " -a ", alpha, 
                " -c ", star.mode$age, 'd9', 
                " -D ", 0, 
                " -g ", 0, 
                " -o ", 0)
            print(command) 
            system(command) 
        } 
    } 
}
